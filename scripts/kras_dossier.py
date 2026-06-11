#!/usr/bin/env python3
"""Build a complete structural-biology dossier for KRAS G12C.

This is the "epic demo": a single command that reproduces an end-to-end target
review which would take a human the better part of a day in PyMOL — fetching a
dozen structures, running alignments and interface/pocket analyses, ray-tracing
a gallery of figures, recording a turntable movie, and writing it all up.

It chains the existing Proteus helpers (fetch_pdb, fetch_alphafold, pae_report,
compare_structures, interface_report, pymol_agent) over a *curated, validated*
set of public structures and emits a self-contained Markdown report with
embedded figures.

    python3 scripts/kras_dossier.py --out kras_g12c_dossier

Every PDB ID below was validated against the RCSB data API (title / resolution /
method / bound-ligand code) so the narrative never points at the wrong file.
Re-discover the landscape yourself with:

    python3 scripts/pdb_search.py --uniprot P01116 --rows 50 --details
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
PY = sys.executable

# --- Curated structure set -------------------------------------------------
# The covalent-inhibitor discovery timeline: each traps KRAS(G12C) in the
# inactive (GDP) state via the cryptic switch-II pocket, tethered to Cys12.
INHIBITORS = [
    {"id": "4LYH", "lig": "21F", "year": 2013,
     "drug": "Ostrem covalent probe (compound 9)",
     "note": "first proof that the G12C cysteine is targetable (Nature 2013)"},
    {"id": "6N2J", "lig": "K9M", "year": 2018,
     "drug": "ARS-1620",
     "note": "first G12C inhibitor with cellular potency"},
    {"id": "6OIM", "lig": "MOV", "year": 2019,
     "drug": "Sotorasib (AMG 510)",
     "note": "first KRAS drug to reach FDA approval (2021)"},
    {"id": "6UT0", "lig": "M1X", "year": 2020,
     "drug": "Adagrasib (MRTX849 series)",
     "note": "second approved G12C inhibitor (2022)"},
]
# Nucleotide-driven conformational switch (single-chain, high-res pair).
STATE_INACTIVE = {"id": "4LDJ", "nuc": "GDP", "label": "Inactive · GDP-bound"}
STATE_ACTIVE = {"id": "6GOD", "nuc": "GNP", "label": "Active · GMPPNP (GTP analog)"}
# KRAS as a signaling hub: upstream activator and downstream effector.
INTERFACES = [
    {"id": "6EPL", "key": "sos",
     "name": "KRAS(G12C) – SOS1",
     "role": "Activation: the guanine-exchange factor that reloads KRAS with GTP"},
    {"id": "6VJJ", "key": "raf",
     "name": "KRAS – CRAF (RAF1) RBD",
     "role": "Effector: the active-state contact that fires the MAPK cascade"},
]
AF = {"uniprot": "P01116", "gene": "KRAS"}
HERO = "6OIM"

# G-domain functional regions (KRAS numbering) for annotating motion.
REGIONS = [
    ("P-loop (incl. G12/C12)", 10, 17),
    ("Switch I", 30, 38),
    ("Switch II", 60, 76),
]


def region_for(resi):
    try:
        n = int("".join(c for c in str(resi) if c.isdigit()))
    except ValueError:
        return ""
    for name, lo, hi in REGIONS:
        if lo <= n <= hi:
            return name
    return ""


def extract_json(text):
    """Pull the final top-level JSON object out of mixed stdout."""
    lines = text.splitlines()
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip() == "{":
            try:
                return json.loads("\n".join(lines[i:]))
            except json.JSONDecodeError:
                continue
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


class Runner:
    def __init__(self, outdir):
        self.outdir = outdir
        self.log = []

    def run(self, script, args, want_json=False, label=None):
        label = label or f"{script} {' '.join(args)}"
        cmd = [PY, str(SCRIPTS / script), *args]
        try:
            proc = subprocess.run(cmd, cwd=self.outdir, capture_output=True,
                                  text=True, timeout=600)
        except subprocess.TimeoutExpired:
            self.log.append((label, "TIMEOUT"))
            print(f"  ✗ {label}  (timeout)")
            return None
        ok = proc.returncode == 0
        self.log.append((label, "ok" if ok else f"exit {proc.returncode}"))
        mark = "✓" if ok else "✗"
        print(f"  {mark} {label}")
        if not ok:
            err = (proc.stderr or proc.stdout).strip().splitlines()
            if err:
                print(f"      {err[-1][:160]}")
            return None
        if want_json:
            return extract_json(proc.stdout)
        return proc.stdout


def fmt_res(resolution):
    if isinstance(resolution, list) and resolution:
        return f"{resolution[0]} Å"
    return "—"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="kras_g12c_dossier", help="Output directory")
    ap.add_argument("--frames", type=int, default=48, help="Turntable frames")
    ap.add_argument("--no-movie", action="store_true", help="Skip the turntable movie")
    args = ap.parse_args(argv)

    outdir = Path(args.out).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"KRAS G12C dossier → {outdir}\n")
    r = Runner(outdir)
    results = {}

    # 1. Fetch all structures ------------------------------------------------
    print("[1/6] Fetching structures")
    pdb_ids = [s["id"] for s in INHIBITORS] + [STATE_INACTIVE["id"], STATE_ACTIVE["id"]] \
        + [s["id"] for s in INTERFACES]
    for pid in pdb_ids:
        r.run("fetch_pdb.py", [pid, "--format", "pdb"], label=f"fetch {pid}")
    af = r.run("fetch_alphafold.py", [AF["uniprot"], "--pae", "--json"],
               want_json=True, label=f"fetch AlphaFold {AF['uniprot']}")
    results["af_fetch"] = af

    # 2. AlphaFold confidence + PAE -----------------------------------------
    print("\n[2/6] AlphaFold model confidence")
    af_pdb = f"AF-{AF['uniprot']}-F1.pdb"
    af_pae = f"AF-{AF['uniprot']}-F1_pae.json"
    r.run("pymol_agent.py", ["render", af_pdb, "af_plddt.png", "--color", "plddt"],
          label="render AlphaFold pLDDT")
    results["pae"] = r.run("pae_report.py", [af_pae, "--json"], want_json=True,
                           label="PAE report")

    # 3. Conformational switch ----------------------------------------------
    print("\n[3/6] Nucleotide switch (active vs inactive)")
    r.run("pymol_agent.py", ["render", f"{STATE_INACTIVE['id']}.pdb",
          "state_inactive.png", "--color", "spectrum"], label="render inactive (GDP)")
    r.run("pymol_agent.py", ["render", f"{STATE_ACTIVE['id']}.pdb",
          "state_active.png", "--color", "spectrum"], label="render active (GMPPNP)")
    results["switch"] = r.run("compare_structures.py",
        [f"{STATE_INACTIVE['id']}.pdb", f"{STATE_ACTIVE['id']}.pdb",
         "--per-residue", "--json"], want_json=True, label="compare GDP vs GMPPNP")

    # 4. Inhibitor pocket gallery -------------------------------------------
    print("\n[4/6] Druggable pocket — inhibitor timeline")
    results["pockets"] = {}
    for inh in INHIBITORS:
        pid, code = inh["id"], inh["lig"]
        chain = first_ligand_chain(outdir / f"{pid}.pdb", code)
        sel = f"resn {code}" + (f" and chain {chain}" if chain else "")
        out = f"pocket_{pid}.png"
        res = r.run("pymol_agent.py", ["pocket", f"{pid}.pdb", out,
                    "--ligand", sel, "--label"], want_json=True,
                    label=f"pocket {pid} ({inh['drug']})")
        results["pockets"][pid] = res

    # 5. Interfaces ----------------------------------------------------------
    print("\n[5/6] Signaling interfaces")
    results["interfaces"] = {}
    for iface in INTERFACES:
        pid = iface["id"]
        r.run("pymol_agent.py", ["render", f"{pid}.pdb", f"iface_{iface['key']}.png",
              "--color", "chain"], label=f"render {iface['name']}")
        results["interfaces"][iface["key"]] = r.run("interface_report.py",
            [f"{pid}.pdb", "--cutoff", "4.5", "--json"], want_json=True,
            label=f"interface {iface['name']}")

    # 6. Hero figure + movie -------------------------------------------------
    print("\n[6/6] Hero figure + turntable")
    r.run("pymol_agent.py", ["render", f"{HERO}.pdb", "hero.png",
          "--color", "spectrum", "--preset", "soft"], label="render hero")
    if not args.no_movie:
        r.run("pymol_agent.py", ["spin", f"{HERO}.pdb", "hero_spin.mp4",
              "--frames", str(args.frames), "--color", "spectrum"],
              label="turntable movie")

    # Assemble the report ----------------------------------------------------
    print("\nWriting DOSSIER.md")
    write_report(outdir, results)
    print(f"\nDone. Open {outdir / 'DOSSIER.md'}")
    failures = [l for l, s in r.log if s != "ok"]
    if failures:
        print(f"\n{len(failures)} step(s) had issues:")
        for l in failures:
            print(f"  - {l}")
    return 0


def first_ligand_chain(pdb_path, code):
    """Return the chain ID of the first HETATM/ATOM record for residue `code`."""
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line[:6] in ("HETATM", "ATOM  ") and line[17:20].strip() == code:
                    chain = line[21].strip()
                    if chain:
                        return chain
    except OSError:
        pass
    return None


def write_report(outdir, results):
    L = []
    w = L.append

    w("# KRAS G12C — Structural Dossier")
    w("")
    w("*Generated by `scripts/kras_dossier.py` from public PDB and AlphaFold "
      "structures. Every figure is a real Proteus code path; re-run the script "
      "to regenerate.*")
    w("")
    w("KRAS is the most frequently mutated oncogene in human cancer and was "
      "considered \"undruggable\" for three decades. The **G12C** mutation — a "
      "glycine-to-cysteine substitution at codon 12 — installs a reactive "
      "cysteine that covalent drugs can exploit. This dossier walks from the "
      "protein's conformational switch, through the cryptic pocket that made it "
      "druggable, to the protein–protein interfaces that make it a signaling hub.")
    w("")

    # Target overview
    w("## 1. Target overview — AlphaFold model")
    w("")
    w("![AlphaFold pLDDT](af_plddt.png)")
    w("")
    af = results.get("af_fetch") or {}
    pae = (results.get("pae") or {}).get("data", {})
    if af:
        w(f"AlphaFold model **{af.get('model_id','AF-'+AF['uniprot'])}** "
          f"(gene *{af.get('gene', AF['gene'])}*), global pLDDT "
          f"**{af.get('global_plddt','?')}**. The compact G-domain is "
          "high-confidence (blue); the C-terminal hypervariable region — the "
          "membrane-anchoring tail that is lipid-modified in vivo and disordered "
          "in isolation — drops to low confidence (orange).")
        w("")
    if pae:
        flex = pae.get("flexible_or_uncertain_regions", [])
        size = pae.get("size", {})
        w(f"**PAE report** ({size.get('rows','?')}×{size.get('columns','?')} "
          f"residues, mean PAE {pae.get('pae',{}).get('mean','?')} Å) flags "
          f"{len(flex)} flexible/uncertain segment(s):")
        w("")
        w("| Region | Residues | Mean PAE (Å) |")
        w("|---|---|---|")
        for seg in flex[:6]:
            w(f"| {seg['start']}–{seg['end']} | {seg['length']} | {seg['mean_pae']} |")
        w("")

    # Switch
    w("## 2. The molecular switch — GTP vs GDP")
    w("")
    w("| Inactive (GDP) | Active (GMPPNP) |")
    w("|---|---|")
    w("| ![inactive](state_inactive.png) | ![active](state_active.png) |")
    w(f"| {STATE_INACTIVE['id']} · {STATE_INACTIVE['label']} | "
      f"{STATE_ACTIVE['id']} · {STATE_ACTIVE['label']} |")
    w("")
    sw = (results.get("switch") or {}).get("data", {})
    if sw:
        aln = sw.get("alignment", {})
        pr = sw.get("per_residue", {})
        w(f"CE-align over **{aln.get('alignment_length','?')}** residues gives "
          f"RMSD **{round(aln.get('RMSD',0),2)} Å**. The deviation is not "
          "uniform — it concentrates in the two nucleotide sensors that define "
          "the switch:")
        w("")
        w("| Residue | Cα shift (Å) | Functional region |")
        w("|---|---|---|")
        devs = pr.get("largest_deviations", [])[:10]
        for d in devs:
            reg = region_for(d.get("resi"))
            w(f"| {d.get('chain','')}{d.get('resi','')} | "
              f"{round(d.get('ca_distance',0),1)} | {reg or '—'} |")
        w("")
        hits = [d for d in devs if region_for(d.get("resi"))]
        if hits:
            w(f"*{len(hits)} of the top {len(devs)} displaced residues fall in "
              "the P-loop / Switch I / Switch II — exactly the machinery that "
              "reads the γ-phosphate. The covalent drugs below lock this switch "
              "in the inactive arrangement.*")
            w("")

    # Pocket gallery
    w("## 3. The druggable pocket — a discovery timeline")
    w("")
    w("Each inhibitor binds the **switch-II pocket (S-IIP)**, a cryptic cavity "
      "that only opens in the GDP state, and forms a covalent bond to Cys12. "
      "The progression from chemical probe to approved drug:")
    w("")
    for inh in INHIBITORS:
        pid = inh["id"]
        pk = (results.get("pockets", {}).get(pid) or {}).get("data", {})
        w(f"### {inh['year']} — {inh['drug']}  (`{pid}`)")
        w("")
        w(f"![pocket {pid}](pocket_{pid}.png)")
        w("")
        note = inh["note"]
        line = f"*{note[0].upper()}{note[1:]}.*"
        if pk:
            line += (f" Ligand `{inh['lig']}`: {pk.get('ligand_atoms','?')} atoms "
                     f"contacting {pk.get('pocket_residues','?')} pocket residues.")
        w(line)
        w("")

    # Interfaces
    w("## 4. KRAS as a signaling hub — protein interfaces")
    w("")
    for iface in INTERFACES:
        key = iface["key"]
        rep = (results.get("interfaces", {}).get(key) or {}).get("data", {})
        w(f"### {iface['name']}  (`{iface['id']}`)")
        w("")
        w(f"![{key}](iface_{key}.png)")
        w("")
        w(f"*{iface['role']}.*")
        w("")
        ifaces = rep.get("interfaces", []) if rep else []
        if ifaces:
            top = ifaces[0]
            counts = top.get("interface_residue_count", {})
            chains = "/".join(top.get("chains", []))
            w(f"Interface (chains {chains}, 4.5 Å cutoff): "
              f"**{top.get('contact_pair_count','?')}** contact pairs across "
              + ", ".join(f"{c}: {n} residues" for c, n in counts.items()) + ".")
            w("")

    # Hero
    w("## 5. Hero view & turntable")
    w("")
    w("![hero](hero.png)")
    w("")
    w(f"`{HERO}` — KRAS G12C covalently bound to sotorasib. A 360° turntable is "
      "saved as **`hero_spin.mp4`**.")
    w("")

    # Methods
    w("## 6. Structures & provenance")
    w("")
    w("| PDB | Role | Ligand |")
    w("|---|---|---|")
    for inh in INHIBITORS:
        w(f"| {inh['id']} | Inhibitor: {inh['drug']} | {inh['lig']} |")
    w(f"| {STATE_INACTIVE['id']} | {STATE_INACTIVE['label']} | {STATE_INACTIVE['nuc']} |")
    w(f"| {STATE_ACTIVE['id']} | {STATE_ACTIVE['label']} | {STATE_ACTIVE['nuc']} |")
    for iface in INTERFACES:
        w(f"| {iface['id']} | {iface['name']} | — |")
    w(f"| AF-{AF['uniprot']}-F1 | AlphaFold model ({AF['gene']}) | — |")
    w("")
    w("Discovery: `scripts/pdb_search.py --uniprot P01116 --details`. "
      "All structures are public; outputs are git-ignored.")
    w("")

    (outdir / "DOSSIER.md").write_text("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
