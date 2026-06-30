# Structure File Formats - Proteus Skill

Use this reference when choosing file formats or when a parser/tool rejects a
structure file.

## Format Choices

| Format | Use for | Notes |
|---|---|---|
| `.cif` / `.mmcif` | Default for modern structures | Preferred by RCSB/PDBx; handles large structures better than legacy PDB |
| `.pdb` | Legacy compatibility | Easy to inspect; can fail for large structures or newer ID schemes |
| `.bcif` | Compact binary transfer | Good for web/model delivery, not ideal for simple stdlib parsing |
| `.sdf` | Small-molecule ligands | Use RDKit/Open Babel when chemical perception matters |
| `.mol2` | Small molecules with atom types/charges | Often useful for docking prep; parser support varies |
| `.mrc` / `.map` | Cryo-EM density maps | Use ChimeraX for visualization and map fitting |

## Practical Defaults

- Download experimental structures as mmCIF unless a downstream tool specifically
  needs legacy PDB.
- Use `scripts/structure_info.py` for quick PDB/mmCIF inspection.
- Use PyMOL or ChimeraX for robust parsing before making biological claims.
- Use ChimeraX for cryo-EM maps; PyMOL is not the right primary map viewer.
- Treat ligand files as chemistry data, not just coordinates. Use RDKit or Open
  Babel when bond orders, charges, protonation, or conformers matter.

## Gotchas

1. Legacy PDB format has fixed-width fields and limited chain/atom naming.
2. mmCIF is more verbose but is the safer default for modern PDB entries.
3. Biological assembly coordinates can differ from asymmetric unit coordinates.
4. AlphaFold stores pLDDT in the B-factor field for PDB/mmCIF outputs.
5. SDF/MOL2 ligand files may encode protonation and bond order differently from
   what appears in a protein structure file.
