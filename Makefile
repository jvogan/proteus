.PHONY: test release-check

test:
	python3 -m unittest -v

# Pre-publish hygiene sweep: no tracked structures/maps/media/secrets, no large
# files, no obvious secret strings. Run before pushing public changes.
release-check:
	@echo "Checking for tracked data/media/secret files..."
	@! git ls-files | grep -iE '\.(pdb|cif|mmcif|bcif|ent|map|mrc|mrcs|ccp4|mp4|mov|gif|webm|fasta|fa|fastq|pt|pth|npy|npz|pem|key)$$' | grep -vE '^tests/fixtures/' \
		|| { echo "ERROR: forbidden files tracked (see above)"; exit 1; }
	@echo "Checking for generated/private directories accidentally tracked..."
	@! git ls-files | grep -E '^(\.proteus-cache/|banners/|kras_g12c_dossier/|[^/]+_dossier/|proteus_batch_out/|proteus_report_out/|design_run_out/|dock_prep_out/|rosetta_score_out/|structures/|ligands/)' \
		|| { echo "ERROR: generated/private directories tracked (see above)"; exit 1; }
	@echo "Checking for local machine paths in tracked/untracked files..."
	@! git grep --untracked -nIE '(/Users/|/home/[^/ ]+)' -- . ':(exclude)Makefile' \
		|| { echo "ERROR: local machine paths found in tracked files (see above)"; exit 1; }
	@echo "Checking for large files (>25MB)..."
	@! find . -path ./.git -prune -o -type f -size +25M -print | grep . \
		|| { echo "ERROR: large files present (see above)"; exit 1; }
	@echo "Scanning for secret-like strings..."
	@! git grep --untracked -nIE '(AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{36}|sk-[A-Za-z0-9]{20,}|xox[baprs]-|BEGIN (RSA|OPENSSH|EC|PGP) PRIVATE KEY)' -- . ':(exclude)Makefile' \
		|| { echo "ERROR: possible secrets (see above)"; exit 1; }
	@echo "release-check passed."
