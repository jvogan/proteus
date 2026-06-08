.PHONY: test release-check

test:
	python3 -m unittest -v

# Pre-publish hygiene sweep: no tracked structures/maps/media/secrets, no large
# files, no obvious secret strings. Run before pushing public changes.
release-check:
	@echo "Checking for tracked data/media/secret files..."
	@! git ls-files | grep -iE '\.(pdb|cif|mmcif|bcif|ent|map|mrc|mrcs|ccp4|mp4|mov|gif|webm|fasta|fa|fastq|pt|pth|npy|npz|pem|key)$$' | grep -vE '^tests/fixtures/' \
		|| { echo "ERROR: forbidden files tracked (see above)"; exit 1; }
	@echo "Checking for large files (>25MB)..."
	@! find . -path ./.git -prune -o -type f -size +25M -print | grep . \
		|| { echo "ERROR: large files present (see above)"; exit 1; }
	@echo "Scanning for secret-like strings..."
	@! git grep -nIE '(AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{36}|sk-[A-Za-z0-9]{20,}|xox[baprs]-|BEGIN (RSA|OPENSSH|EC|PGP) PRIVATE KEY)' -- . ':(exclude)Makefile' \
		|| { echo "ERROR: possible secrets (see above)"; exit 1; }
	@echo "release-check passed."
