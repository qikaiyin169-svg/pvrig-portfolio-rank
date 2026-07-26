# Data provenance and redistribution

This repository contains the generated candidate sequences and the fixed
stage-one result snapshot needed to audit the final ranking.

It does **not** vendor the organizer's 48-positive reference library or the PDB
coordinate files. `scripts/fetch_inputs.py` retrieves those inputs from their
authoritative upstream locations and verifies their pinned versions:

- `clickmab-bio/ab-data-validator`, commit
  `97df17aa09bc576a861cf0d8242de97af379fd80`;
- RCSB PDB entries `8X6B` and `9E6Y`, verified by SHA-256.

The code is MIT licensed. Third-party structures, reference sequences, model
weights and patent-derived sequence disclosures retain their own terms. The
MIT license does not grant patent rights or relicense third-party data.
