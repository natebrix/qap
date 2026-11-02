# Handoff validation — 2 November 2026

Built against cd43c7645d8af01a2ec22aa0ae9bba9fafea7022 using the SciPy OpenBLAS link configuration documented in README.md.

- `verify_witness.py`: exact PSD spectrum, all atomic triangles, 2,880 permutation-cut checks, 2,000 mixed-support checks, and Boolean SOS identity passed.
- `validate_instances.py`: all six grids and integer threshold decompositions passed; hashes in grid_validation.json.
- Rebuilt probe: toy baseline certificate 8.9999998013203442; cut certificate 9.9999999999801314.
- `make test` with matching BLAS/LAPACK variables: 4,536 brute-force runs passed; five resume cases identical; all parallel checks passed; seven regression instances matched exactly.

The original wil50 scalar logs are retained, not rerun for this packaging update. No strengthened Nugent experiment or regional M2 implementation is included.
