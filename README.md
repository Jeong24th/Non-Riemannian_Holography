# Non-Riemannian Holography Verification

Reproducibility scripts supporting the manuscript **“Long Strings and Non-Riemannian Hair”** by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

This archive contains the symbolic and analytic checks used to verify the Riemannian and non-Riemannian holographic response, doubled-yet-gauged worldsheet reduction, effective Gomis–Ooguri limit, covariant charges, radial branches, BRST closure, and supersymmetry statements.

## Environment

The strict verification environment is:

- Python 3.12
- SymPy 1.14.0
- Node.js for the standalone `.mjs` check

Install the pinned Python dependency with:

```bash
python -m pip install -r requirements-verification.txt
```

## Core checks

```bash
python checks/verify_sm_nr_linearization.py
python checks/verify_sm_riemannian_falloff.py
python checks/verify_dyg_reduction.py --strict-pin
python checks/verify_lambda_limit_ws.py --strict-pin
python checks/verify_10d_killing_spinor.py
python checks/verify_hairy_killing_spinor.py
python checks/verify_n2_mirror_killing_spinor.py
python checks/verify_brst_w1.py
python checks/verify_gamma2_action.py
```

Additional independent/domain/negative-control checks are under `evidence/`.

## Mathematica suite

`mathematica/` contains an independent exact-symbolic verification of the displayed
equations of the Letter and the Supplemental Material, organized in paper order and
implemented in Mathematica (tested with 13.2; the free Wolfram Engine also works):

```bash
wolframscript -file mathematica/NRH00_RunAll.wl
```

runs all 169 checks (about ten minutes) and exits nonzero on any failure.  Every file
is also provided as a double-clickable `.nb` notebook with identical content — download
the `mathematica/` folder, open `NRH00_RunAll.nb`, and use *Evaluation → Evaluate
Notebook*.  See `mathematica/README.md` for the file-by-file coverage table and method
notes.

Some contract-aware guards also compare verified expressions against `NR_Holography.tex`. For those checks, place the manuscript source at the repository root. Their independent algebraic checks remain documented in the scripts.

## Scope

This public archive contains reproducibility software only. Internal companion-paper notes, review deliberations, and unpublished working documents are intentionally excluded.

## Versioning

The manuscript Data Availability Statement should cite a tagged release or immutable commit of this repository. File hashes for the initial archive are recorded in `MANIFEST.sha256`.
