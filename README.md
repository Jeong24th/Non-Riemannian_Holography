# Non-Riemannian Holography Verification

Reproducibility scripts supporting the manuscript **“Non-Riemannian Hair in Long-String Holography”** by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

The scripts check the Riemannian and non-Riemannian response formulas, the doubled-yet-gauged worldsheet reduction, the Gomis–Ooguri limit, covariant charges, radial branches, the linear Virasoro condition, and supersymmetry.

## Mathematica suite (current manuscript)

The Mathematica sources share background and frame definitions and use the current
A/M notation for particular two-point kernels. They include frame variations and
logarithmic response coefficients.

`mathematica/` contains exact-symbolic checks for the Letter and Supplemental
Material, with coverage and limitations stated in the equation ledger. The suite is
implemented in Mathematica (tested with 13.2; the free Wolfram Engine also works):

```bash
wolframscript -file mathematica/NRH00_RunAll.wl
```

runs all 270 checks (a few minutes) and exits nonzero on any failure.  Every file
is also provided as a double-clickable `.nb` notebook with identical content — download
the `mathematica/` folder, open `NRH00_RunAll.nb`, and use *Evaluation → Evaluate
Notebook*.  See `mathematica/README.md` for the file-by-file coverage table and method
notes.  `mathematica/EQUATION_LEDGER.md` walks through every numbered equation of the
Letter and the Supplemental Material in order and names the check that verifies each one
(or states that it is a definition).  The neutral reference-execution record is in
`mathematica/REFERENCE_RUN.md`, and `mathematica/MANUSCRIPT_MAP.md` records the
current manuscript SHA-256 and LaTeX-label-to-equation-number mapping.

The manuscript source is not included in this software archive.  When
`NR_Holography.tex` is absent, the core scripts run their algebraic checks and
print that the LaTeX string comparison was skipped.  To check the displayed
formulas as well, place the manuscript source at the repository root.  The
source checks run inside the core Python scripts listed below.

## Python environment

The strict verification environment is:

- Python 3.12
- SymPy 1.14.0

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

The Riemannian falloff check also runs
`checks/verify_exact_projected_fluctuations.py`.  The charge calculations can be
run separately:

```bash
python checks/dft_asymptotic_charge.py
python checks/dft_covariant_phase_space.py
python checks/dft_translation_charge.py
python checks/dft_zero_mode_symplectic.py
```

Additional symbolic domain checks and negative controls are under `evidence/`.
The numerical `nr_corr_core.mjs`, sampled `go_ceff_contract_check.py`, and
`framing_contract.py` were retired in the 2026-09-06 audit.
The first used obsolete fluctuation and stress normalizations; the last checked
a fixed proposal without reading the manuscript. Their relevant mathematical coverage is supplied by
the current symbolic checks and the Mathematica suite.

## Scope

This public archive contains reproducibility software only. Internal companion-paper notes, review deliberations, and unpublished working documents are intentionally excluded.

## Versioning

The manuscript Data Availability Statement should cite a tagged release or immutable commit of this repository. `MANIFEST.sha256` records the current payload's Git-blob hashes.
