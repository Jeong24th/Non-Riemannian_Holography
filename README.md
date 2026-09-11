# Non-Riemannian Holography Verification

Reproducibility scripts supporting the manuscript **“Non-Riemannian Hair in Long-String Holography”** by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

The scripts check the Riemannian and non-Riemannian response formulas, the doubled-yet-gauged worldsheet reduction, the Gomis–Ooguri limit, covariant charges, radial branches, the linear Virasoro condition, and supersymmetry.

## Mathematica suite

The Mathematica sources share background and frame definitions and use the current
A/M notation for particular two-point kernels. They include frame variations and
logarithmic response coefficients.

`mathematica/` contains exact-symbolic checks for the Letter and Supplemental
Material, with coverage and limitations stated in the equation ledger. The suite is
implemented in Wolfram Language and was tested with Mathematica 13.2:

```bash
wolframscript -file mathematica/NRH00_RunAll.wl
```

runs all 270 checks (a few minutes) and exits nonzero on any failure.  Every file
is also provided as a double-clickable `.nb` notebook with identical content — download
the `mathematica/` folder, open `NRH00_RunAll.nb`, and use *Evaluation → Evaluate
Notebook*.  See `mathematica/README.md` for the file-by-file coverage table and method
notes.  `mathematica/EQUATION_LEDGER.md` walks through every numbered equation of the
Letter and the Supplemental Material in order and states the public coverage of each one
(or states that it is a definition, cited statement, or currently uncovered).  The execution record is in
`mathematica/REFERENCE_RUN.md`, and `mathematica/MANUSCRIPT_MAP.md` records the
current manuscript SHA-256 and LaTeX-label-to-equation-number mapping.

The equation map is synchronized with the current Letter (1)-(22) and SM (1)-(202).
The September 11, 2026 source pin includes the latest prose and citation edits;
all 226 displayed equations are unchanged from the September 10 snapshot.
The expanded general-source SM1 derivation and current real Killing-spinor basis
are not fully covered by the public suite. In particular, a passing reduced jet
system is not a construction of nonzero supercharges. See the ledger for scope.

The manuscript source is not included. Python scripts run their algebraic checks
without it and report skipped LaTeX string comparisons. Three legacy comparisons
(`verify_sm_nr_linearization.py`, `verify_sm_riemannian_falloff.py`, and
`verify_gamma2_action.py`) still target the pre-rewrite SM1 snapshot identified in
[REFERENCE_RUN.md](mathematica/REFERENCE_RUN.md); they fail against the current source because formulas and labels
were replaced. They must not be used to validate current SM1. The doubled-yet-gauged
worldsheet source comparison still passes. Keep the current private manuscript
outside the repository when running the historical algebra regressions below.

## Python environment

The strict verification environment is:

- Python 3.12
- SymPy 1.14.0

Install the pinned Python dependency with:

```bash
python -m pip install -r requirements-verification.txt
```

## Algebraic regression checks

```bash
python checks/verify_sm_nr_linearization.py
python checks/verify_sm_riemannian_falloff.py
python checks/verify_dyg_reduction.py --strict-pin
python checks/verify_lambda_limit_ws.py --strict-pin
python checks/verify_10d_killing_spinor.py
python checks/verify_hairy_killing_spinor.py most-general
python checks/verify_n2_mirror_killing_spinor.py
python checks/verify_brst_w1.py
python checks/verify_gamma2_action.py
```

The `most-general` mode checks the reduced one-sided jet system in its own
frame convention; it does not verify the current full real Majorana basis.
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

The manuscript Data Availability Statement should cite a tagged release or immutable commit of this repository. `MANIFEST.sha256` records SHA-256 hashes of the tracked payload bytes (after Git line-ending normalization).
