# Non-Riemannian Holography Verification

Reproducibility scripts supporting the manuscript **“Non-Riemannian Holography: Long Strings and Soft Hair”** by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

The scripts check the two exact saddles and their boundary data, the linearized dynamics and
holographic renormalization (near-boundary solutions, counterterms, one- and two-point
functions), the covariant charges, the doubled-yet-gauged worldsheet reduction and the
Gomis–Ooguri limit, the ten-dimensional uplift and Killing spinors, and the boundary
candidate action.

## Mathematica suite

`mathematica/` contains exact-symbolic checks organized in the order of the manuscript: one
file for the Letter and one file per Supplemental Material section (SM 3 is split into three
files). The suite is implemented in Wolfram Language and was tested with Mathematica 13.2:

```bash
wolframscript -file mathematica/NRH00_RunAll.wl
```

runs all 336 checks (about five minutes) and exits nonzero on any failure. Every file is
also provided as a double-clickable `.nb` notebook with identical content — download the
`mathematica/` folder, open `NRH00_RunAll.nb`, and use *Evaluation → Evaluate Notebook*.
See `mathematica/README.md` for the file-by-file coverage table, the conventions and the
method notes. `mathematica/EQUATION_LEDGER.md` walks through every numbered equation of
the Letter (1)–(22) and the Supplemental Material (1)–(209) in order and states the public
coverage of each one (or states that it is a definition, a cited statement, or currently
uncovered). Four displays added to the manuscript on 2026-09-18 (SM 105–107, the
boundary-curvature counterterm candidate, and SM 124, the linearized Weyl anomaly) are not
yet covered by the suite and are marked accordingly. The execution record is in
`mathematica/REFERENCE_RUN.md`, and
`mathematica/MANUSCRIPT_MAP.md` records the current manuscript SHA-256 and the
LaTeX-label-to-equation-number mapping.

The computations follow the methods stated in the manuscript. In particular, the two-point
functions are computed from the ordered second variation of the on-shell action with two
independent bulk solutions inserted (SM 3.8), not by differentiating one-point functions; the
near-boundary solutions are obtained from the coupled hierarchy of the linearized EDFE (SM
3.4); and the counterterm cancellation is tested modulo total tangential derivatives with the
stress constraints imposed. Backgrounds carry arbitrary chiral functions L±(x±) and, on the
non-Riemannian branch, an arbitrary hair function W₁(x⁺,x⁻); no sample data are used.

## Python environment

The strict verification environment for the historical Python checks is:

- Python 3.12
- SymPy 1.14.0

Install the pinned Python dependency with:

```bash
python -m pip install -r requirements-verification.txt
```

## Algebraic regression checks (historical)

The Python scripts under `checks/` and `evidence/` predate the current Mathematica suite and
are retained as historical regressions. Their algebraic predicates run without the manuscript
source and report skipped LaTeX string comparisons; three of them
(`verify_sm_nr_linearization.py`, `verify_sm_riemannian_falloff.py`, and
`verify_gamma2_action.py`) target a pre-rewrite SM snapshot and fail against the current
source because formulas and labels were replaced. They must not be used to validate the
current manuscript; the corresponding coverage is supplied by the Mathematica suite. Keep the
private manuscript outside the repository when running them.

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

The `most-general` mode checks the reduced one-sided jet system in its own frame convention;
it does not verify the current full real Majorana basis. The Riemannian falloff check also
runs `checks/verify_exact_projected_fluctuations.py`. The charge calculations can be run
separately:

```bash
python checks/dft_asymptotic_charge.py
python checks/dft_covariant_phase_space.py
python checks/dft_translation_charge.py
python checks/dft_zero_mode_symplectic.py
```

Additional symbolic domain checks and negative controls are under `evidence/`.

## Scope

This public archive contains reproducibility software only. Internal companion-paper notes,
review deliberations, and unpublished working documents are intentionally excluded.

## Versioning

The manuscript Data Availability Statement should cite a tagged release or immutable commit of
this repository. `MANIFEST.sha256` records SHA-256 hashes of the tracked payload bytes (after
Git line-ending normalization).
