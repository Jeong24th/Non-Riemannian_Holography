# Non-Riemannian Holography Verification

Reproducibility scripts supporting **Non-Riemannian Holography: Long Strings and Soft Hair** by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

This archive contains symbolic checks of the exact saddles and boundary data,
linearized dynamics, holographic renormalization, covariant charges, worldsheet
reduction, Killing-spinor sectors, and a classical boundary action. Coverage
varies by equation; definitions, cited results, partial checks, and untested
results are distinguished in the [equation ledger](mathematica/EQUATION_LEDGER.md).

## Manuscript version and verification status

The documentation and equation map target the 2026-09-30 manuscript source:

```text
SHA-256: 57960AAD42F2A9150C0958D1C253F08B0939B76C0E646EE250A7D824D6319B68
```

The Letter has equations (1)-(21), comprising 23 numbered displays, and the
Supplemental Material has (SM1)-(SM129), comprising 130 numbered displays:
153 in total. Current SM3 has six subsections and 49 displays. Its compact
presentation replaces a much longer derivation; the corresponding expanded
calculations remain useful regression checks in the archive.

The latest recorded Wolfram execution is **336/336 checks passed on 2026-09-17**
with Mathematica 13.2.1, against an earlier manuscript. The 2026-09-30 update
changes the mapping, coverage descriptions, comments, and printed check names;
it leaves the executable algebra unchanged. No Wolfram runtime was available
for a new execution. This is not a new 336/336 verification of the current
manuscript. See the [execution record](mathematica/REFERENCE_RUN.md) and
[revision notes](mathematica/SM_REVISION_NOTES.md).

In particular, the general off-shell Codazzi identity, exact radial metric
variation and compact metric Hessian, and the new fixed-flux Bessel solution
and interior response are not fully checked by the retained suite. Older
conditional hair-correlator normalization tests do not verify the new
specified-vacuum derivation. The current [equation ledger](mathematica/EQUATION_LEDGER.md)
and [manuscript map](mathematica/MANUSCRIPT_MAP.md) give the scope of each display.

## Mathematica suite

Run all section files from the repository root:

```bash
wolframscript -file mathematica/NRH00_RunAll.wl
```

A failed check returns a nonzero exit status. The recorded full run took about
five minutes. Each `.wl` file also has a `.nb` notebook with the same input
expressions. Keep the `mathematica/` folder together, open `NRH00_RunAll.nb`,
and choose **Evaluation > Evaluate Notebook**. Section files run independently.
[mathematica/README.md](mathematica/README.md) explains the file organization,
conventions, and methods.

The suite directly linearizes the bulk curvature and solves the coupled
near-boundary hierarchy. Its two-point calculations use the ordered second
variation and source linear response, with two independent bulk solutions and
the stated source normalization. Counterterm cancellation is checked modulo
total tangential derivatives with stress constraints imposed in both slots.
These calculations determine particular kernels and leave general interior
responses and contact terms unresolved. General-background checks use
arbitrary chiral functions and arbitrary NR hair; restricted sectors are
identified in the ledger.

## Python environment

The pinned environment for the historical Python checks is:

- Python 3.12
- SymPy 1.14.0

```bash
python -m pip install -r requirements-verification.txt
```

## Algebraic regression checks (historical)

The scripts under `checks/` and `evidence/` predate the Mathematica suite.
Their algebraic predicates can run without the manuscript source and report
skipped LaTeX string comparisons. Three scripts target a pre-rewrite snapshot:
`verify_sm_nr_linearization.py`, `verify_sm_riemannian_falloff.py`, and
`verify_gamma2_action.py`. Their LaTeX comparisons fail against the current
source and must not be used as current-manuscript verification. Current public
coverage, including its gaps, is recorded in the Mathematica ledger. The
manuscript source is not included in this repository.

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
The Riemannian falloff script also runs
`checks/verify_exact_projected_fluctuations.py`. Charge calculations can be run
separately:

```bash
python checks/dft_asymptotic_charge.py
python checks/dft_covariant_phase_space.py
python checks/dft_translation_charge.py
python checks/dft_zero_mode_symplectic.py
```

Additional symbolic domain checks and negative controls are under `evidence/`.

## Scope and versioning

This public archive contains reproducibility software. Internal companion-paper
notes, review deliberations, and unpublished working documents are excluded.
The manuscript Data Availability Statement should cite a tagged release or an
immutable commit. `MANIFEST.sha256` records hashes of the tracked payload bytes
after Git line-ending normalization.
