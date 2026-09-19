# Verification execution record

## Recorded Mathematica run of the reorganized suite

Executed 2026-09-17 against manuscript SHA-256
`09B04C8BCF63EBDB4A879E1C893D28F43455756535340161D1138123B55DDA84`
(Letter (1)-(22), SM (1)-(205) of that source; 231 numbered displays).

- Environment: 13.2.1 for Microsoft Windows (64-bit) (January 27, 2023); `Windows-x86-64`.
- Command: `wolframscript -file mathematica/NRH00_RunAll.wl`.
- Command-line result: **336/336**, exit 0; 303.6 seconds.
- Section counts: NRH02 49, NRH03 13, NRH04 20, NRH05 57, NRH06 67, NRH07 29,
  NRH08 32, NRH09 55, NRH10 14. Every section file also passes when run on its own.
- The `.nb` notebooks were generated from the `.wl` sources by a kernel script (Title cell
  followed by one Input cell per blank-line-separated block) and hold the same input
  expressions.

The heavy steps are the direct linearization of the exact EDFE on the two general saddles
through order e^{-2y/l} (NRH05, about 35 s and 47 s) and the assembly of the ordered cutoff
bilinear with the computed momenta (NRH06, about 90 s in total).

## What the run establishes

All checks are exact symbolic identities in arbitrary chiral functions L±(x±) and, on the
non-Riemannian branch, an arbitrary hair function W₁(x⁺,x⁻), unless a check name states a
restriction (constant L sectors, the one-sided L₋ = 0 family, the common vacuum). The
near-boundary solutions are solved from the coupled hierarchy with free responses R±, H_s and
zero mode c_s; the two-point functions are obtained from the ordered second variation with two
independent bulk solutions; counterterm cancellation is tested modulo total tangential
derivatives with the stress constraints imposed in both slots, with a negative control.

Not evaluated: the undetermined interior responses and the hair self-response, the finite
local terms induced by the cutoff subtraction, complete Green functions, the general variation
formulas of SM 3.5 (the linearized EDFE are computed directly), and the full ten-dimensional
fermionic equations in the current real Majorana basis (the reduced one-sided jet system is
verified in the printed frame). See EQUATION_LEDGER.md for per-equation status.

## Manuscript revisions after the recorded run

The recorded run used source `09B04C8BCF63EBDB4A879E1C893D28F43455756535340161D1138123B55DDA84` (2026-09-17). The current
manuscript source is `EFB2D75953E9306C92165B532E855C590CA1FF8929FC61BF8EC8531BC522CD94` (2026-09-19). The revisions in between:

- Prose-only edits (2026-09-18 wording changes and the 2026-09-19 proofreading pass); no
  displayed equation body changed apart from the two notation renames listed next.
- Notation: the generalized-metric momentum is written A instead of 𝒜 in SM 3, and the worldsheet
  auxiliary gauge potential is written with a fraktur a instead of 𝒜 in SM 4 (2026-09-18). The
  checks are symbolic and print manuscript labels, so they are unaffected.
- Four new displays (2026-09-18): SM 105–107 (intrinsic boundary-curvature counterterm candidate,
  end of SM 3.7) and SM 124 (linearized Weyl anomaly, SM 3.8). Equation numbers shift accordingly:
  old SM 105–120 become SM 108–123 and old SM 121–205 become SM 125–209; SM 1–104 and the Letter
  numbers are unchanged. The four new displays have no public check yet (**N**, or **D** for the
  definitions); every other ledger entry keeps its recorded coverage.
- The suite was not re-executed for this documentation update; the check labels printed by the
  recorded run remain valid for the renumbered ledger.

## Historical runs

- 2026-09-10: the previous suite (NRH00–NRH08, 270 checks, 253.7 s) against source
  `FCBE00570616A89820C3F969E8E97BA74454B8555DDA5EBC25ED56C268C5BE6F`. That suite is superseded by
  the present files; its coverage is contained in the new files (the response checks now follow
  the SM 3.4–3.8 route with general sources).
- The Python scripts under `checks/` and `evidence/` are historical regressions; the three
  LaTeX-contract comparisons that target the pre-rewrite SM 1 snapshot fail against the current
  source and are not counted. No Python verification code was changed in this release.
