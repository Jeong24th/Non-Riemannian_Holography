# Verification execution record

## Current documentation update: 2026-09-30

The current equation map and coverage audit target manuscript SHA-256
`57960AAD42F2A9150C0958D1C253F08B0939B76C0E646EE250A7D824D6319B68`.
The Letter has (1)-(21), with 23 numbered displays; the Supplemental Material
has (SM1)-(SM129), with 130 displays. The total is 153. Current SM3 contains
six subsections and 49 displays.

This update revises the documentation, equation mapping, comments and printed
check descriptions. The executable Wolfram algebra is unchanged. No Wolfram
runtime was available in the update environment, so **no new Wolfram execution
is reported**. The 336/336 result below belongs to the 2026-09-17 run and its
source, not to all formulas in the current manuscript.

Affected notebooks were serialized from their canonical Wolfram inputs and checked
for input parity and complete, balanced cells. The toolbox notebook was also
repaired: four function definitions previously split at internal blank lines now
remain in complete Input cells. Its Wolfram source and the runner are unchanged;
no notebook frontend or kernel execution is claimed for this update.

The compact SM3 replaces many expanded displays; their existing calculations
remain as auxiliary regression checks. Coverage is reassessed by content, not
inherited solely from matching labels. The [equation ledger](EQUATION_LEDGER.md)
and [revision notes](SM_REVISION_NOTES.md) identify the current correspondence.

The general off-shell Codazzi identity, universal tensor-Box variation identities,
exact radial metric variation and compact metric Hessian are not fully verified
by the old saddle/component checks. The fixed-flux K3 solution and interior H_R
matching are new to the manuscript; the old conditional normalization tests do
not verify the specified-vacuum hair-correlator derivation. The Einstein-tensor
assembly in the old constraint calculation assumes its defining decomposition.
New compact integrated coefficients have supporting differential-hierarchy checks,
but their entire rewritten presentation is not independently re-executed here.

## Recorded Mathematica run: 2026-09-17

Executed against manuscript SHA-256
`09B04C8BCF63EBDB4A879E1C893D28F43455756535340161D1138123B55DDA84`
(Letter (1)-(22), SM (1)-(205) of that source; 231 numbered displays).

- Environment: Mathematica 13.2.1 for Microsoft Windows (64-bit), January 27,
  2023; `Windows-x86-64`.
- Command: `wolframscript -file mathematica/NRH00_RunAll.wl`.
- Command-line result: **336/336**, exit 0; 303.6 seconds.
- Section counts: NRH02 49, NRH03 13, NRH04 20, NRH05 57, NRH06 67, NRH07 29,
  NRH08 32, NRH09 55, NRH10 14. Each section file also passed independently.
- The `.nb` notebooks were generated from the `.wl` sources by a kernel script
  with a Title cell and one Input cell per blank-line-separated block, holding
  the same input expressions.

The heavy steps were direct linearization of the two general saddle curvatures
through order e^(-2y/l), taking about 35 and 47 seconds in NRH05, and the ordered
cutoff bilinear with computed momenta, taking about 90 seconds in NRH06.

## Scope of the recorded calculations

The checks use symbolic chiral functions and NR hair where their inputs are
general; individual tests also include constant-profile, vacuum, one-sided and
other explicitly restricted sectors. Some predicates test definitions,
normalizations or algebraic consequences rather than a full physical claim.
The current ledger makes these distinctions explicit.

The near-boundary hierarchy retains stress responses, H_s and the constant mode
as free data. The two-point calculation uses the ordered second variation and
source linear response with two independent bulk solutions. Counterterm
cancellation is tested modulo tangential total derivatives after imposing the
stress constraints in both slots; omitting the counterterm is a negative control.
The particular nonlocal kernels follow from the Lorentzian inverse-derivative
prescription. This calculation does not determine general interior responses,
all finite local terms, or complete Green functions.

The full universal curvature-variation formulas were not separately checked:
linearized equations were obtained by direct curvature differentiation. The
reduced one-sided Killing-spinor jet system was checked in its frame convention;
it did not verify all ten-dimensional fermionic equations in the current real
Majorana basis or construct integrable fermionic charges. Worldsheet and
classical boundary-action checks likewise retain their stated sector restrictions.

## Earlier documentation update: 2026-09-19

The 2026-09-19 equation map targeted source
`EFB2D75953E9306C92165B532E855C590CA1FF8929FC61BF8EC8531BC522CD94`:
Letter (1)-(22), SM (1)-(209), 235 numbered displays. It followed prose changes,
the A notation for generalized-metric momentum, the fraktur notation for the
worldsheet auxiliary potential, and four added counterterm/Weyl-anomaly displays.
Those added displays had no public check at that time. No Wolfram re-execution
was recorded for that documentation update.

This is historical numbering. The 2026-09-30 manuscript has additional changes,
including replacement of the long SM3, so the old map and its coverage counts
must not be used as current coverage.

## Earlier execution and Python regressions

- 2026-09-10: the previous NRH00-NRH08 suite passed 270 checks in 253.7 seconds
  against source
  `FCBE00570616A89820C3F969E8E97BA74454B8555DDA5EBC25ED56C268C5BE6F`.
  That suite was superseded by the reorganized files used on 2026-09-17.
- Python scripts under `checks/` and `evidence/` are historical regressions.
  The three LaTeX comparisons in `verify_sm_nr_linearization.py`,
  `verify_sm_riemannian_falloff.py`, and `verify_gamma2_action.py` target an
  earlier snapshot and fail against the current source; they are not counted
  as current manuscript verification. Without the manuscript, string comparisons
  are skipped while the separate algebraic predicates can run.
- No Python verification code is changed by this documentation update.
