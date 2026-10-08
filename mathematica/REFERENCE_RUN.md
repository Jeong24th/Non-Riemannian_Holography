# Mathematica execution — October 8, 2026

**451/451 checks passed in 14 independent programs**, using Mathematica 13.2.1 on Windows (64-bit).
Pinned manuscript SHA-256: `2E6A7A8A20D581043D6457865662E76F7D46433AF29847393F1DBA4BF753F822`.

Each program was copied into its own isolated directory and executed by a fresh `wolframscript -file` process. No sibling source, package, or manuscript input was available there. All exit codes were zero; no Wolfram diagnostic messages or failed assertions appeared. Source hashes were rechecked after execution. Each matching notebook's complete embedded program was parsed and compared with its WL source, normalizing line endings and the final newline.

The [machine-readable receipt](RUN_RESULTS.json) records execution times, counts, exact executed hashes and LF-normalized hashes. The [notebook validation log](runs/2026-10-08/notebooks.txt) covers all 14 pairs. Notebook payload equivalence was checked in the kernel; this does not claim a separate front-end execution of every notebook.

| Program | Checks | Seconds | Log |
|---|---:|---:|---|
| [01_ExactBackgrounds_EDFE.wl](01_ExactBackgrounds_EDFE.wl) | 49/49 | 39.4 | [output](runs/2026-10-08/01_ExactBackgrounds_EDFE.txt) |
| [02_RenormalizedOnShellAction.wl](02_RenormalizedOnShellAction.wl) | 15/15 | 16.8 | [output](runs/2026-10-08/02_RenormalizedOnShellAction.txt) |
| [03_RadialSolutionDerivation.wl](03_RadialSolutionDerivation.wl) | 20/20 | 18.5 | [output](runs/2026-10-08/03_RadialSolutionDerivation.txt) |
| [04_GeneralBackground_LEDFE.wl](04_GeneralBackground_LEDFE.wl) | 61/61 | 76.3 | [output](runs/2026-10-08/04_GeneralBackground_LEDFE.txt) |
| [05_QuadraticRenormalization_TwoPoint.wl](05_QuadraticRenormalization_TwoPoint.wl) | 67/67 | 76.9 | [output](runs/2026-10-08/05_QuadraticRenormalization_TwoPoint.txt) |
| [06_NoetherCharges.wl](06_NoetherCharges.wl) | 31/31 | 3.1 | [output](runs/2026-10-08/06_NoetherCharges.txt) |
| [07_Worldsheet.wl](07_Worldsheet.wl) | 49/49 | 3.5 | [output](runs/2026-10-08/07_Worldsheet.txt) |
| [08_Uplift_Killing.wl](08_Uplift_Killing.wl) | 55/55 | 61.7 | [output](runs/2026-10-08/08_Uplift_Killing.txt) |
| [09_BoundaryCandidate.wl](09_BoundaryCandidate.wl) | 14/14 | 2.7 | [output](runs/2026-10-08/09_BoundaryCandidate.txt) |
| [10_AsymptoticSymmetries.wl](10_AsymptoticSymmetries.wl) | 7/7 | 3.6 | [output](runs/2026-10-08/10_AsymptoticSymmetries.txt) |
| [11_ZeroL_ExactLEDFE.wl](11_ZeroL_ExactLEDFE.wl) | 5/5 | 7.2 | [output](runs/2026-10-08/11_ZeroL_ExactLEDFE.txt) |
| [12_RiemannianVacuum_ExactLEDFE.wl](12_RiemannianVacuum_ExactLEDFE.wl) | 8/8 | 10.8 | [output](runs/2026-10-08/12_RiemannianVacuum_ExactLEDFE.txt) |
| [13_TwoPointFunctions.wl](13_TwoPointFunctions.wl) | 19/19 | 7.0 | [output](runs/2026-10-08/13_TwoPointFunctions.txt) |
| [14_SM29_IntrinsicCounterterm.wl](14_SM29_IntrinsicCounterterm.wl) | 51/51 | 4.2 | [output](runs/2026-10-08/14_SM29_IntrinsicCounterterm.txt) |

This is a fresh October 8 run. The retained [October 6 logs](runs/2026-10-06/) are historical evidence for their earlier sources, not the current execution receipt. The totals include retained ancillary calculations, including seven checks of a classical winding-energy passage deleted from the current manuscript.

File 14 derives the intrinsic SM29 quadratic expansion. It verifies the NR-vacuum divergent and finite match, the R-vacuum divergent match with an explicit finite local contact difference, and a nonzero divergent obstruction for generic constant L_pm on both branches. The program prints the intermediate curvatures, candidate densities, finite difference and obstruction. Passing its obstruction checks establishes the stated failure away from the vacua, not a universal nonlinear counterterm.

Exact L_pm=0 and order-z general-background checks remain distinct. These tests do not settle all interior conditions, homogeneous responses, zero modes or quotient prescriptions. The [equation ledger](EQUATION_LEDGER.md) records the individual scopes and uncovered identities.
