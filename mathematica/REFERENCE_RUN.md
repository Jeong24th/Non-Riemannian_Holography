# Fresh Mathematica execution — October 6, 2026

**377/377 checks passed in 13 independent programs**, using Mathematica 13.2.1 on Windows (64-bit).
Pinned manuscript SHA-256: `F8AC68BB80B6C3D572160FD83705FE7506FCF277B00FB95021D39C2AF1E8E239`.

Each program was copied into an otherwise empty directory and executed by a fresh `wolframscript -file` process. No sibling tools, source files or packages were available in that directory. All process exit codes were zero. Source hashes were rechecked after execution. Each matching notebook's embedded program was parsed and compared with its WL source (normalizing line endings and the final newline).

The complete [machine-readable receipt](RUN_RESULTS.json) contains execution times, per-file totals, exact executed hashes and LF-normalized source hashes. These distinguish Windows line endings from Git blob normalization. The [notebook validation log](runs/2026-10-06/notebooks.txt) covers all 13 pairs.

| Program | Checks | Seconds | Log |
|---|---:|---:|---|
| [01_ExactBackgrounds_EDFE.wl](01_ExactBackgrounds_EDFE.wl) | 49/49 | 42.2 | [output](runs/2026-10-06/01_ExactBackgrounds_EDFE.txt) |
| [02_RenormalizedOnShellAction.wl](02_RenormalizedOnShellAction.wl) | 15/15 | 18.1 | [output](runs/2026-10-06/02_RenormalizedOnShellAction.txt) |
| [03_RadialSolutionDerivation.wl](03_RadialSolutionDerivation.wl) | 20/20 | 18.5 | [output](runs/2026-10-06/03_RadialSolutionDerivation.txt) |
| [04_GeneralBackground_LEDFE.wl](04_GeneralBackground_LEDFE.wl) | 61/61 | 76.8 | [output](runs/2026-10-06/04_GeneralBackground_LEDFE.txt) |
| [05_QuadraticRenormalization_TwoPoint.wl](05_QuadraticRenormalization_TwoPoint.wl) | 67/67 | 76.8 | [output](runs/2026-10-06/05_QuadraticRenormalization_TwoPoint.txt) |
| [06_NoetherCharges.wl](06_NoetherCharges.wl) | 31/31 | 5.1 | [output](runs/2026-10-06/06_NoetherCharges.txt) |
| [07_Worldsheet.wl](07_Worldsheet.wl) | 35/35 | 5.1 | [output](runs/2026-10-06/07_Worldsheet.txt) |
| [08_Uplift_Killing.wl](08_Uplift_Killing.wl) | 55/55 | 63.7 | [output](runs/2026-10-06/08_Uplift_Killing.txt) |
| [09_BoundaryCandidate.wl](09_BoundaryCandidate.wl) | 14/14 | 3.4 | [output](runs/2026-10-06/09_BoundaryCandidate.txt) |
| [10_AsymptoticSymmetries.wl](10_AsymptoticSymmetries.wl) | 7/7 | 3.8 | [output](runs/2026-10-06/10_AsymptoticSymmetries.txt) |
| [11_ZeroL_ExactLEDFE.wl](11_ZeroL_ExactLEDFE.wl) | 5/5 | 7.4 | [output](runs/2026-10-06/11_ZeroL_ExactLEDFE.txt) |
| [12_RiemannianVacuum_ExactLEDFE.wl](12_RiemannianVacuum_ExactLEDFE.wl) | 8/8 | 9.4 | [output](runs/2026-10-06/12_RiemannianVacuum_ExactLEDFE.txt) |
| [13_TwoPointFunctions.wl](13_TwoPointFunctions.wl) | 10/10 | 3.0 | [output](runs/2026-10-06/13_TwoPointFunctions.txt) |

This is a new run, not a reuse of the historical September 17 result. That earlier 336/336 run remains in Git history. The total includes both retained regressions and new checks; the per-equation [ledger](EQUATION_LEDGER.md) records their scope. Exact L_pm=0 and order-z general-background results are intentionally distinguished. Passing checks do not settle all interior/zero-mode prescriptions or prove the uncovered identities listed in the README.
