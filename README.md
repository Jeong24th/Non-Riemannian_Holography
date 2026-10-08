# Non-Riemannian Holography verification

Reproducibility calculations supporting *Non-Riemannian Holography: Long Strings and Soft Hair* by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

The October 8 update provides **14 standalone Mathematica programs and matching standalone notebooks**, aligned with the revised manuscript. Each program embeds its tensor tools and can run alone in an empty directory. Each notebook embeds the full program and needs no companion file.

Source SHA-256: `2E6A7A8A20D581043D6457865662E76F7D46433AF29847393F1DBA4BF753F822` (2026-10-08). The settled source has 136 numbered displays, including subequations. The archive contains verification code, not the manuscript or protected coauthor material.

To reproduce the quadratic expansion of the nonperturbative intrinsic counterterm candidate in **SM29**, run:

```sh
wolframscript -file mathematica/14_SM29_IntrinsicCounterterm.wl
```

The calculation includes the intrinsic D=2 DFT Christoffel connection, boundary curvatures, measure and source-to-cutoff fields. It compares the candidate with the volume term plus the quadratic source counterterms. At quadratic order, the NR vacuum matches through finite terms modulo tangential total derivatives; the R vacuum cancels divergences but needs a finite local contact adjustment. Constant-state residuals demonstrate why this candidate does not replace the general-state counterterms. Its nonperturbative expression is not a proven universal nonlinear renormalization scheme.

The rest of the suite covers exact R/NR EDFE, asymptotic symmetries and Noether charges, exact L_pm=0 linearized solutions, general-background asymptotic LEDFE, renormalized actions, worldsheet reductions, and specified two-point functions. The latest SM2 organization places common normalization in SM2.2 and branch source-response derivations in SM2.3/SM2.4. The worldsheet conventions use real Lorentzian chiral sectors, X/bar-X non-Riemannian covectors and the Kalb–Ramond two-form.

Use the [file guide](mathematica/README.md), [SM29 walkthrough](mathematica/README.md#sm29-calculation), [revision notes](mathematica/SM_REVISION_NOTES.md), [equation ledger](mathematica/EQUATION_LEDGER.md), and [current numbering map](mathematica/MANUSCRIPT_MAP.md). Recorded kernel versions, execution totals, logs and notebook validation are in [REFERENCE_RUN.md](mathematica/REFERENCE_RUN.md) and [RUN_RESULTS.json](mathematica/RUN_RESULTS.json).

General L_pm solutions are checked through exp(-2y/l); the specified L_pm=0 solutions have full-radius residual checks. The regular R-vacuum K3 solution fixes its logarithmic hair response and separated-point self-correlator. General-state homogeneous completion, contacts, the spatial quotient and NR interior conditions remain separate inputs. Explicit source-dependent terms in H_R are retained; only the remaining integration-function dependence requires further interior data. A passing suite does not establish every manuscript equation or interpretation.

Historical calculations under `checks/` and `evidence/` remain for provenance. Their manuscript-string guards can target older sources and are not the current release entry point. Removed manuscript passages, including the static long-string energy probe, can survive as ancillary checks. No new Python verification program is published in this update. `MANIFEST.sha256` hashes tracked payload bytes after Git normalization.


Fresh execution: **451/451 checks passed**, 14 programs, 14 standalone notebook payloads validated.
