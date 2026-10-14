# Mathematica file guide

Run a standalone file with `wolframscript -file FILE.wl`, or open its `.nb` and evaluate the single input cell in a fresh kernel. Every notebook embeds its complete source using `ToExpression`; it does not load a sibling file. To run the complete archive with logging, use `python run_verification.py --mathematica` from the repository root.

The [current receipt](REFERENCE_RUN.md) distinguishes fresh Python execution, matching historical Mathematica code, and unexecuted Mathematica 09. Historical assertion names and old SM numbers inside unchanged programs are identifiers; use the [current map](MANUSCRIPT_MAP.md) for locations.

Source SHA-256: `1DADD2A6F4A65E2D142A1086CE907176B8C4677951D877F3076A6077D03C9B55`. Version `submission-2026-10-14`; 139 numbered displays including subequations, ending at SM115. Numbers follow the settled build.

| Standalone file | Calculation and limits |
|---|---|
| [01_ExactBackgrounds_EDFE.wl](01_ExactBackgrounds_EDFE.wl) | Full R/NR EDFE, arbitrary chiral L_pm, the NR equation for arbitrary W and its exact integrated solution, and Letter responses. |
| [02_RenormalizedOnShellAction.wl](02_RenormalizedOnShellAction.wl) | Gamma-squared identity, endpoint contribution, volume counterterm and finite cutoff limit; negative control without the counterterm. |
| [03_RadialSolutionDerivation.wl](03_RadialSolutionDerivation.wl) | Supplemental derivation of the Letter NR radial solution, retained after removal of the former radial-branch SM section. |
| [04_GeneralBackground_LEDFE.wl](04_GeneralBackground_LEDFE.wl) | Derives all ten linearized equations and solves their hierarchy; substitutes the current R/NR compact expansions through exp(-2y/l). |
| [05_QuadraticRenormalization_TwoPoint.wl](05_QuadraticRenormalization_TwoPoint.wl) | Radial momenta, ordered second variation, quadratic derivative counterterm cancellation modulo tangential integration by parts, and particular two-point responses. |
| [06_NoetherCharges.wl](06_NoetherCharges.wl) | Noether potential, field-dependent parameter subtraction, surface one-form, integrated Q[epsilon], circle normalization and central terms. |
| [07_Worldsheet.wl](07_Worldsheet.wl) | Exact classical finite-chi Lorentzian auxiliary reduction, X/bar-X algebra, `SNCdualB`/`SNCreconstruction`, both B-transformed frames, and stated vertex/probe sectors. The removed static winding probe is ancillary. |
| [08_Uplift_Killing.wl](08_Uplift_Killing.wl) | Specified uplift/Killing and Clifford sectors, including the reduced one-sided NR system. |
| [09_BoundaryCandidate.wl](09_BoundaryCandidate.wl) | Current free abelian RR model: nilpotency, d cancellation, Cartan/Witt identities and chiral SUSY variation. Prepared but not executed; use the fresh Python companion as current evidence. |
| [10_AsymptoticSymmetries.wl](10_AsymptoticSymmetries.wl) | Generalized Lie derivatives produce both R/NR asymptotic transformations and print delta L and delta W1. |
| [11_ZeroL_ExactLEDFE.wl](11_ZeroL_ExactLEDFE.wl) | Exact NR L_pm=0 solution, arbitrary W1 and sources, c_NR=0, free H_NR/chiral modes; solves the -W1^2 r z^2/16 coefficient and checks ten full-radius equations. |
| [12_RiemannianVacuum_ExactLEDFE.wl](12_RiemannianVacuum_ExactLEDFE.wl) | Exact R L_pm=0 reconstruction with Bessel solutions and two exact quadratures in the fixed NS-flux, nonzero-momentum sector; ten full-radius residuals, K3 regularity and H_R matching. |
| [13_TwoPointFunctions.wl](13_TwoPointFunctions.wl) | H_R decomposition and R stress/source-response checks, NR source-generator substitution, particular kernels and normalization, reverse ordering, R-vacuum hair self-correlator, and the regular exchange defect. |
| [14_SM29_IntrinsicCounterterm.wl](14_SM29_IntrinsicCounterterm.wl) | Quadratic expansion of current SM30 (historical filename SM29) from intrinsic D=2 geometry; NR/R vacuum comparisons including finite terms, and constant-state failure of the candidate. |

## Reading the results

SM2 uses S_ren throughout, with the connected Hessian factor 1/(4i) set by the +2 j O source coupling. Historical variable names for quadratic expressions do not introduce another factor. Files 05 and 13 retain the same numerical normalization.

File 14 expands the intrinsic SM30 candidate (historical filename SM29), including its volume term, measure and source-to-cutoff fields, to quadratic order. On the NR vacuum it matches the reference counterterms through finite terms modulo tangential total derivatives. On the R vacuum a finite local contact adjustment is needed. Constant-state residuals rule out interpreting this candidate as a universal nonlinear subtraction scheme. Current source counterterms are SM28/SM29; the common Hessian normalization is SM27.

File 08 retains the same Clifford representation (historical Gamma now gamma, and the three-dimensional gamma now lambda in SM4). Its reduced unprimed calculations do not independently prove the full current real internal product basis, every primed solution or all global polarization counts. File 07 contains ancillary checks of a removed winding-energy example; they are not current manuscript claims.

File 09 uses the SM5 gamma representation, which is displayed separately from the SM4 three-dimensional representation. There is one unrestricted bispinor B, no non-Abelian gauge connection and no zeta symmetry. Its Mathematica implementation checks a subset of the Python companion's identities. Neither implementation claims generic curved-background SUSY or a complete quantum background-response prescription.
