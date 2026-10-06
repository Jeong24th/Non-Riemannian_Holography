# Standalone Mathematica verification

Each `.wl` file contains every definition it uses. Copy **one file** into an empty directory and run it in a fresh Mathematica kernel. No shared tools, companion files, external packages, downloads or manuscript source are needed. Each `.nb` embeds that entire program and can likewise be used alone. Notebook evaluation clears `Global` and the suite's `NRH` context.

```sh
wolframscript -file 11_ZeroL_ExactLEDFE.wl
```

In Mathematica, open the matching notebook and choose **Evaluation > Evaluate Notebook**. The single input cell embeds the complete source as a string evaluated by `ToExpression`; it does not load its `.wl`. Code and notebook payloads are checked for exact agreement. A failed check exits command-line execution with status 1.

To run every file with separate kernels in PowerShell:

```powershell
Get-ChildItem *.wl | Sort-Object Name | ForEach-Object {
    wolframscript -file $_.FullName
    if ($LASTEXITCODE -ne 0) { throw "Verification failed: $($_.Name)" }
}
```

Source SHA-256: `F8AC68BB80B6C3D572160FD83705FE7506FCF277B00FB95021D39C2AF1E8E239` (2026-10-06). 139 numbered displays, including subequations.

The current SM has five sections: action; linearized dynamics and renormalization; worldsheet; uplift and symmetries; boundary candidate. Former SM2 and SM3.6 were removed from the manuscript, but the supporting radial and charge calculations remain here.

| Standalone file | Calculation and limits |
|---|---|
| [01_ExactBackgrounds_EDFE.wl](01_ExactBackgrounds_EDFE.wl) | Full R and NR EDFE; arbitrary chiral L_pm; NR equation for arbitrary W and its exact integrated solution; Letter responses. |
| [02_RenormalizedOnShellAction.wl](02_RenormalizedOnShellAction.wl) | Gamma-squared identity, endpoint contribution, volume counterterm and finite cutoff limit; negative control without CT. |
| [03_RadialSolutionDerivation.wl](03_RadialSolutionDerivation.wl) | Supplemental derivation of the Letter NR radial solution, retained after removal of former SM2. |
| [04_GeneralBackground_LEDFE.wl](04_GeneralBackground_LEDFE.wl) | Derive all ten linearized equations; solve the hierarchy; directly substitute current R/NR compact solutions through exp(-2y/l). |
| [05_QuadraticRenormalization_TwoPoint.wl](05_QuadraticRenormalization_TwoPoint.wl) | Radial momenta, ordered second variation, quadratic derivative CT cancellation modulo tangential integration by parts, particular two-point response. |
| [06_NoetherCharges.wl](06_NoetherCharges.wl) | Noether potential, field-dependent parameter subtraction, surface one-form, path integral to Q[epsilon], circle normalization and central terms. |
| [07_Worldsheet.wl](07_Worldsheet.wl) | First-order worldsheet and clock algebra; current SNCdualB/SNCreconstruction and both B-transformed frames; stated vertex/probe sectors. |
| [08_Uplift_Killing.wl](08_Uplift_Killing.wl) | Specified uplift/Killing and Clifford sectors, including the reduced one-sided NR system. |
| [09_BoundaryCandidate.wl](09_BoundaryCandidate.wl) | Classical boundary action; non-Abelian bosons and specified Grassmann realizations. |
| [10_AsymptoticSymmetries.wl](10_AsymptoticSymmetries.wl) | Generalized Lie derivatives give both R/NR asymptotic transformations and print delta L and delta W1. |
| [11_ZeroL_ExactLEDFE.wl](11_ZeroL_ExactLEDFE.wl) | Exact NR L_pm=0 solution, arbitrary W1 and sources, c_NR=0, free H_NR/chiral modes; solves the -W1^2 r z^2/16 coefficient and checks ten full-radius equations. |
| [12_RiemannianVacuum_ExactLEDFE.wl](12_RiemannianVacuum_ExactLEDFE.wl) | Exact R L_pm=0 reconstruction using Bessel solutions and two exact quadratures, fixed NS flux/nonzero momentum; ten full-radius residuals; K3 regularity and H_R matching. |
| [13_TwoPointFunctions.wl](13_TwoPointFunctions.wl) | Particular diagonal/mixed kernels, reverse ordering, R-vacuum hair self-correlator, normalization and regular exchange defect. |

## What the checks establish

The EDFE implementation computes the semi-covariant connection, curvature, scalar and mixed projector equation from the fields. It does not just compare copied answers. The general-background LEDFE are obtained by differentiating these tensors; the current compact expressions are then substituted separately. Their zero residuals are through z=exp(-2y/l), not an assertion of exact full-radius solutions at general L_pm.

For L_pm=0, file 11 checks an exact NR solution without a radial series truncation. File 12 derives exact R equations including z^2 terms, solves the reduced Bessel equation and verifies reconstruction via A'=phi and B'=A. Primitive constants remain free. The regular K3 branch on nonzero Euclidean momentum fixes the specified R-vacuum logarithmic hair response. Neither calculation silently imposes a unique NR interior condition or deletes all homogeneous data.

The on-shell action limit assumes no extra action at the interior endpoint. Quadratic counterterm checks in file 05 use tangential integration by parts and the stress constraints in both slots; they do not establish a universal nonlinear intrinsic CT. Particular two-point kernels use the manuscript's joint Lorentzian Green-function prescription; local contacts, spatial quotient and general homogeneous completion remain to be fixed. File 13 computes the R-vacuum hair coefficient from the differentiated logarithmic response, with file 12 providing its exact interior matching.

## Conventions

Doubled coordinate order: (dual x+, dual x-, dual y; x+, x-, y). Dual derivatives vanish. Lowered Gamma_CAB has unit-weight antisymmetrization and Gamma^B_BA=-2 partial_A d. The EDFE are PSPbar=0 and S_(0)=-4/l^2, equivalently G_MN=2 J_MN/l^2. Boundary eta_(+-)=-1 and barred eta_(+-)=+1. These are the manuscript's core lecture conventions.

`u` in the older exact-background code means exp(+2y/l); `z` always means exp(-2y/l). `yy` is the explicit y and DyZY=d_yy-(2z/l)d_z. `a0,b0,r0,c0,v0` are the five sources (h--,h++,h-+,h+-,delta d). `Rp,Rm,Hr,Hn` are response functions. In file 12 p,m denote eigenvalues of partial+/- with k^2=2pm>0 on the Euclidean branch. That nonzero-momentum reconstruction does not fix chiral/zero-mode sectors.

## Remaining coverage limits

The generic off-shell Codazzi identity, complete universal Box identities, exact compact metric first variation and polarized Hessian are not newly proven by the existing on-shell or ordered-asymptotic tests. The complete current real ten-dimensional spinor basis and nonzero integrable fermionic charges are not constructed. Some worldsheet checks are limited to the indicated radial, boundary or constant-profile sector. Boundary candidate Grassmann tests use the stated zero-connection realization. The [equation ledger](EQUATION_LEDGER.md) distinguishes these gaps instead of equating a passing test total with full manuscript verification.

See [REFERENCE_RUN.md](REFERENCE_RUN.md) for executable hashes, fresh-kernel results and logs, and [MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) for settled current equation numbers.

Fresh execution: **377/377 checks passed**, 13 programs, 13 standalone notebook payloads validated.
