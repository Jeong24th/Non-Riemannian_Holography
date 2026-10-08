# Standalone Mathematica verification

Each `.wl` file contains every definition it uses. Copy one file into an empty directory and run it in a fresh Mathematica kernel. No shared tools, companion files, external packages, downloads or manuscript source are needed. Each `.nb` embeds the entire matching program and can likewise run alone. Notebook evaluation clears `Global` and the suite's `NRH` context.

```sh
wolframscript -file 14_SM29_IntrinsicCounterterm.wl
```

In Mathematica, open the matching notebook and choose **Evaluation > Evaluate Notebook**. Its single input cell evaluates the embedded source string with `ToExpression`; it does not load the `.wl`. The execution record includes code/notebook payload validation. A failed check exits command-line execution with status 1.

To run every file with separate kernels in PowerShell:

```powershell
Get-ChildItem *.wl | Sort-Object Name | ForEach-Object {
    wolframscript -file $_.FullName
    if ($LASTEXITCODE -ne 0) { throw "Verification failed: $($_.Name)" }
}
```

Source SHA-256: `2E6A7A8A20D581043D6457865662E76F7D46433AF29847393F1DBA4BF753F822` (2026-10-08). There are 136 numbered displays, including subequations; the map was checked against the settled build of this source.

The SM has five sections: action; linearized dynamics and renormalization; worldsheet; uplift and symmetries; boundary candidate. Common two-point normalization is in SM2.2, followed by the R derivation in SM2.3 and NR derivation in SM2.4. The former independent SM2.5 and its duplicate kernel display were removed. The final kernels remain in Letter (20), `Rcorrelators`.

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
| [09_BoundaryCandidate.wl](09_BoundaryCandidate.wl) | Classical boundary action, non-Abelian bosons and specified Grassmann realizations. |
| [10_AsymptoticSymmetries.wl](10_AsymptoticSymmetries.wl) | Generalized Lie derivatives produce both R/NR asymptotic transformations and print delta L and delta W1. |
| [11_ZeroL_ExactLEDFE.wl](11_ZeroL_ExactLEDFE.wl) | Exact NR L_pm=0 solution, arbitrary W1 and sources, c_NR=0, free H_NR/chiral modes; solves the -W1^2 r z^2/16 coefficient and checks ten full-radius equations. |
| [12_RiemannianVacuum_ExactLEDFE.wl](12_RiemannianVacuum_ExactLEDFE.wl) | Exact R L_pm=0 reconstruction with Bessel solutions and two exact quadratures in the fixed NS-flux, nonzero-momentum sector; ten full-radius residuals, K3 regularity and H_R matching. |
| [13_TwoPointFunctions.wl](13_TwoPointFunctions.wl) | H_R decomposition and R stress/source-response checks, NR source-generator substitution, particular kernels and normalization, reverse ordering, R-vacuum hair self-correlator, and the regular exchange defect. |
| [14_SM29_IntrinsicCounterterm.wl](14_SM29_IntrinsicCounterterm.wl) | Quadratic expansion of SM29 from intrinsic D=2 geometry; NR/R vacuum comparisons including finite terms, and constant-state failure of the candidate. |

## SM29 calculation

SM29 has stable label `SMboundarycurvaturecandidate`. Its nonperturbative expression contains the volume term, intrinsic scalar curvature, two curvature-squared combinations and the scalar box term. File 14 expands this expression through second order in boundary fluctuations. This is the order required for the quadratic action and linear response of its variation.

The curvatures are built intrinsically from the induced tangential generalized metric and dilaton using the **D=2 DFT Christoffel connection** and tangential derivatives. Restricting the bulk curvature to the cutoff would give a different calculation. The scalar box is the O(D,D) symmetric divergence operator stated in the manuscript. The measure and full source-to-cutoff relation are included before the asymptotic comparison.

The reference quadratic counterterms are SM27/SM28, `SMderivativectNR` and `SMderivativectR`, together with the volume term SM23. SM29 already includes that volume term; it must not be added a second time. In the R comparison, the NR counterterm polynomial is evaluated on the R fluctuation fields, as prescribed in the manuscript.

At quadratic order on the NR vacuum L_pm=W1=0, the intrinsic candidate and the reference counterterms agree in every divergent and finite term as Y tends to infinity, modulo tangential total derivatives. On the R vacuum L_pm=0, divergences cancel but a finite local adjustment is needed to use the reference contact convention. Such a local adjustment affects contacts, not separated-point kernels. The constant-state checks retain arbitrary constant L_pm (and constant W1 on the NR branch) and expose uncanceled divergences on both branches; they do not establish a general nonlinear completion. For general states the manuscript retains SM27/SM28.

The explicit results use the boundary sources

$$
a=h^{(0)}_{\ominus\bar\ominus},\quad
b=h^{(0)}_{\oplus\bar\oplus},\quad
r=h^{(0)}_{\ominus\bar\oplus},\quad
c=h^{(0)}_{\oplus\bar\ominus},\quad v=\delta d^{(0)}.
$$

These are `a0,b0,r0,c0,v0` in the code. Define the retained nondecaying density difference by

$$
S_{\rm ct}^{\partial,(2)}-S_{\rm ct,ref}^{(2)}
=-\frac{1}{32\pi G}\int d^2x\,\Delta\mathcal L+o(1),\qquad Y\to\infty,
$$

where the reference includes the volume term and the appropriate SM27/SM28 source counterterm. Modulo tangential total derivatives, the NR-vacuum nondecaying difference is zero. The R-vacuum difference is the finite local polynomial

$$
\begin{aligned}
\Delta\mathcal L_{{\rm R},{\rm vac}}
=\frac l8\big[&l^2a\,\partial_+^3\partial_-r
+l^2b\,\partial_+\partial_-^3r
-8l^2r\,\partial_+^2\partial_-^2v\\
&-16a\,\partial_+^2v-16b\,\partial_-^2v
+64v\,\partial_+\partial_-v\big].
\end{aligned}
$$

For constant states, the divergent part is the same on both branches:

$$
\begin{aligned}
D_{\rm rem}={}&-2Y(L_+a+L_-b)\,\partial_+\partial_-r\\
&+\left(\frac{lY^2}{2}-\frac{l^2Y}{4}\right)
r\left(L_+\partial_+\partial_-^3+L_-\partial_+^3\partial_-\right)r.
\end{aligned}
$$

There is no remaining power divergence in the difference, but these Y and Y^2 terms are generically nonzero. A periodic cosine-source witness computed from the residual has nonzero integral, so the obstruction is not a discarded tangential total derivative. This constant-state calculation does not assert the displayed residual for nonconstant L_pm or W1.

File 14 prints the linear and quadratic scalar curvatures, linear mixed Ricci tensor, free-cutoff density, vacuum difference, finite R contact mismatch and constant-state divergent residual. Its densities omit the common factor -1/(32 pi G). The common volume term is checked separately and canceled identically when taking the difference. Here `u=Exp[-2 Y/l]`; the source-amplitude expansion is performed before the small-u expansion, retaining u^-1 and u^0 with every power of Y.

After evaluating the notebook, the results remain accessible as associations:

```wolfram
nrResult["VacuumDifference"]   (* zero, modulo tangential total derivatives *)
rResult["VacuumDifference"]    (* finite local contact mismatch *)
rResult["Divergence"]          (* constant-state obstruction; same on NR *)
rResult["ScalarLinear"]
rResult["ScalarQuadratic"]
rResult["MixedRicciLinear"]
```

`CandidateFree` is the quadratic candidate density after the common volume term is removed; `SourceCountertermFree` is the branch source-counterterm density. To match the R-vacuum contact convention, subtract the printed finite mismatch from the measure-inclusive quadratic density, or equivalently add the finite action term $+(32\pi G)^{-1}\int d^2x\,\Delta\mathcal L_{{\rm R},{\rm vac}}$. All comparisons by integration by parts assume compactly supported or periodic tangential variations.

## What the other checks establish

The EDFE implementation computes the DFT Christoffel connection, curvature, scalar and mixed projector equation from the fields. The general-background LEDFE are obtained by differentiating these tensors; current compact expressions are then substituted separately. Their vanishing residuals are through z=exp(-2y/l), without a claim of exact full-radius solutions at general L_pm.

For L_pm=0, file 11 checks an exact NR solution without a radial-series truncation. File 12 derives exact R equations including z^2 terms, solves the reduced Bessel equation and verifies reconstruction through A'=phi and B'=A. Primitive constants remain free. On the chosen nonzero Euclidean-momentum branch, K3 regularity fixes the specified R-vacuum logarithmic hair response; Lorentzian correlators follow by the manuscript's continuation and prescription. These calculations do not select a unique NR interior condition.

The on-shell action limit assumes no extra action at the interior endpoint. File 05 checks quadratic counterterms using tangential integration by parts and stress constraints in both slots. File 14 supplies the separate intrinsic-candidate comparison and its limitations.

The two-point calculation differentiates one source at a time. In the R branch, the explicit b_pm+omega_pm and local terms in H_R determine particular responses; the remaining F_R and other integration functions require an interior prescription for any additional source dependence. In the NR branch, the chosen generalized diffeomorphism gives the stated mixed response, but does not fix the full hair self-response. File 13 applies the joint Lorentzian Green-function prescription and excludes local contacts. General homogeneous completion and the spatial quotient remain to be specified. The diagonal particular kernels have no singular exchange defect; a regular defect remains for nonconstant L_pm and must be canceled by the completed response.

## Conventions

Doubled bulk coordinate order is (dual x+, dual x-, dual y; x+, x-, y), with vanishing dual derivatives. Lowered Gamma_CAB denotes the DFT Christoffel connection, with unit-weight antisymmetrization and Gamma^B_BA=-2 partial_A d. The EDFE are PSPbar=0 and S_(0)=-4/l^2, equivalently G_MN=2 J_MN/l^2. Boundary eta_(+-)=-1 and barred eta_(+-)=+1. SM29 uses the induced four-dimensional doubled boundary with D=2 and its own intrinsic curvatures.

`u` in older exact-background code means exp(+2y/l); `z` in those files means exp(-2y/l). `yy` is the explicit y and DyZY=d_yy-(2z/l)d_z. `a0,b0,r0,c0,v0` denote the five sources (h--,h++,h-+,h+-,delta d). `Rp,Rm,Hr,Hn` are responses. File 14 defines its own cutoff variable and source symbols explicitly. In file 12, p,m are eigenvalues of partial+/- with k^2=2pm>0 on the Euclidean branch; this reconstruction does not fix chiral/zero-mode sectors.

SM3 uses real Lorentzian worldsheet null coordinates, epsilon^(01)=1 and epsilon^(+-)=-1; bars distinguish real chiral sectors. The non-Riemannian covectors are X and bar-X. The stable equation label `SNCtau` and legacy tauP/tauM variable names refer to these forms; Pauli tau_i in SM4 are unrelated. B denotes the Kalb–Ramond two-form. File 07 expresses radial fusion per real chiral separation; these chiral-algebra checks do not construct the full Lorentzian path integral.

## Remaining coverage limits

The generic off-shell Codazzi identity, complete universal box identities, exact compact metric first variation and polarized Hessian are not proved by the on-shell or ordered-asymptotic tests. The complete current real ten-dimensional spinor basis and nonzero integrable fermionic charges are not constructed. The classical finite-chi worldsheet auxiliary reduction is checked, while quantum equivalence, the full causal prescription and an all-order dressed W1 deformation are not derived. Some vertex/probe checks apply only in the stated radial, boundary or constant-profile sector. Boundary candidate Grassmann checks use the stated zero-connection realization. The NR two-point source-generator substitution does not independently rederive the full radial-gauge completion. The [equation ledger](EQUATION_LEDGER.md) records these restrictions individually.

[REFERENCE_RUN.md](REFERENCE_RUN.md) and [RUN_RESULTS.json](RUN_RESULTS.json) contain execution evidence and notebook payload validation. [MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) gives current equation numbers and identifies removed labels. Historical Python files under `checks/` and `evidence/` are ancillary and can contain guards for earlier manuscript versions.


Fresh execution: **451/451 checks passed**, 14 programs, 14 standalone notebook payloads validated.
