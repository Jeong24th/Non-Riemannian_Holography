# Equation ledger

Current source SHA-256: `C733463D16EA773C6DDA09143E44C1823F7C3CD46F7D771F14E53C8ED722C2DE` (2026-09-11).

Letter (1)-(22); Supplemental Material (1)-(202).

The labels and numbers below refer to the current source. Coverage refers only to
the public software, not to unpublished calculation files. The SM1 rewrite and the
current real Killing-spinor basis postdate the recorded 270-check Mathematica run.
An unchanged or reused label is not by itself evidence that a new formula is checked.

**V**: symbolic check of the stated relation; **V°**: only a stated consequence,
background or asymptotic order; **D**: definition; **S**: cited statement;
**N**: no current public check. The actual execution record is in
[REFERENCE_RUN.md](REFERENCE_RUN.md).

Quoted check names retain the executable files' historical wording. Any equation
numbers inside those quotes refer to the older source, not to the current map.

The general-source EDFE coefficients and constrained finite-slice recursion are not
implemented by the public suite. Its response checks cover the common-vacuum
linearization, moving-frame finite momenta on both backgrounds, and the particular
Ward kernels. They do not determine the remainder, hair self-response, local
counterterms, interior conditions or source zero modes.

The current one-sided Killing spinors have arbitrary chiral functions, including
infinitely many periodic modes. These counts refer to real polarizations, not Fourier modes.
The public reduced jet checks do not construct the current full real basis or any
fermionic phase-space charge. The older NRH07 label mentioning an S3 line and zeta+
must not be read as a check of the current `SMinternalprojectors` formula.


## Letter

| No. | Label | Content | Public coverage | Status |
|---:|---|---|---|:---:|
| 1 | `Rfields` | NS–NS Bañados family (frame ϑ±, metric, B, dilaton) | NRH02 "H J H = J (O(3,3) constraint)" (3 checks) (built from (1)); the EDFE and horizon checks; NRH06 "the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+" (2 checks); NRH07 "R(AdS3, Banados) = -6/l^2 and H^2 = -24/l^2 for arbitrary chiral L_pm (pairwise cancellation with S3…" | V |
| 2 | `RDFTfields` | DFT variables H_MN and e^{-2d} for the Riemannian family | NRH02 verifies the O(3,3) constraint, dilaton determinant and EDFE. | V |
| 3 | `Rboundary` | asymptotic falloffs | NRH02 "H − H^infty = O(e^{−2y/l})", "d + y/l = O(e^{−4y/l})" | V |
| 4a | `boundaryH` | H^∞, d^∞ | NRH02 "boundaryH…" (two checks) | V |
| 4b | `Rboundaryframe` | aligned D = 2 boundary frame | NRH02 "Rboundaryframe: boundary vielbeins reproduce P^(0) and Pbar^(0)" | V |
| 5 | `NRvariables` | (Π, q, e^σ, χ) and e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) | NRH03 "e^{sigma} sinh(chi/2) = …", "NRradialchange: d chi/d q…" | V |
| 6 | `NRHcompact` | the everywhere non-Riemannian H_MN | NRH03 "H J H = J for the exact non-Riemannian matrix (generic W)", "type (1,1) at every radius"; EDFE checks; NRH07 "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)" | V |
| 7 | `NRdilaton` | e^{−2d} = e^{2y/l}(1 − q²), d = −y/l + ln cosh(χ/2√2) | NRH03 "e^{-2d} e^{-2y/l} = 1 − q^2…"; "d chi/dy…" | V |
| 8 | `NRWgeneral` | exact hair profile W | NRH03 "EDFE scalar: S_(0) = −4/l^2 for ARBITRARY W", "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F", "G_MN = 2 l^-2 J_MN on the ODE shell", "NRWgeneral] solves d^2W/dchi^2 = F", "W_0 and W_1 multiply the two homogeneous modes", "near the boundary…", "falloff bookkeeping…" | V |
| 9 | `MainGprofile` | radial profile G(χ) = 4√2∫₀^χ[sinh t/sinh(t/√2) − √2]dt, G(0) = G′(0) = 0 | NRH03 "NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)", "NRGprofile: G'(0) = 0" | V |
| 10 | `RGKPW` | GKPW relation Z_DFT = Z_CFT | — (defining relation) | D |
| 11 | `Rrenvariation` | On-shell boundary variation of the renormalized action | NRH04 verifies the bulk cancellation and leading counterterm on constant BTZ for arbitrary tangential h and delta d; NRH02 checks one-point normalization. | V° |
| 12 | `Mainmomenta` | Generalized-metric and scalar radial momenta A and B | NRH04 checks their role in the Gamma-squared variation on constant BTZ with arbitrary tangential h and delta d; NRH05 checks B^y on R/NR. | V° |
| 13a | `RKdef` | ⟨K⟩, ⟨T₍₀₎⟩ from the rescaled momenta | NRH02 "RKdef: coefficient matching -2K = (16 pi G)^{-1} A^y and 2T_(0) = (16 pi G)^{-1} 2(B^y + 4/l)"; **nonlinear evaluation on both exact saddles**: NRH02 "on the exact family: −(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L−},{1, L−}}/(16 pi G l)", "…e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0"; NRH04 "on the exact NR family: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, W1/4},{0, L-}}/(16 pi G l)…" (2 checks), "…hence <T_(0)> = 0"; NRH04 "SMscalarcutoffresponse: 2(16 pi G)^{-1} e^{-2d}(B^y + 4/l) with B^y = 4 d_y d = -4/l + 4 d_y delta d…" (2 checks) | V |
| 13b | `RDFTconservation` | T^CFT_AB and ∇^A T_AB = 0 | NRH02 "div T = {…} exactly"; "SMinvariantstress…" | V |
| 14 | `Rcontinuity` | Ward identities | NRH02 "div T…", "no local condition on K_{op bom}" | V |
| 15 | `Rkilling` | common R/NR asymptotic-symmetry generator (ξ^±_R with the l²e^{−2y/l}∂²ε/4 tail, ξ^±_NR = ε^±; common radial and dual components) | NRH02 "Lhat_xi d = O(e^{−4y/l})", "Lhat_xi H − delta_eps H = O(e^{−4y/l})"; NRH05 full generator checks; the NR form ξ^± = ε^± underlies the NRH03 checks of (17) | V |
| 16 | `RVirasoro` | δ_εL± with −(l²/4)∂³ε | NRH02 "Lhat_xi H - delta_eps H = O(e^{-4y/l})", "c = 3l/2G reproduces delta_eps T = eps T' + 2 T eps' - (c/12) eps'''"; NRH05 "Virasoro cocycle: (charge bracket density) + (l^2/4) e1 e2''' is a total derivative", "(iii) R bracket − adjoint − central…" | V |
| 17 | `NRasympt` | δ_εL±, δ_εW₁ at W₀ = 0 | NRH03 "NRasympt…" (four checks) | V |
| 18 | `Mainonepoints` | R/NR one-point matrices and vanishing scalar responses | NRH02 and NRH04 evaluate the radial momenta on the two backgrounds. | V |
| 19 | `Mainframevariation` | Frame variations in the gauge with no local Lorentz rotation | NRH04 checks both projector variations for arbitrary mixed h. | V |
| 20 | `Mainfivebyfive` | connected 5x5 matrix: particular A/M entries + R_s + local contacts | Definition; NRH04 checks the A/M entries. R_s includes hair self-response and other Ward-undetermined terms. | D |
| 21 | `Rcorrelators` | particular A_NR, A_R, M_NR, M_R kernels, with coefficients at x | NRH04 "Rcorrelators:" checks: Ward operators applied to the chosen inverse derivative; source-sign negative control. | V |
| 22 | `Mainworldsheet` | ℒ_GO = (1/2πα′)(β∂̄x⁺ + β̄∂x⁻ + ∂y∂̄y), V_𝒲 = (𝒲/4πα′)∂x⁺∂̄x⁻ | NRH06 "at chi -> 0 the W coupling reduces to…", "with the 1/(2 pi alpha') prefactor this is V_W = …" | V |

## Supplemental Material

| No. | Label | Content | Public coverage | Status |
|---:|---|---|---|:---:|
| 1 | `SMprojectors` | P, P̄ | NRH01 (used throughout); NRH02 "Gamma compatibility: nabla_C P_AB = 0" | D |
| 2 | `SMmixedfluctuation` | h_{pq̄} := δH_MN V^M_p V̄^N_q̄ | — | D |
| 3 | `SMcosetreconstruction` | δH = 2V_(M^p V̄_N)^q̄ h_pq̄ | NRH04 "delta(V eta V^T) = delta H/2…"; the constrained families of NRH04 (H(r,w), H4) | V |
| 4 | `SMframevariation` | gauge-fixed δV_{Mp} = ½V̄_M^q̄ h_pq̄, δV̄ = −½V h | NRH04 "delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2 for generic h_{p qbar}" | V |
| 5 | `SMresponsedef` | K_{pq̄}, T₍₀₎ definitions | — | D |
| 6 | `SMbackgroundconnection` | Torsionless DFT connection at D=3 | NRH01 GammaDFT; used in all curvature and momentum calculations. | V |
| 7 | `SMconnectionvariation` | δΓ_KMN in terms of δP and δd (Eq. (2.56) of Park:2025core) | — (quoted from Ref. Park:2025core; not re-derived in this archive) | S |
| 8 | `SMsixprojectors` | six-index projectors 𝒫, 𝒫̄ (Eq. (2.52) of Park:2025core) | — (definition) | D |
| 9 | `SMflatmetrics` | η, η̄ | NRH04 "V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty"; NRH07 "rep: {gamma^p, gamma^q} = 2 eta^{pq} (3d lightcone blocks)" | V |
| 10 | `SMinfinityvielbein` | limiting D = 3 frame | NRH04 "V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty", "V and Vbar are mutually orthogonal"; NRH07 "at hh = 0 and W_0 = 0, lowering the local indices with eta, etabar gives the limiting frame" | V |
| 11 | `SMFG` | Fefferman–Graham gauge | imposed in every NRH04 linearization | D |
| 12 | `SMsources` | Four crossed lower-index source couplings | NRH04 source-index raising and factor-two normalization; no extra source symbol is needed. | V |
| 13 | `SMW0completion` | Transport of fields and frames by the closed W0 B shift | NRH03/NRH06 verify the metric B shift; the full transported general-source system is not evaluated. | V° |
| 14 | `SMW0variation` | Product rule when the W0 source varies | Definition of the tangent to the B-transformed family. | D |
| 15 | `SMbackgroundexpansion` | General R/NR dilatons and derivative-dependent W expansion | NRH03 checks the radial expansion; NRH04 uses the finite-response orders. | V° |
| 16 | `SMexactboxEDFE` | Exact linearized scalar and mixed EDFE via the universal box | Cited universal-box identity; the public suite does not evaluate the full current general-source system. | S |
| 17 | `SMboxdefinition` | Universal box and section-condition identities | Definition from the cited universal-box paper. | D/S |
| 18 | `SMmixedboxdefinition` | Full mixed-tensor box with curvature and connection terms | Cited operator, not an independent public component verification. | S |
| 19 | `SMboxcurvature` | Connection field strength and its Ricci contraction | Definition. | D |
| 20 | `SMexactRdata` | Exact Riemannian fields and lower-index saddle frames | Background fields are checked in NRH02; this complete current frame recipe is not separately checked. | V° |
| 21 | `SMexactNRdata` | Exact NR substitution with full derivative-dependent W | NRH03 checks the nonlinear solution and NRH07 checks its frame; the full linearized operator is not evaluated here. | V° |
| 22 | `SMexactcomponentrecipe` | Definition of the ten component equations | Definition. | D |
| 23 | `SMexactfieldorder` | Five-field and five-evolution-equation order | Definition. | D |
| 24 | `SMexactnormalform` | Radial normal form and its principal matrix | No general all-radius component check in the public suite. | N |
| 25 | `SMexactcoefficientpolynomial` | Polynomial extraction of exact operator coefficients | Definition; its substitution has not been evaluated publicly. | D/N |
| 26 | `SMexactTaylor` | Taylor coefficients at a regular finite slice | Definition. | D |
| 27 | `SMexactTaylorSolution` | Constrained local radial Taylor recursion | No convergence or global boundary-to-interior solution is verified by this suite. | N |
| 28 | `SMexactCauchyconstraints` | Five radial constraints | No current general-background public check. | N |
| 29 | `SMexactconstraintpropagation` | Linearized Bianchi identity and constraint propagation | Cited identity; no current general radial propagation check. | S |
| 30 | `GammaDFT` | L_Γ² = e^{−2d}S₍₀₎ − ∂(e^{−2d}B) = e^{−2d}(PP − P̄P̄)(ΓΓ…) | NRH05 "on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)", "on NR (arbitrary W)…" (the toolbox Γ² density is the displayed quadratic form) | V |
| 31 | `variation` | Variation of the Gamma-squared Lagrangian | NRH04 verifies the identity on constant BTZ in its Riemannian frame, for arbitrary tangential h and delta d in radial gauge. Radial h components are not varied. | V° |
| 32 | `defB` | B^K and A^K_{p qbar} | NRH04 checks the Gamma-squared boundary variation on constant BTZ with arbitrary tangential h and delta d; NRH05 checks B^y on R/NR. | V° |
| 33 | `SMrenvariation` | On-shell variation of the action and leading counterterm at the cutoff | NRH04 checks the on-shell cancellation and counterterm variation on constant BTZ. | V° |
| 34 | `SMradialdictionary` | Finite-part prescription for the radial one-point responses | NRH02/NRH04 check coefficient matching; derivative counterterms and finite local terms are not determined. | D/V° |
| 35 | `SMsecondvariation` | Independent tangent variations, frame and measure terms | NRH04 checks the independent h arguments and moving-frame momentum; derivative counterterms are unspecified. | V° |
| 36 | `SMmomentumvariation` | δA^K_{pq̄}, δB^K including the frame-variation terms; 𝒜^K_MN | NRH04 "SMmomentumvariation:" differentiates both frame projections of A and all three terms of B on the vacuum; a frozen-frame negative control is nonzero. | V° |
| 37 | `SMbackgroundmomenta` | Leading R/NR mixed momenta and scalar momentum | NRH02/NRH04 background one-point evaluations. | V° |
| 38 | `SMsamechiralitymomenta` | Same-chirality tangential momenta | NRH04 checks the moving-frame consequence; not a separate general all-radius identity. | V° |
| 39 | `SMonept` | Riemannian one-point matrix | NRH02 nonlinear radial-momentum evaluation. | V |
| 40 | `SMNRonept` | Non-Riemannian one-point matrix | NRH04 nonlinear radial-momentum evaluation. | V |
| 41 | `SMstatefluctuations` | Source-free tangent variations of the two saddle families | The legacy fixed-frame finite differences are different objects; no dedicated public tangent check. | N |
| 42 | `SMoperatorseries` | Operator expansion through exp(-2y/l) | Definition of the asymptotic expansion. | D |
| 43 | `SMradialansatz` | Independent sources and radial logarithmic coefficients | Definition of the field ansatz; not a solved boundary-value problem. | D |
| 44 | `SMNRradialintegration` | Radial derivatives of u times a polynomial | NRH04 checks the elementary radial operator and logarithmic modes. | V° |
| 45 | `SMtotalorderhierarchy` | Total-order operator/field hierarchy | The new general-background hierarchy is not implemented in the public suite. | N |
| 46 | `SMcoupledorders` | Combined equations at orders u0 and u1 | No current general-background public check. | N |
| 47 | `SMmixingexample` | R type-changing coefficient example | No current general-background public check. | N |
| 48 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 49 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 50 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 51 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 52 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 53 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 54 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 55 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 56 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 57 | — | Common leading EDFE component | No check of this current formula in the public suite. | N |
| 58 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 59 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 60 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 61 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 62 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 63 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 64 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 65 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 66 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 67 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 68 | `SMRoperatorrule` | Relation between the displayed R and NR first-order operators | No current general-background public check. | N |
| 69 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 70 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 71 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 72 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 73 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 74 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 75 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 76 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 77 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 78 | — | State-dependent EDFE component | No check of this current formula in the public suite. | N |
| 79 | `SMleadingdilaton` | Leading dilaton slope and type-changing logarithms | NRH04 common-vacuum equations check the leading coefficients; no full general-source completion. | V° |
| 80 | `SMleadingstress` | Leading diagonal logarithmic slopes | NRH04 common-vacuum radial solution. | V° |
| 81 | `SMNRsol` | Leading hair-channel polynomial coefficients | NRH04 common-vacuum solution; this label no longer denotes the complete old radial solution. | V° |
| 82 | `SMlogexample` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 83 | `SMRh2L_omopb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 84 | `SMRdd2L` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 85 | `SMRh2L_opopb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 86 | `SMRh2L_omomb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 87 | `SMRdd2` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 88 | `SMRh2L_opomb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 89 | `SMNRh2L_omopb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 90 | `SMNRdd2L` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 91 | `SMNRdd2` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 92 | `SMNRh2L_opopb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 93 | `SMNRh2L_omomb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 94 | `SMNRh2L_opomb` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 95 | `SMtypeconstraint` | Integrated type-changing constraints with source zero modes | The arbitrary-source general-background constraints are not in the public suite. | N |
| 96 | `SMRconstraint0` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 97 | `SMRconstraint1` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 98 | `SMNRconstraint0` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 99 | `SMNRconstraint1` | General-background logarithmic coefficient or constraint | Not implemented in the public suite. | N |
| 100 | `SMgeneralradialsolution` | Remaining integration functions and inverse derivatives | Definition of the unresolved response data; no interior prescription. | D/N |
| 101 | `SMlogfreeconditions` | Restricted fixed-dilaton ensemble on the common vacuum | NRH04 checks the vacuum logarithmic coefficients; not a restriction on general source derivatives. | V° |
| 102 | `SMframemomentumcheck` | Contribution of both moving-frame projections | NRH04 SMmomentumvariation checks, including the frozen-frame negative control. | V° |
| 103 | `SMconnectionmomentumcheck` | Projected connection contribution through u | NRH04 finite-state momentum variation through exp(-2y/l). | V° |
| 104 | `SMlinearizedmomenta` | Off-shell tangential momentum variations through u | NRH04 retains arbitrary tangential fluctuations and both finite-state backgrounds. | V° |
| 105 | `SMfinitepartexplicit` | Finite coefficients after radial differentiation | NRH04 logarithmic finite-part extraction; omitted local counterterms are not fixed. | V° |
| 106 | `SMgeneralresponsevariation` | Linear responses including normalizable logarithms | NRH04 finite-part extraction; not the complete source-to-response map. | V° |
| 107 | `SMstresshessian` | Diagonal and scalar connected-response normalization | NRH04 source factor and momentum normalization; full kernels remain conditional. | V° |
| 108 | `SMradialhessian` | Five source derivatives defining the connected matrix | Definition at separated points. | D |
| 109 | `SMstresssourceexample` | Independent stress-source constraints | NRH04 verifies the Ward-derived particular consequence, not the new bulk component derivation. | V° |
| 110 | `SMgeneralWardoperators` | Boundary differential operators for stress and hair | NRH02/NRH03 asymptotic transformation laws, NRH04 particular kernels. | V |
| 111 | `SMstressparticular` | Particular stress kernels from an inverse derivative | NRH04 Rcorrelators checks. | V |
| 112 | `SMindependentsourcegenerator` | Radial compensator keeping all other leading sources zero | The new nonchiral bulk completion is not checked by the public suite. | N |
| 113 | `SMhairparticularsolution` | Full nonchiral hair response including its local term | NRH04 checks the separated-point Ward result only; the local completion is not covered. | V° |
| 114 | `SMgeneralparticularkernels` | R and NR particular mixed kernels | NRH04 Rcorrelators checks use M directly. | V |
| 115 | `SMPBHkernel` | Local-plane Green-function convention | Chosen normalization; no Lorentzian/state or cylinder zero-mode prescription is selected. | D |
| 116 | `SMgeneralpositionkernels` | Position-space particular A and M poles | NRH04 Rcorrelators checks and source-sign negative control. | V |
| 117 | `SMRtwopt` | Riemannian plane-vacuum stress normalization | NRH02/NRH04 reproduce c=3l/(2G). | V |
| 118 | `SMNRhairhessian` | Hair self-response as the derivative of an undetermined function | Normalization only; the source dependence is not computed. | D/V° |
| 119 | `SMCPSform` | surface-charge one-form | NRH05 `chargeOneForm` (all three terms), "lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary" (3 checks), "(i) integrability…" | V |
| 120 | — | Q_R | NRH05 "lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''" (2 checks) with "(iii) normalization chain: (16 pi G)^{-1} (4/l) = 1/(4 pi G l), the Letter's charge normalization" | V |
| 121 | `SMRpotentialcomponents` | lim e^{−2d}K̂ components (the minus-sector check uses K̂^{y+} = −K̂^{+y}) | NRH05 "lim e^{-2d} Khat^{-y}[eps+] = …", "…Khat^{y+}[eps-] = … (mirror)" | V |
| 122 | `SMCPSresult` | Θ̂ → 0, k = (4/l)εδL | NRH05 "lim e^{-2d} Thetahat^{+,-,y} = 0", "k^{-y}[eps+] = (4/l) eps+ dL+", "…k^{+y}[eps-…" | V |
| 123 | `SMNRchargefalloffs` | state-dependent falloffs | NRH05 "state-dependent falloffs…" | V |
| 124 | `SMNRchargecancellation` | componentwise W₁ cancellation | NRH05 "W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise" | V |
| 125 | `SMNRcharge` | δQ, Q | NRH05 "(i) k^{-y}[eps+] = delta[(4/l) eps+ L+]", "(i) mirror…" | V |
| 126 | `SMCPSalgebra` | {Q,Q} = Q[[ε,η]], c_charge = 0 (+ footnote) | NRH05 "the same-chirality C-bracket closes up to a closed B-gauge parameter (slot x~+ only)" (4 checks) (four checks), "(ii) {Q[e1+], Q[e2+]} - Q[[e1,e2]] is a total derivative: centerless plus sector" (3 checks), "(iv) the central cocycle is antisymmetric modulo total derivatives" (3 checks), "(v) the adjusted mixed bracket has a leftover (the raw closure fails, as it should)" (2 checks), "(vi) minus-sector C-bracket closes up to a closed B-gauge parameter (slot x~- only)" (2 checks), "footnote: with delta T = eps T' + 2 T eps' - l^2 eps''' the combination L + T/4 obeys the law with -…" | V |
| 127 | `SMgamma2` | Γ² identity | NRH05 "on R…", "on NR (arbitrary W)…" | V |
| 128 | `SMgamma2flux` | B^y = 4∂_y d, e^{−2d}B^y = −(4/l)(e^{2y/l} + μe^{−2y/l}) | NRH05 "on R…", "on NR…" (three checks) | V |
| 129 | `SMmudefinition` | μ dichotomy | NRH05 "on NR: −2 d_y e^{-2d} = … [mu-dichotomy]" | V |
| 130 | `SMgamma2cutoff` | S_ren(Y) | NRH05 "the regulated combination equals…" | V |
| 131 | `SMgamma2value` | S_ren = −8√μ/(16πGl)∫d²x | NRH05 "Y -> Infinity limit gives…" | V |
| 132 | — | Killing horizon e^{4y/l} = L₊L₋ | NRH02 "horizon: e^{-2d} = 0 at e^{4y/l} = L+ L−"; NRH05 "endpoints…" | V |
| 133 | `NRhill` | ψ±, A± = ¼(∂lnL)² − ½∂²lnL = ψ″/ψ | NRH03 "NRhill: …"; A± enters the verified source F | V |
| 134 | `NRradialchange` | q, χ, ∂_y, ∂_q | NRH03 "NRradialchange: d chi/d q = …", "d chi/dy = …" | V |
| 135 | `NRradialoperator` | radial operator identity | NRH03 "NRradialoperator…" | V |
| 136 | `NRchiODE` | ∂²_χ W = F | NRH03 "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F" | V |
| 137 | `NRsource` | the source F | same check (F as displayed) | V |
| 138 | `NRg` | ρ(χ) | NRH03 "NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)" | V |
| 139 | `NRGprofile` | G(χ) normalization | NRH03 "NRGprofile: G'(0) = 0", "…(e^{s chi} − 1 − s chi)…" | V |
| 140 | — | q = √(L₊L₋/2) μ_RG^{−2}, χ, μ_RG = e^{y/l} | — (definitions; used in the next line) | D |
| 141 | — | μ_RG dχ/dμ_RG = −2√2 sinh(χ/√2) | NRH03 "SM: RG rapidity mu d chi/d mu = −2 Sqrt[2] Sinh[chi/Sqrt[2]]" | V |
| 142 | `SMBtransform` | finite B-transformation | NRH03 "SMBtransform, SMWshift…", "Omega_b is O(3,3)" | V |
| 143 | `SMWshift` | Ω_bH(W)Ω_bᵀ = H(W − 2b_{+−}) | NRH03 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| 144 | `SMdyg` | doubled-yet-gauged action | — (defining action; its reductions are verified below) | D |
| 145 | `SMphysicalsectionA` | Â_αμ, D_αx^M on the section | — | D |
| 146 | `SMRfirstorder` | Riemannian first-order form | NRH06 "the auxiliary equations give…", "eliminating the auxiliaries reproduces E_{mu nu}…" | V |
| 147 | `SMceff` | c_eff² = 2F | NRH06 "−det g_par = …", "c_eff^2 = 2F -> 4 Sqrt[L+L−…"; NRH02 "c_eff^2 := 2F…" | V |
| 148 | `SMlongstringE` | winding-string energy of the static probe (SM 5.1) | NRH06 "before : … d_t is Killing…", "before : gamma^{tau a} g_{t nu} d_a X^nu = … = 1", "before : B_{t phi} = …", "before : E = −Int d sigma P_t … reproduces E(y) = …", "−det g_(t,phi) = l^2 (e^{2y/l} − L+L− e^{−2y/l})^2 per winding…", "E(y) = …", "remark: with phi_0 = 0 the area density equals l e^{-2d}…" | V |
| 149 | `SNCtau` | SNC clock forms | NRH03 "SNCtau: tau+ = dx+ − L_− e^{−2y/l} dx− + O(e^{−4y/l})", "unit clock determinant", "W-part of the lower-right block = W(tau+ tau− + tau− tau+)"; NRH07 "x rows = Sqrt[2] tau^pm…" | V |
| 150 | `SMdygconstraints` | τ⁺·∂̄x = 0 = τ⁻·∂x | NRH06 "SM: H^{mu nu} tau^pm_nu = 0", "SM: the dual vectors Y, Ybar exist at every radius" (the kernel and multiplier structure) | V° |
| 151 | `SMdygGO` | exact reduced Lagrangian | NRH06 "symmetric-block route − antisymmetric-clock route = W (tau− . dx)(tau+ . dbar x)", "at chi -> 0 the W coupling reduces to…" | V° |
| 152 | `SMGO` | Gomis–Ooguri limit | NRH06 "at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)"; contraction bookkeeping | V |
| 153 | `SMvertex` | V_W = (1/4πα′)W∂x⁺∂̄x⁻ | NRH06 "with the 1/(2 pi alpha') prefactor this is V_W = …" | V |
| 154 | `SMWisB` | H(W) = Ω_bH(0)Ω_bᵀ, b_{+−} = −W/2 | NRH03 "SMWisB…" | V |
| 155 | `SMdeltaL` | antisymmetric clock coupling | NRH06 "symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)" | V |
| 156 | `SMFTsector` | S_y, S_FT | — (definitions; used in SM 158) | D |
| 157 | `SMWZWlevel` | k = l²/α′ | NRH06 "k = \|Int_{S^3} H\| / (4 pi^2 alpha') = l^2/alpha'" | V |
| 158 | `SMFTweight` | T_y, h_y(a) | NRH06 "the two OPE contributions assemble to h_y(a)/(z−w)^2" | V |
| 159 | `SMBRSTradial` | (∂²_y + (2/l)∂_y)f = 0 | NRH06 "the marginal roots…", "(d_y^2 + (2/l) d_y) f = 0…" | V |
| 160 | `SMBRSTmomentum` | h_y(p_y), roots p_y = 0, 2i/l | NRH06 "h_y(i p_y) = … P_y = p_y − i/l…" | V |
| 161 | `SMBRSTvertex` | U_{W₁} is a weight-(1,1) primary | SM 158–160 plus NRH06 "<x^+ x^−> = 0…" (W₁ has weight (0,0)) | V° |
| 162 | `SMBRSTcentral` | c_y = 1 + 6α′/l², total 3(k+2)/k | NRH06 "c_y/2 = …", "c_{beta gamma} + c_y = …" | V |
| 163 | `SMWgaugeobstruction` | W₀ gauge, e^{−2y/l}W₁ not | NRH06 "Lhat_xi H^infty = D H^infty + H^infty D^T…", "the conditions force…", "(db)_{+−y} = …", "v^y = … so closure forces d_y varpi = 0" | V |
| 164 | `SMBRSTfusion` | self-contraction, h_n, resonances | NRH06 "e^{a y(z)} e^{a y(0)} ~ \|z\|^{-alpha' a^2} = \|z\|^{-4 alpha'/l^2} for a = -2/l (from <y y> = -(alpha'/…" (3 checks) (three checks) | V |
| 165 | `SMupliftblocks` | H₁₀ = H₃ ⊕ H_{S³} ⊕ H_{R⁴}, d₁₀ = d₃ + d_{S³} | NRH07 `d10Assemble` and the symbolic ten-dimensional identities "S_(0)^{(10)} = 0", "(P S Pbar)^{(10)} = 0" for both uplifts (arbitrary chiral L±; one-sided hairy family with arbitrary W₀, W₁) | V |
| 166 | `SMDFTKilling` | generalized Killing equations | — (definition; `GenLieH`/`GenLieD` implement it: NRH07 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm", "Lhat_xi d = 0 (the radial component compensates the divergence)" (2 checks)) | D |
| 167 | `SMtypeIIKS` | type-II Killing-spinor systems | — (definition of the systems solved in NRH07) | D |
| 168 | `SMsusyclosure` | closure on ĥL_X, X = iε̄₂Γε₁ | NRH07 "C = i sigma_2 is the Majorana conjugation…", "X^M = …", "Lhat_X H^infty = 0 and Lhat_X d = 0…", "the bilinear is symmetric…" (three-dimensional calibration) | V° |
| 169 | `SMsemicov` | semi-covariant connection | NRH01 `GammaDFT`; NRH02 "Gamma compatibility: nabla_C P_AB = 0", "Gamma dilaton trace…", "Gamma torsionless…" | V |
| 170 | `SMuplift` | AdS₃×S³×R⁴ with flux | NRH07 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2", "R_{mu nu} = (1/4) H H on the S3 factor", "R(AdS3, Banados) = −6/l^2 and H^2 = −24/l^2…", "R_{mu nu} = (1/4) H H on the AdS3 factor", "S_(0)(S3 block) = +4/l^2", "S_(0)(R4 block) = 0" | V |
| 171 | `SMRDFTKilling` | Ordinary-field form of the generalized Killing equations | Cited decomposition; the ordinary Killing vectors are not enumerated in the suite. | S |
| 172 | `SMRlocaliso` | k^{(ij)} = s_is_j stabilizers | NRH07 "k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' − (l^2/4) k''' = 0" | V |
| 173 | — | vol₃, H₃ = −(2/l)vol₃, H_{y−+} = +(2/l)√\|g₃\| | NRH07 "flux orientation: H_{y−+} = +(2/l) Sqrt[\|g_3\|…" | V |
| 174 | — | ∇_με± = ±(1/2l)γ_με± | — (the torsionful Killing-spinor equation of the Riemannian branch; its consequence SM 175–176 is checked) | S |
| 175 | `SMRcomponentHill` | first-order pair → Hill equation | NRH07 "the first-order pair … closes into (l^2/2) s'' = L s" | V |
| 176 | `SMRlocalKS` | Hill equations for s± | NRH07 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s", "global AdS3 (L = −1/4)…", "massless BTZ (L = 0): the (u, v) pair shifts by exactly 2 pi v over one circuit (unipotent)", "constant L > 0: monodromy multipliers e^{pm 2 Sqrt[L0] pi} are real and not +-1: (0,0) kernel" | V |
| 177 | — | N_local = 16 = 8 + 8 | Standard local AdS3 x S3 result; arithmetic bookkeeping alone is not an independent Killing-spinor derivation. | S |
| 178 | `SMvielbein` | exact non-Riemannian double vielbein | NRH07 "V_M^p V_Np = P_MN and Vbar_M^pbar Vbar_Npbar = Pbar_MN for the exact frame , any hh, sigma, W", "V^M_p Vbar_{M qbar} = 0 and P + Pbar = J", "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)", "at hh = 0 and W_0 = 0, lowering the local indices with eta, etabar gives the limiting frame", "x rows = Sqrt[2] tau^pm (), x~ rows = the dual vectors -Y/Sqrt[2], Ybar/Sqrt[2], and the W entries =…" | V |
| 179 | `SMvielbeincheck` | defining relations | same checks | V |
| 180 | `SMNRlocalstabilizer` | local stabilizer system | NRH07 "eps = c/Sqrt[L] solves…", "exact one-sided identity Lhat_xi H = delta H, with weights (1,2)"; NRH03 "line 2: … delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i" | V |
| 181 | `SMexactiso` | vacuum isometries | NRH07 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm"; "X^M = (omega_+, 0, -(l/2) v'; v^+, 0, -(l/2) v') with v^+ = 2 f1 f2, omega_+ = -2 l^2 f1' f2'" | V |
| 182 | `SMweighteddilaton` | weighted dilaton condition | NRH07 "Lhat_xi d = 0…", "equivalently d_M(e^{-2d} xi^M) = 0" | V |
| 183 | `SMkillingspinor` | vacuum Killing spinor E = (√2f, l∂₊f) | NRH07 "SM: vacuum spin connection = displayed…", "D_{pbar} E = 0…", "gamma^p D_p E = E/(Sqrt[2] l)" | V |
| 184 | `SMcomplexblocks` | σ_i, τ_i, ρ_m | NRH07 (matrices used in SM 185–187 checks); "rep: {gamma^p, gamma^q} = 2 eta^{pq}" | V |
| 185 | `SMgammaten` | ten-dimensional Γ^p̂, Γ₁₁ | NRH07 "{Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32", "Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0" | V |
| 186 | `SMgammaMajorana` | Majorana intertwiner | NRH07 "BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, …" | V |
| 187 | `SMgammabarred` | barred Clifford algebra | NRH07 "{Gammabar, Gammabar} = −2 eta, Gammabar_{pq} = −Gamma_{pq}, same Majorana intertwiner" | V |
| 188 | `SMreducedDirac` | reduced two-component system | NRH07 "reduced system…", "the opposite channel E = (0, g(x+)) has eigenvalue −1/(Sqrt[2] l)", "both channels also satisfy D_{pbar} E = 0" | V |
| 189 | `SMinternalprojectors` | Complex product spinors with auxiliary/R4 signs (-,-) and (+,+) | NRH07 checks Clifford identities and S3 integrability, but its historical projection-rank check does not verify this updated product basis or the following real Majorana combinations. | V° |
| 190 | `SMspinorcountchain` | Majorana, Weyl and background-projector count | NRH07 verifies 32 real Majorana and 16 real Weyl components. Its historical complex-rank count is not a verification of the final four-dimensional real kernel. | V° |
| 191 | `SMhairyKS` | Four arbitrary real chiral functions on each one-sided NR branch | NRH07 checks the displayed reduced jet system; full D=10 fermionic equations with the current internal basis and nonzero supercharges are not verified by that check. | V° |
| 192 | `SMcandidateD` | gauge derivative 𝔻 | — (definition; used in NRH08) | D |
| 193 | `SMcandidateaction` | candidate action and its reduction | NRH08 "the P^(0) term contains only D_+…", "each term selects a single independent component of phi_A", "gamma^oplus … rank one" | V |
| 194 | `SMcandidategauge` | internal gauge invariance | NRH08 "delta_Lambda L = 0 for the non-abelian bosonic sector" | V |
| 195 | `SMcandidateWitt` | Witt transformations | NRH08 "delta_v L = d_+(v^+ L) + d_−(v^− L)"; "Lhat_xi H^(0) = 0 <=> d_− v^+ = 0 = d_+ v^−" | V |
| 196 | `SMcandidateWittL` | δ_vL = ∂(vL) | same | V |
| 197 | `SMcandidatefermionic` | chiral fermionic transformations | NRH08 "delta_{eps+} L = …", "delta_{eps−} L = …", "the fermionic parameters carry arbitrary chiral profiles" | V |
| 198 | `SMcandidatefermionicL` | their total-derivative variations | same | V |
| 199 | `SMcandidateextra` | Grassmann-even chiral transformation | NRH08 "delta_zeta L = …" | V |
| 200 | `SMcandidateextraL` | its variation | same | V |
| 201 | `SMcandidatetrivial` | equation-of-motion redundancies | NRH08 "the transformations are built from the fermion equations of motion…" | V |
| 202 | `SMcandidatetrivialL` | δ_αL total derivative | NRH08 "delta_alpha L = d_−(psi^+ delta_alpha psi^+) + d_+(psi^− delta_alpha psi^−) identically (off shell)", "the two-component fermions used here are genuinely Grassmann…" | V |

## Additional scope

- NRH05 verifies the bosonic covariant charge and its centerless NR algebra. It
  does not verify fermionic charges. The common charge is now inline in the Letter.
- NRH04 checks the Gamma-squared variation on constant BTZ with arbitrary
  tangential h and delta d in radial gauge; radial h components are not varied.
- The on-shell action value with sqrt(mu) outside the boundary integral assumes
  constant positive L+ and L-, and no additional interior boundary functional.
- The massless-scalar falloff analogy does not identify a boundary CFT dimension
  with the worldsheet (1,1) weight. The doubled source pairing remains crossed.
- NRH07 checks ten-dimensional bosonic equations for arbitrary chiral R data and
  the one-sided NR family. It does not establish the full two-sided fermionic
  classification, a global interior completion, or alpha-prime exactness of NR.
- The scalar response and the displayed zeros in the particular 5x5 matrix must
  be read with the remainder and contact terms specified in the manuscript.
