# Equation ledger

Current source SHA-256: `57960AAD42F2A9150C0958D1C253F08B0939B76C0E646EE250A7D824D6319B68` (2026-09-30).

Letter (1)–(21); SM (1)–(129); **153 numbered displays**.
Status counts: V 77, V° 45, D 22, S 5, N 4.

This is a static coverage audit of the existing public Wolfram calculations against the
revised manuscript. **It is not a new execution of the suite.** The recorded 336/336 run
belongs to the 2026-09-17 source; see [REFERENCE_RUN.md](REFERENCE_RUN.md).

**V**: an existing symbolic check of the relation at the stated branch/order scope;
**V°**: only components, a specialization, earlier equivalent ingredients or normalization
are checked, without a direct check of the full current display; **D**: definition or
prescription; **S**: cited/stated identity not checked here; **N**: no relevant direct public check.
A reused LaTeX label does not by itself transfer verification to a changed formula.

In particular, the generic off-shell Codazzi relation, the compact metric first variation,
the full new metric Hessian, and the K3/interior-response chain are not established by
the old boundary Ward, coset-coordinate or conditional-normalization tests.
Existing ordered second-variation calculations supply particular asymptotic kernels;
they do not determine the full state-dependent Green function or the hair self-response.
Tensor Box curvature and connection terms are retained in the manuscript, but their
complete covariance is cited, not newly machine-verified by this archive.

File prefixes such as NRH05 refer to the `.wl` sources. Quoted check identifiers remain
searchable even when their earlier displayed equations were compressed into prose.
The detailed restrictions in the coverage column take precedence over a status letter.

## Letter

| No. | Label | Content | Public coverage and limitations | Status |
|---:|---|---|---|:---:|
| 1 | `Rfields` | NS–NS Bañados family (frame ϑ±, metric, B, dilaton) | NRH02 "Rfields: dy^2 - 2 theta+ theta- reproduces the displayed metric components", "Rfields: B_{-+} = e^{2y/l} + e^{-2y/l} L+ L- (B = (...) dx^- ^ dx^+)" | V |
| 2 | `RDFTfields` | DFT variables H_MN and e^{-2d} for the Riemannian family | NRH02 "RDFTfields: H J H = J (O(3,3) constraint)", "RDFTfields: e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l}) = sqrt(-g) at phi_0 = 0" | V |
| 3 | `Rboundary` | asymptotic falloffs | NRH02 "Rboundary: H - H^infty = O(e^{-2y/l}) and d + y/l = O(e^{-4y/l})" | V |
| 4a | `boundaryH` | H^∞, d^∞ | NRH02 "boundaryH: H^infty is the displayed 6x6 matrix", "boundaryH: the induced boundary H^(0) is of type (1,1): vanishing upper-left block" | V |
| 4b | `Rboundaryframe` | aligned D = 2 boundary frame | NRH02 "Rboundaryframe: the aligned D = 2 double vielbein reproduces P^(0), Pbar^(0)" | V |
| 5 | `NRvariables` | (Π, q, e^σ, χ) and e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) | NRH02 "NRvariables: chi = 2 Sqrt[2] arctanh q and e^{pm sigma} sinh(chi/2) = e^{-2y/l} L_pm + O(e^{-6y/l})" | V |
| 6 | `NRHcompact` | the everywhere non-Riemannian H_MN | NRH02 "NRHcompact: type (1,1) at every radius (upper-left block diag(0,0,1))", "NRHcompact: H J H = J for the exact non-Riemannian matrix (generic W)" | V |
| 7 | `NRdilaton` | e^{−2d} = e^{2y/l}(1 − q²), d = −y/l + ln cosh(χ/2√2) | NRH02 "NRdilaton: e^{-2d} = e^{2y/l}(1 - q^2) with q = tanh(chi/(2 Sqrt[2]))" | V |
| 8 | `NRWgeneral` | exact hair profile W | NRH02/NRH04: exact radial equation and homogeneous modes, using GG for calligraphic I; NRH05: order-u^2 expansion of the hair profile. | V |
| 9 | `MainGprofile` | Radial profile calligraphic I(chi), its integral normalization and initial conditions | NRH02: derivative limit, derivative series and I_radial second derivative = rho. The value I_radial(0)=0 is imposed by the integral definition. Executable GG/Gp names are retained; see README. | V° |
| 10 | `RGKPW` | GKPW relation Z_DFT = Z_CFT | — (defining relation) | D |
| 11 | `Rrenvariation` | On-shell boundary variation of the renormalized action | NRH06 "SMrenvariation, SMvolumecounterterm: on shell the bulk term cancels delta(-2 Lambda e^{-2d}) and the counte…" | V° |
| 12 | `Mainmomenta` | Generalized-metric and scalar radial momenta A and B | NRH06 "SMmetricradialdictionary, SMscalarradialdictionary: coefficient matching -2K = (16 pi G)^{-1} A^y and 2 T_(…", "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" (+1 more) | V° |
| 13 | `Codazzi` | Off-shell DFT analogue of the contracted Codazzi identity in radial gauge | No generic off-shell intrinsic-divergence versus bulk-curvature check. Exact-saddle and constant-boundary Ward tests do not establish this identity. | N |
| 14a | `RKdef` | ⟨K⟩, ⟨T₍₀₎⟩ from the rescaled momenta | NRH02 "RKdef: -2K = (16 pi G)^{-1} A^y and 2 T_(0) = (16 pi G)^{-1} 2 (B^y + 4/l) fix the coefficients -1/(32 pi G…" | V |
| 14b | `RDFTconservation` | On-shell conservation of the boundary DFT energy-momentum tensor | NRH02: constant limiting-frame divergence components and flat-boundary Ward equations only; no derivation from the generic off-shell Codazzi identity. | V° |
| 15 | `Rkilling` | common R/NR asymptotic-symmetry generator (ξ^±_R with the l²e^{−2y/l}∂²ε/4 tail, ξ^±_NR = ε^±; common radial and dual components) | NRH02 "Rkilling, RVirasoro: Lhat_xi H - delta_eps H = O(e^{-4y/l}) with delta L = eps L' + 2 L eps' - (l^2/4) eps'''", "Rkilling: Lhat_xi d = O(e^{-4y/l}) (the radial component compensates the weight)" (+1 more) | V |
| 16 | `RVirasoro` | δ_εL± with −(l²/4)∂³ε | NRH02 "Rkilling, RVirasoro: Lhat_xi H - delta_eps H = O(e^{-4y/l}) with delta L = eps L' + 2 L eps' - (l^2/4) eps'''", "charges: with T = 8 pi K = L+/(2 G l), RVirasoro is delta T = eps T' + 2 T eps' - (c/12) eps''' with c = 3l…" | V |
| 17 | `NRasympt` | δ_εL±, δ_εW₁ at W₀ = 0 | NRH02 "Rkilling (NR form), NRasympt: Lhat_xi H - delta_(L, W1) H = O(u^-2) componentwise", "NRasympt: delta psi is equivalent to delta L = eps L' + 2 L eps' (no central term)" (+1 more) | V |
| 18 | `Mainonepoints` | R/NR one-point matrices and vanishing scalar responses | NRH02 "Mainonepoints on R: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) (exact fam…", "Mainonepoints on R: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0" (+2 more) | V |
| 19 | `Mainfivebyfive` | Connected five-by-five two-point response matrix with interior and contact contributions | NRH02/NRH06: particular stress and mixed kernels and the relation M_+^R = L_- A_+^R. The full response matrix, undetermined interior functions and contact terms are not fixed. | V° |
| 20 | `Rcorrelators` | particular A_NR, A_R, M_NR, M_R kernels, with coefficients at x | NRH02 "Rcorrelators: A_+^NR = [L+'/D+ - 2 L+/D+^2]/(128 pi^2 G l)", "Rcorrelators: A_+^R = A_+^NR + 3 l/(256 pi^2 G D+^4)" (+1 more); NRH06 "SMgeneralpositionkernels, Rcorrelators: A_+^NR = [L+'/D+ - 2L+/D+^2]/(128 pi^2 G l)", "SMgeneralpositionkernels, Rcorrelators: A_+^R = A_+^NR + 3l/(256 pi^2 G D+^4)" (+1 more) | V |
| 21 | `Mainworldsheet` | ℒ_GO = (1/2πα′)(β∂̄x⁺ + β̄∂x⁻ + ∂y∂̄y), V_𝒲 = (𝒲/4πα′)∂x⁺∂̄x⁻ | NRH02 "Mainworldsheet: eliminating beta, betabar from the first-order action reproduces E_{mu nu} dx dbar x (E = g…", "Mainworldsheet: V_W = (1/(4 pi alpha')) W dx+ dbar x- with the 1/(2 pi alpha') prefactor of the reduced action" | V |

## Supplemental Material

| No. | Label | Content | Public coverage and limitations | Status |
|---:|---|---|---|:---:|
| SM 1 | `SMgamma2` | Gamma-squared action identity | NRH03: e^(-2d) S_(0) equals the Gamma-squared density plus the divergence of e^(-2d) B^M on both saddles; the flux W-independence is also checked. Tangential integration uses the stated periodic or decay conditions. | V |
| SM 2 | `SMgamma2flux` | B^y = 4∂_y d, e^{−2d}B^y = −(4/l)(e^{2y/l} + μe^{−2y/l}) | NRH03 "SMgamma2flux on R: B^y = 4 d_y d and e^{-2d} B^y = -(4/l)(e^{2y/l} + L+L- e^{-2y/l})", "SMgamma2flux, SMmudefinition: the same formula e^{-2d} B^y = -(4/l)(e^{2y/l} + mu e^{-2y/l}) holds on both …" (+1 more) | V |
| SM 3 | `SMmudefinition` | μ dichotomy | NRH03 "SMmudefinition on NR: e^{-2d} = e^{2y/l}(1 - q^2) = u - (L+L-/2)/u with q = Sqrt[L+L-/2] e^{-2y/l}, so -2 d…", "SMgamma2flux, SMmudefinition: the same formula e^{-2d} B^y = -(4/l)(e^{2y/l} + mu e^{-2y/l}) holds on both …" | V |
| SM 4 | `SMgamma2cutoff` | S_ren(Y) | NRH03 "SMgamma2cutoff: the regulated flux plus the volume counterterm equals (16 pi G)^{-1}[(8 mu/l) e^{-2Y/l} - (…", "SMgamma2cutoff: the interior endpoint e^{2y*/l} = Sqrt[mu] contributes (4/l)(Sqrt[mu] + Sqrt[mu]) = (8/l) S…" | V |
| SM 5 | `SMgamma2value` | S_ren = −8√μ/(16πGl)∫d²x | NRH03 "SMgamma2value: Y -> Infinity gives S_ren = -(8 Sqrt[mu])/(16 pi G l) Int d^2x" | V |
| SM 6 | — | Killing horizon of the Bañados family, e^{4y/l} = L₊L₋ (where e^{−2d} = 0) | NRH04 "SM2.1: exterior terminates at the Killing horizon e^{4y/l} = L+ L- where e^{-2d} = 0" | V |
| SM 7 | `NRhill` | ψ±, A± = ¼(∂lnL)² − ½∂²lnL = ψ″/ψ | NRH04 "NRhill: A = (1/4)(d ln L)^2 - (1/2) d^2 ln L = psi''/psi with psi = L^{-1/2}" | V |
| SM 8 | `NRradialchange` | q, χ, ∂_y, ∂_q | NRH04 "NRradialchange: d chi/d q = 2 Sqrt[2]/(1 - q^2) and d_y = -(2q/l) d_q for q = e^{-2y/l} Sqrt[Pi/2]" | V |
| SM 9 | `NRradialoperator` | radial operator identity | NRH04 "NRradialoperator: (1-q^2)^2/8 (d_q^2 - 2q/(1-q^2) d_q) = d_chi^2" | V |
| SM 10 | `NRchiODE` | ∂²_χ W = F | NRH04 "NRchiODE: the tensor EDFE contains d^2 W/d chi^2", "NRchiODE, NRsource: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F with the displayed source" (+1 more) | V |
| SM 11 | `NRsource` | the source F | NRH04 "NRchiODE, NRsource: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F with the displayed source" | V |
| SM 12 | `NRg` | Universal radial kernel rho(chi) | NRH04: NRg/NRGprofile second-derivative identity for the radial calligraphic I profile. | V |
| SM 13 | `NRGprofile` | Normalized radial profile calligraphic I and elementary exponential particular solution | NRH04: derivative identity, derivative limit/series and exponential second derivatives. The zero integration constant is specified by the integral definition. | V° |
| SM 14 | — | constant-L radial variables q = √(L₊L₋/2) μ_RG^{−2}, χ = 2√2 arctanh q, μ_RG = e^{y/l} | definition; used by the NRH04 check of SM15 | D |
| SM 15 | — | μ_RG dχ/dμ_RG = −2√2 sinh(χ/√2) | NRH04 "constant L: q = Sqrt[L+L-/2] mu_RG^{-2}; mu d chi/d mu = -2 Sqrt[2] sinh(chi/Sqrt[2])" | V |
| SM 16 | `SMBtransform` | finite B-transformation | NRH04 "SMBtransform: Omega_b is O(3,3)", "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| SM 17 | `SMWshift` | Ω_bH(W)Ω_bᵀ = H(W − 2b_{+−}) | NRH04 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched"; NRH05 "SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly" | V |
| SM 18 | `SMcosetreconstruction` | Mixed fluctuation and reconstruction | NRH05: `SMcosetreconstruction`, `SMmixedfluctuation`, `constraint: delta(H J H) = 0` Scope: Algebra is checked with the limiting frames and generic mixed h; not a newly added differential identity. | V |
| SM 19 | `SMframevariation` | Lorentz-gauge frame variation | NRH05: `SMframevariation` Scope: Same unchanged variation algebra. | V |
| SM 20 | `SMinfinityvielbein` | Limiting bulk double vielbeins | NRH05: `SMflatmetrics, SMinfinityvielbein`, `both saddle frames tend to the limiting frame SMinfinityvielbein at u = 0` Scope: Exact boundary algebra and asymptotic saddle limit. | V |
| SM 21 | `SMbackgroundconnection` | Torsionless DFT connection | NRH05: `SMbackgroundconnection on R/NR: nabla_C P_AB = 0 through z^2`, `Gamma^B_BA = -2 d_A d through z^2`, `torsionless cyclic sum` Scope: The helper implements this connection. Its properties are checked on the two saddle series through z^2, not a general proof of the displayed formula. | V° |
| SM 22 | `SMfullscalarvariation` | Full off-shell scalar curvature variation | No direct public check. Scope: Cited Lee-Park formula; NRH05 directly linearizes exact curvatures instead of checking this universal Box identity. | S |
| SM 23 | `SMfullcurvaturevariation` | Full off-shell projected curvature variation | No direct public check. Scope: Cited Lee-Park formula; no public universal Box implementation/check of the full off-shell variation. | S |
| SM 24 | `SMfullOmega` | Omega divergence definition and flat-frame expression | No direct public check. Scope: No direct public check of equality between the two displayed forms. | D |
| SM 25 | `SMboxdefinition` | Box = Delta - barDelta and section condition identities | No direct public check. Scope: Declared universal-operator definitions/identities; public direct component linearization does not verify them. | D |
| SM 26 | `SMmixedboxdefinition` | Tensor Box with curvature and connection terms | No direct public check. Scope: No direct public implementation/comparison of the displayed mixed Delta formula; the new prose is not tested by component EDFE checks. | D |
| SM 27 | — | Scalar Box as dilaton-weighted divergence | No direct public check. Scope: Not separately checked; do not confuse scalar Laplacian with the full tensor operator. | D |
| SM 28 | `SMboxcurvature` | Connection field strength and Ricci contraction | No direct public check. Scope: Definition; curvature helpers use these operations but this is not a separate independent theorem check. | D |
| SM 29 | `GammaDFT` | Gamma-squared density and curvature/divergence identity | NRH06: `variation, SMBdefinition, SMAdefinition` Scope: NRH06 uses Gamma2Density and verifies its first variation on constant BTZ. Full density identity is additionally checked on saddle families by NRH03 SMgamma2, not generically here. | V° |
| SM 30 | `variation` | First variation of Gamma-squared action | NRH06: `variation, SMBdefinition, SMAdefinition: ... generic tangential h and delta d on BTZ` Scope: Arbitrary tangential fluctuations, but only the constant-BTZ R background; not a universal off-shell-background check. | V° |
| SM 31a | `SMBdefinition` | Definition of B vector | NRH06: `variation, SMBdefinition, SMAdefinition` Scope: Implemented as GammaBVector and used in checked variation; the two defining expressions are not separately proven for arbitrary backgrounds here. | D |
| SM 31b | `SMAdefinition` | Definition of generalized metric momentum | NRH06: `variation, SMBdefinition, SMAdefinition` Scope: Implemented explicitly and used in the BTZ variation check; definition rather than a general identity proof. | D |
| SM 32 | `SMvolumecounterterm` | Volume counterterm | NRH06: `SMrenvariation, SMvolumecounterterm` Scope: Coefficient shift B^y -> B^y+4/l and bulk cosmological cancellation checked; not nonlinear derivative-counterterm completeness. | V° |
| SM 33 | `SMrenvariation` | On-shell radial first variation in h,A,B | NRH06: `variation, SMBdefinition, SMAdefinition`, `SMrenvariation, SMvolumecounterterm` Scope: Supported by BTZ first variation plus algebraic cosmological/counterterm shift, not a general exact radial proof. | V° |
| SM 34 | `SMradialmetricvariation` | Exact radial first variation in delta H and delta d | No direct public check. Scope: New contraction h A^y = (1/4) delta H^{AB} partial_y H_AB is not directly checked by NRH05/06/07. Old asymptotic momentum checks are insufficient for this exact finite-cutoff formula. | N |
| SM 35 | `SMsecondvariation` | Polarized metric/dilaton Hessian | NRH06: `SMsecondvariation: mixed coset constraint at second order`, `SMsecondvariation: h_1 projection and delta_2 h_1 = 0`, `SMdirectcutoffbilinear on R/NR` Scope: The label survived a substantive change. Existing tests cover coset-coordinate auxiliaries and the OLD ordered asymptotic bilinear through Y^0. They do not directly test the NEW symmetric metric formula, its 1/(32 pi G) factor, or interior-flux symmetry assumptions. | V° |
| SM 36 | `SMderivativectNR` | Quadratic NR derivative counterterm | NRH06: `SMderivativect on NR`, `negative control on NR: without the counterterm the divergences remain`, `SMctlocality on NR` Scope: Displayed ctQ matches; divergence cancellation is modulo tangential total derivatives with stress constraints in both slots. No claim of universal nonlinear intrinsic completion. | V |
| SM 37 | `SMderivativectR` | Quadratic R derivative counterterm | NRH06: `SMderivativect on R`, `negative control on R: without the counterterm the divergences remain`, `SMctlocality on R` Scope: Displayed R additions match ctQ. Same quadratic/asymptotic/IBP restrictions as SM36. | V |
| SM 38 | `SMEinsteintensor` | G = 4(PSbarP)_[MN] - J S0/2 | NRH05: `SMexactconstraintpropagation: delta G_{y~ N}` Scope: Old dG is ASSIGNED using this definition before testing its projected row. That check does not independently derive the Einstein tensor or prove fixed-Lambda equivalence. | D |
| SM 39 | `SMexactboxEDFE` | On-saddle linearized EDFE in Box form | NRH05: `SMexactnormalform`, `SMexactCauchyconstraints`, `SMEzerocomponents`, `SMEoneNRcomponents`, `SMRoperatorrule, SMDeltaEonecomponents` Scope: Public code directly linearizes scalar/projected curvature through u and checks component operators. It does not implement and independently compare the fully covariant Box expression or its full Omega bridge. | V° |
| SM 40 | `SMexactRdata` | Exact R fields and lower-index frames | NRH05: `SMbackgroundexpansion`, `SMexactRdata` Scope: Exact frame/projector relations and background dilaton match unchanged formulas. | V |
| SM 41 | `SMRvacuumdilatoncoefficient` | Vacuum dilaton u coefficient | NRH05: `SMRh2L_*, SMRdd2L, SMRdd2` Scope: L+=L-=0 specialization of the directly checked old dd2target[R]. No separate current-label vacuum prescription test. | V° |
| SM 42 | `SMRvacuumsolution` | Explicit five-field R vacuum expansion | NRH05: `u^0 order`, `SMleadingdilaton, SMleadingstress, SMNRsol`, `SMRconstraint0, SMRconstraint1`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Old generic coefficients and constraints support the local expansion. Current integrated vacuum responses and chosen zero modes are not transcribed/substituted as a complete current-display test. | V° |
| SM 43 | `SMRvacuumregularscalar` | Fixed-flux physical dilaton equation and regular K3 solution | No direct public check. Scope: No K3/Bessel fixed-flux bulk scalar check exists in NRH05/06/07. | N |
| SM 44 | `SMRvacuumresponse` | Interior-determined H_R logarithmic response and C_mu | No direct public check. Scope: No current public K3 boundary matching or mu cancellation test. H_R is free in the old general-source suite. | N |
| SM 45 | `SMRvacuumhessian` | Correlator = (1/4i) source Hessian | NRH06: `SMstressparticular`, `SMNRhairhessian`, `SMonepointequivalence` Scope: Quarter-factor and Lorentzian normalization algebra checked; no full new vacuum interior Hessian computation. | V° |
| SM 46 | `SMLorentziandeltaconvention` | Joint Lorentzian delta prescription | NRH06: `SMLorentziandeltaconvention`, `SMPBHkernel` Scope: Box-integral normalization and inverse-derivative kernel away from light cone checked; not a general distribution-theory proof. | V° |
| SM 47 | `SMRtwopt` | R vacuum stress correlator | NRH06: `SMgeneralpositionkernels, Rcorrelators`, `SMRtwopt` Scope: Plus kernel and vacuum central normalization directly checked. Minus sector follows by relabeling plus/minus; not separately executed there. | V |
| SM 48 | `SMRvacuumhaircorrelator` | R vacuum hair self-correlator | NRH06: `SMRhairconditional` Scope: Only algebraic normalization/re-expression of an OLD conditional completion is tested. It does not establish the new unconditional specified-vacuum result from K3 matching or source Hessian. | V° |
| SM 49 | `SMradialansatz` | General near-boundary polynomial expansion | NRH05: `u^0 order`, `u^1 order on R/NR` Scope: Ansatz definition; the old hierarchy uses equivalent h^(1), h^(1b), h^(2), h^(2L), h^(2LL) naming. | D |
| SM 50 | `SMRstresspluscoefficients` | R plus stress coefficients in compact a,b,omega basis | NRH05: `SMleadingdilaton, SMleadingstress, SMNRsol`, `SMRconstraint0, SMRconstraint1`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Leading/log coefficients and unintegrated stress constraint checked. New explicit b++/omega+/f0 expression is not in the old code. | V° |
| SM 51 | `SMRstressminuscoefficients` | R minus stress coefficients in compact a,b,omega basis | NRH05: `SMleadingdilaton, SMleadingstress, SMNRsol`, `SMRconstraint0, SMRconstraint1`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Leading/log coefficients and unintegrated stress constraint checked. New explicit b--/omega-/f0 expression is not in the old code. | V° |
| SM 52 | `SMRtypecoefficients` | R type-channel coefficients | NRH05: `SMtypeconstraint on R`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Old c_R parameter equals current -4 f0. The two displayed coefficients match after this trivial renaming, but current display is not directly named/tested. | V° |
| SM 53 | `SMRhaircoefficients` | R hair coefficients after stress elimination | NRH05: `u^0 order`, `u^1 order on R: type-changing response H_s is unconstrained`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Old h2Ltarget[R,pm] contains stress-response derivatives; new source-only compact expression requires substitution of stress constraints, not presently tested as a bridge. | V° |
| SM 54 | `SMRdilatoncoefficients` | R dilaton coefficients | NRH05: `SMleadingdilaton, SMleadingstress, SMNRsol`, `SMRh2L_*, SMRdd2L, SMRdd2` Scope: Term-by-term equivalent to old verified coefficient expressions, but the updated (n,k) notation is not directly transcribed in current checks. | V° |
| SM 55 | `SMRintegratedstress` | Integrated R stress constraints and a,b,omega definitions | NRH05: `SMRconstraint0, SMRconstraint1` Scope: Old code checks differential stress constraints. New inverse-derivative integration, a/b/omega parametrization, and chiral zero-mode bookkeeping are not directly tested. | V° |
| SM 56 | `SMNRcompactsolution` | Compact NR five-field expansion | NRH05: `SMleadingdilaton, SMleadingstress, SMNRsol`, `SMtypeconstraint on NR`, `SMNRh2L_*, SMNRdd2L, SMNRdd2` Scope: Equivalent ingredients of the old general hierarchy checked; current compact U form and integrated responses have no direct full-expression bridge test. | V° |
| SM 57 | `SMNRintegratedstress` | Integrated NR stress responses including W1 | NRH05: `SMNRconstraint0, SMNRconstraint1` Scope: Differential equations checked with correct W1 coefficients. Current inverse-derivative solution and chiral zero modes are not directly substituted/tested. | V° |
| SM 58 | `SMNRcompactU` | Compact NR U log coefficient | NRH05: `SMNRh2L_*, SMNRdd2L, SMNRdd2` Scope: Expansion matches the old h2Ltarget[NR,pm], including -l W1 d+d-r/8. No direct equality check of the current compact U grouping. | V° |
| SM 59 | `SMstatefluctuations` | NR state variations | NRH06: `SMstatefluctuations on NR` Scope: Direct projection of exact state variation through u, including delta W1/2 and delta d=O(u^2). | V |
| SM 60 | `SMgeneralpositionkernels` | Particular separated-point A/M kernels | NRH06: `SMgeneralpositionkernels, Rcorrelators`, `SMgeneralparticularkernels` Scope: Particular plus-sector Ward kernels checked with mixed NR factor 1/(512 pi^2 G l). Minus formulas follow by relabeling. Does not fix homogeneous/interior-dependent or contact terms. | V |
| SM 61 | `SMNRhairhessian` | Full hair response via functional derivative of H_s | NRH06: `SMdirectstressrows, SMdirectremainingrows`, `SMNRhairhessian` Scope: Coefficient 1/(64 pi i G l) and finite action-response structure checked; functional H_s[j] remains unspecified, so no general hair-hair kernel is computed. | V° |
| SM 62 | `SMCPSform` | Field-dependent Noether-Wald surface-charge one-form | NRH07: `lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary` Scope: Cited prescription is implemented, not derived. Code subtracts K_delta_xi and includes improved Theta. | D |
| SM 63 | `SMRpotentialcomponents` | R limiting potential components | NRH07: `lim e^{-2d} Khat^{-y}[eps+]`, `lim e^{-2d} Khat^{y+}[eps-]` Scope: Both chiralities checked with equivalent reversed-index sign for y+. | V |
| SM 64 | `SMCPSresult` | NR limiting surface-charge components | NRH07: `k^{-y}[eps+] = (4/l) eps+ dL+`, `k^{+y}[eps-] = (4/l) eps- dL-`, `W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise` Scope: Boundary limits computed using sufficient asymptotic saddle order; no finite-radius charge claim. | V |
| SM 65 | `SMNRcharge` | Integrated charge and normalization | NRH07: `(i) k^{-y}[eps+] = delta[(4/l) eps+ L+]`, `(i) mirror`, `(iii) normalization chain` Scope: Integrability and 1/(4 pi G l) normalization checked; background subtraction is a choice of integration constant. | V |
| SM 66 | `SMdyg` | doubled-yet-gauged action | — (defining action; its reductions are verified below) | D |
| SM 67 | `SMphysicalsectionA` | Â_αμ, D_αx^M on the section | — | D |
| SM 68 | `SMRfirstorder` | Riemannian first-order form | NRH08 "the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+", "eliminating the auxiliaries reproduces E_{mu nu} dx^mu dbar x^nu, E = g - B" | V |
| SM 69 | `SMceff` | c_eff² = 2F | NRH02 "c_eff^2 := 2F = 2(e^{2y/l} + L+ L- e^{-2y/l}) ~ 2 e^{2y/l}: the Gomis-Ooguri limit is reached without tuning"; NRH08 "c_eff^2 = 2F -> 4 Sqrt[L+L-] at the horizon u = Sqrt[L+L-]" | V |
| SM 70 | — | string momentum density P_μ and the long-string energy E = −Q_{∂t} = (1/2πα′)∫dσ(√−γ − wB_{tφ}) | NRH08 "before : gamma^{tau a} g_{t nu} d_a X^nu = gamma^{tau a} gamma_{a tau} = 1 on the static embedding", "before : B_{t phi} = l (e^{2y/l} + L+L- e^{-2y/l}) and B_{t nu} X'^nu = w B_{t phi}" (+1 more) | V |
| SM 71 | `SMlongstringE` | Winding-string energy of the static constant-L probe (SM4.1) | NRH08: static constant-L probe with the specified winding orientation; not a general time-dependent background energy calculation. | V |
| SM 72 | `SNCtau` | SNC clock forms | NRH08 "SM: the dual vectors Y, Ybar exist at every radius (unit clock determinant)"; NRH09 "x rows = Sqrt[2] tau^pm, x~ rows = the dual vectors -Y/Sqrt[2], Ybar/Sqrt[2], and the W entries = -(W/(2 Sq…" | V |
| SM 73 | `SMdygconstraints` | Worldsheet constraints from the auxiliary field | NRH08: null-clock and dual-vector identities supporting the constraints; not a complete variation and elimination of all auxiliary components. | V° |
| SM 74 | `SMdygGO` | Finite-radius first-order non-Riemannian worldsheet form | NRH08: chi-to-zero W coupling and constraint-term algebra; the full finite-radius auxiliary-field elimination is not independently reproduced by that predicate. | V° |
| SM 75 | `SMGO` | Gomis–Ooguri limit | NRH08 "at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)" | V |
| SM 76 | `SMvertex` | V_W = (1/4πα′)W∂x⁺∂̄x⁻ | NRH08 "with the 1/(2 pi alpha') prefactor this is V_W = (1/(4 pi alpha')) W dx+ dbar x-" | V |
| SM 77 | `SMWisB` | H(W) = Ω_bH(0)Ω_bᵀ, b_{+−} = −W/2 | NRH04 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| SM 78 | `SMdeltaL` | antisymmetric clock coupling | NRH08 "symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)" | V |
| SM 79 | `SMFTsector` | Radial kinetic and Fradkin-Tseytlin sectors | Definitions; used in the radial BRST condition, SM82. | D |
| SM 80 | `SMWZWlevel` | k = l²/α′ | NRH08 "k = \|Int_{S^3} H\| / (4 pi^2 alpha') = l^2/alpha'" | V |
| SM 81 | `SMFTweight` | T_y, h_y(a) | NRH08 "the two OPE contributions assemble to h_y(a)/(z-w)^2" | V |
| SM 82 | `SMBRSTradial` | (∂²_y + (2/l)∂_y)f = 0 | NRH08 "the marginal roots of h_y are a = 0 and a = -2/l (modes {1, e^{-2y/l}})", "(d_y^2 + (2/l) d_y) f = 0 for f = W0 + W1 e^{-2y/l}" | V |
| SM 83 | `SMBRSTmomentum` | h_y(p_y), roots p_y = 0, 2i/l | NRH08 "h_y(i p_y) = (alpha'/4) p_y (p_y - 2 i/l) and P_y = p_y - i/l gives (alpha'/4)(P_y^2 + 1/l^2)" | V |
| SM 84 | `SMBRSTvertex` | Radially dressed W1 vertex and its linearized conformal weights | NRH08: radial conformal-weight and longitudinal-contraction checks. The undressed W0 contraction test alone does not establish exact marginality or self-OPE regularity of the dressed W1 vertex. | V° |
| SM 85 | `SMBRSTcentral` | c_y = 1 + 6α′/l², total 3(k+2)/k | NRH08 "c_{beta gamma} + c_y = 2 + (1 + 6/k) = 3(k+2)/k" | V |
| SM 86 | `SMWgaugeobstruction` | W₀ gauge, e^{−2y/l}W₁ not | NRH08 "Lhat_xi H^infty = D H^infty + H^infty D^T with D = ((-(dv)^T, 0), (b, dv)), b = d lambda~", "the conditions force d_y v^mu = 0, d_- v^+ = 0 = d_+ v^-, b_{+-} = -varpi/2, b_{+y} = d_+ v^y, b_{-y} = -d_…" | V |
| SM 87 | `SMBRSTfusion` | self-contraction, h_n, resonances | NRH08 "e^{a y(z)} e^{a y(0)} ~ \|z\|^{-alpha' a^2} = \|z\|^{-4 alpha'/l^2} for a = -2/l (from <y y> = -(alpha'/2) Log\|…", "the n-fold fused weight h_n = n + q n(1-n) equals n + h_y(-2n/l) with q = alpha'/l^2" (+1 more) | V |
| SM 88 | `SMupliftblocks` | Ten-dimensional uplift blocks and sector curvature cancellation | NRH09: arbitrary-chiral R uplift and the one-sided NR uplift L_-=0 with arbitrary W0/W1. No direct full two-sided NR uplift run. | V° |
| SM 89 | `SMDFTKilling` | DFT Killing equations | Defining equations; NRH09 checks their explicit vacuum realization, not arbitrary backgrounds. | D |
| SM 90 | `SMtypeIIKS` | Type-II Killing-spinor equations | Stated framework; NRH09 checks selected reductions and Clifford algebra, not the complete current ten-dimensional system in the real Majorana basis. | D |
| SM 91 | `SMsusyclosure` | Supersymmetry closure and the Killing parameter from a spinor bilinear | NRH09: reduced vacuum bilinear and vanishing generalized Lie derivative; not the full supersymmetry commutator with compensating transformations. | V° |
| SM 92 | `SMsemicov` | semi-covariant connection | NRH05 "SMbackgroundconnection on R: nabla_C P_AB = 0 through z^2", "SMbackgroundconnection on NR: nabla_C P_AB = 0 through z^2" (+4 more); NRH09 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2 (pairwise cancellation with AdS3)", "S_(0)(S3 block) = +4/l^2 via the closed form (cancels -4/l^2 of either 3d saddle)" | V° |
| SM 93 | `SMuplift` | AdS₃×S³×R⁴ with flux | NRH09 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2 (pairwise cancellation with AdS3)", "S_(0)(S3 block) = +4/l^2 via the closed form (cancels -4/l^2 of either 3d saddle)" (+1 more) | V |
| SM 94 | `SMRDFTKilling` | Ordinary-field form of the generalized Killing equations | Cited decomposition; the ordinary Killing vectors are not enumerated in the suite. | S |
| SM 95 | `SMRlocaliso` | k^{(ij)} = s_is_j stabilizers | NRH09 "k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' - (l^2/4) k''' = 0 (the stabilizer equation)" | V |
| SM 96 | — | vol₃, H₃ = −(2/l)vol₃, H_{y−+} = +(2/l)√\|g₃\| (flux orientation) | NRH09 "flux orientation: H_{y-+} = +(2/l) Sqrt[\|g_3\|] = (2/l)(e^{2y/l} - L+L- e^{-2y/l})" | V |
| SM 97 | — | AdS₃ Killing-spinor equation ∇_μ ε± = ±(1/2l)γ_μ ε± | — | S |
| SM 98 | `SMRcomponentHill` | first-order pair → Hill equation | NRH09 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s" | V |
| SM 99 | `SMRlocalKS` | Hill equations for s± | NRH09 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s", "massless BTZ (L = 0): the (u, v) pair shifts by exactly 2 pi v over one circuit (unipotent)" (+1 more) | V |
| SM 100 | — | N_local = 16 real = 8_{Spin(1,9)} + 8_{Spin(9,1)} | — | S |
| SM 101 | `SMvielbein` | exact non-Riemannian double vielbein | NRH05 "SMexactNRdata: the lowered SMvielbein frames reconstruct P, Pbar of the exact hairy metric and are orthogonal" | V |
| SM 102 | `SMvielbeincheck` | defining relations | NRH09 "V_M^p V_Np = P_MN and Vbar_M^pbar Vbar_Npbar = Pbar_MN for the exact frame , any hh, sigma, W", "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)" (+1 more) | V |
| SM 103 | `SMNRlocalstabilizer` | Non-Riemannian local stabilizer conditions | NRH09: epsilon=c/sqrt(L) and the exact one-sided transformation identity. The full arbitrary two-sided W0/W1 stabilizer system is not checked. | V° |
| SM 104 | `SMexactiso` | vacuum isometries | NRH09 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm", "X^M = (omega_+, 0, -(l/2) v'; v^+, 0, -(l/2) v') with v^+ = 2 f1 f2, omega_+ = -2 l^2 f1' f2'" | V |
| SM 105 | `SMweighteddilaton` | weighted dilaton condition | NRH09 "Lhat_xi d = 0 (the radial component compensates the divergence)", "equivalently d_M(e^{-2d} xi^M) = 0" | V |
| SM 106 | `SMkillingspinor` | vacuum Killing spinor E = (√2f, l∂₊f) | NRH09 "SM: vacuum spin connection = displayed Phi_{~+ oplus y} = 1/(2l), Phi_{+ ominus y} = -1/l (raised slots)", "D_{pbar} E = 0 for E = (Sqrt[2] f(x+), l f'(x+)), arbitrary chiral f" (+2 more) | V |
| SM 107 | `SMcomplexblocks` | σ_i, τ_i, ρ_m | NRH09 "rep: {gamma^p, gamma^q} = 2 eta^{pq} (3d lightcone blocks)" | V |
| SM 108 | — | index convention ε_α ≡ ε_{α₁α₂α₃α₄}, α₁,α₂,α₃ = 1,2, α₄ = 1,…,4 | — | D |
| SM 109 | `SMgammaten` | ten-dimensional Γ^p̂, Γ₁₁ | NRH09 "{Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32", "Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0" | V |
| SM 110 | `SMgammaMajorana` | Majorana intertwiner | NRH09 "BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, same for Gamma_11, and BB10 BB10^* = 1" | V |
| SM 111 | `SMgammabarred` | barred Clifford algebra | NRH09 "{Gammabar, Gammabar} = -2 eta, Gammabar_{pq} = -Gamma_{pq}, same Majorana intertwiner" | V |
| SM 112 | `SMreducedDirac` | reduced two-component system | NRH09 "SM text after SMreducedDirac: D_{pbar} E has the components 0, -L+ u e0/l, em0/Sqrt2, em1/Sqrt2 - W1 u e0/(…", "SM text after SMreducedDirac: gamma^p D_p E + E/(Sqrt2 l) has the components ey0/Sqrt2 + Sqrt2 e0/l and ep0…" | V |
| SM 113 | `SMinternalprojectors` | Internal projectors and complex product-spinor representatives | NRH09: historical complex rank count only; it does not verify the current product representatives or their real Majorana combinations. | V° |
| SM 114 | — | real Majorana combinations Ψ + B₁₀^{−1}Ψ*, i(Ψ − B₁₀^{−1}Ψ*) | NRH09 "the Majorana condition leaves 32 real components" | D |
| SM 115 | `SMspinorcountchain` | Majorana, Weyl and background-projector polarization counts | NRH09: algebraic ranks and arithmetic consistency. Arbitrary chiral functions are not discarded; this is not a construction of global spinor spaces or nonzero integrable supercharges. | V° |
| SM 116 | `SMhairyKS` | four arbitrary real chiral functions on each one-sided NR branch (hairy Killing spinors) | the reduced one-sided jet system in the printed frame is verified; the full ten-dimensional fermionic equations with the current internal basis and nonzero supercharges are not; NRH09 "the displayed system has rank six (seven relations, one dependent) and spans the derived equations", "the solution space is exactly {e_1, d_+ e_1}: e_0 = 0, e_1 = F_+(x^+) with d_+ e_1 unconstrained" (+1 more) | V° |
| SM 117 | — | Fourier expansion of the chiral functions F^r_±(x^±) with reality f_{−n} = f_n* | — | D |
| SM 118 | — | formal pairing δQ[F] = ∫dφ Σ_r F^r δS_r = 2π Σ f_n^r δS_{r,−n} | — | D |
| SM 119 | `SMcandidateD` | Gauge-covariant derivative for the boundary candidate | Definition; used in NRH10. | D |
| SM 120 | `SMcandidateaction` | candidate action and its reduction | NRH10 "the P^(0) term contains only D_+ and the Pbar^(0) term only D_-", "each term selects a single independent component of phi_A" | V |
| SM 121 | `SMcandidategauge` | Non-Abelian gauge transformation of the classical candidate | NRH10: non-Abelian bosonic sector; the entire covariant Grassmann system is not checked by this predicate. | V° |
| SM 122 | `SMcandidateWitt` | Witt transformations | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
| SM 123 | `SMcandidateWittL` | δ_vL = ∂(vL) | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
| SM 124 | `SMcandidatefermionic` | chiral fermionic transformations | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
| SM 125 | `SMcandidatefermionicL` | their total-derivative variations | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
| SM 126 | `SMcandidateextra` | Grassmann-even chiral transformation | NRH10 "delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]" | V |
| SM 127 | `SMcandidateextraL` | its variation | NRH10 "delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]" | V |
| SM 128 | `SMcandidatetrivial` | equation-of-motion redundancies | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
| SM 129 | `SMcandidatetrivialL` | δ_αL total derivative | NRH10: Grassmann realization with ordinary partial derivatives and vanishing gauge connection (A=0); not a full non-Abelian covariant calculation. | V° |
