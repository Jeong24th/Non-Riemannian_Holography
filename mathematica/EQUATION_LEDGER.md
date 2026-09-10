# Equation ledger

Current manuscript SHA-256: `FCBE00570616A89820C3F969E8E97BA74454B8555DDA5EBC25ED56C268C5BE6F` (2026-09-10).
Letter (1)-(23), SM (1)-(141). LaTeX labels are stable identifiers.
Check locators below quote descriptive text in the named Mathematica file;
obsolete printed equation numbers have been removed from the executable checks.

**V**: exact symbolic verification in the stated domain; **V°**: consequences,
restricted background or stated asymptotic order only; **D**: definition;
**S**: cited statement or assumption; **N**: no check. A passing suite does not
turn definitions, assumptions or uncomputed completions into verified results.

## Letter

| No. | Label | Content | Where verified | Status |
|---:|---|---|---|:-:|
| 1 | `Rfields` | NS–NS Bañados family (frame ϑ±, metric, B, dilaton) | NRH02 "H J H = J (O(3,3) constraint)" (3 checks) (built from (1)); the EDFE and horizon checks; NRH06 "the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+" (2 checks); NRH07 "R(AdS3, Banados) = -6/l^2 and H^2 = -24/l^2 for arbitrary chiral L_pm (pairwise cancellation with S3…" | V |
| 2 | `RDFTfields` | DFT variables H_MN, e^{−2d} of (1) | NRH02 "H J H = J", "(e^{-2d})^2 = -det g", "e^{-2d} = e^{2y/l}(1 − L+L− e^{−4y/l})"; EDFE "G_MN = 2 l^-2 J_MN" | V |
| 3 | `Rboundary` | asymptotic falloffs | NRH02 "H − H^infty = O(e^{−2y/l})", "d + y/l = O(e^{−4y/l})" | V |
| 4a | `boundaryH` | H^∞, d^∞ | NRH02 "boundaryH…" (two checks) | V |
| 4b | `Rboundaryframe` | aligned D = 2 boundary frame | NRH02 "Rboundaryframe: boundary vielbeins reproduce P^(0) and Pbar^(0)" | V |
| 5 | `NRvariables` | (Π, q, e^σ, χ) and e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) | NRH03 "e^{sigma} sinh(chi/2) = …", "NRradialchange: d chi/d q…" | V |
| 6 | `NRHcompact` | the everywhere non-Riemannian H_MN | NRH03 "H J H = J for the exact non-Riemannian matrix (generic W)", "type (1,1) at every radius"; EDFE checks; NRH07 "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)" | V |
| 7 | `NRdilaton` | e^{−2d} = e^{2y/l}(1 − q²), d = −y/l + ln cosh(χ/2√2) | NRH03 "e^{-2d} e^{-2y/l} = 1 − q^2…"; "d chi/dy…" | V |
| 8 | `NRWgeneral` | exact hair profile W | NRH03 "EDFE scalar: S_(0) = −4/l^2 for ARBITRARY W", "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F", "G_MN = 2 l^-2 J_MN on the ODE shell", "NRWgeneral] solves d^2W/dchi^2 = F", "W_0 and W_1 multiply the two homogeneous modes", "near the boundary…", "falloff bookkeeping…" | V |
| 9 | `MainGprofile` | radial profile G(χ) = 4√2∫₀^χ[sinh t/sinh(t/√2) − √2]dt, G(0) = G′(0) = 0 | NRH03 "NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)", "NRGprofile: G'(0) = 0" | V |
| 10 | `RGKPW` | GKPW relation Z_DFT = Z_CFT | — (defining relation) | D |
| 11 | `Rrenvariation` | δS_ren on shell | NRH04 "delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) + d_K(h e^{-2d} A^K + 2 delta d e^{…" (the exact identity), "on shell (S_{p qbar} = 0, S_(0) = 2 Lambda = -4/l^2) the bulk term 2e^{-2d}(-delta d S_(0)) cancels …" (on-shell reduction), "delta S_ct = -(16 pi G)^{-1}(4/l) delta e^{-2d} = (16 pi G)^{-1} e^{-2d} 2 delta d (4/l): the +4/l s…"; NRH02 "RKdef: coefficient matching -2K = (16 pi G)^{-1} A^y and 2T_(0) = (16 pi G)^{-1} 2(B^y + 4/l)" | V |
| 12 | `Mainmomenta` | A^K_{pq̄} and B^K: generalized-metric and scalar radial momenta | NRH05 "on R/NR: B^y = 4 d_y d"; A^K enters the SM (13) check | V |
| 13a | `RKdef` | ⟨K⟩, ⟨T₍₀₎⟩ from the rescaled momenta | NRH02 "RKdef: coefficient matching -2K = (16 pi G)^{-1} A^y and 2T_(0) = (16 pi G)^{-1} 2(B^y + 4/l)"; **nonlinear evaluation on both exact saddles**: NRH02 "on the exact family: −(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L−},{1, L−}}/(16 pi G l)", "…e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0"; NRH04 "on the exact NR family: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, W1/4},{0, L-}}/(16 pi G l)…" (2 checks), "…hence <T_(0)> = 0"; NRH04 "SMscalarcutoffresponse: 2(16 pi G)^{-1} e^{-2d}(B^y + 4/l) with B^y = 4 d_y d = -4/l + 4 d_y delta d…" (2 checks) | V |
| 13b | `RDFTconservation` | T^CFT_AB and ∇^A T_AB = 0 | NRH02 "div T = {…} exactly"; "SMinvariantstress…" | V |
| 14 | `Rcontinuity` | Ward identities | NRH02 "div T…", "no local condition on K_{op bom}" | V |
| 15 | `Rkilling` | common R/NR asymptotic-symmetry generator (ξ^±_R with the l²e^{−2y/l}∂²ε/4 tail, ξ^±_NR = ε^±; common radial and dual components) | NRH02 "Lhat_xi d = O(e^{−4y/l})", "Lhat_xi H − delta_eps H = O(e^{−4y/l})"; NRH05 full generator checks; the NR form ξ^± = ε^± underlies the NRH03 checks of (17) | V |
| 16 | `RVirasoro` | δ_εL± with −(l²/4)∂³ε | NRH02 "Lhat_xi H - delta_eps H = O(e^{-4y/l})", "c = 3l/2G reproduces delta_eps T = eps T' + 2 T eps' - (c/12) eps'''"; NRH05 "Virasoro cocycle: (charge bracket density) + (l^2/4) e1 e2''' is a total derivative", "(iii) R bracket − adjoint − central…" | V |
| 17 | `NRasympt` | δ_εL±, δ_εW₁ at W₀ = 0 | NRH03 "NRasympt…" (four checks) | V |
| 18 | `Maincharge` | Q[ε] = (1/4πGl)∮(L₊ε⁺dx⁺ − L₋ε⁻dx⁻), common to both branches | NRH05 "lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''" (2 checks) (components) and "(iii) normalization chain (16 pi G)^{-1}(4/l) = 1/(4 pi G l)" | V |
| 19 | `Mainonepoints` | 16πGl⟨K⟩_R = ((L₊, L₊L₋),(1, L₋)), 16πGl⟨K⟩_NR = ((L₊, W₁/4),(0, L₋)), ⟨T₍₀₎⟩ = 0 | NRH02 "on the exact family: −(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L−},{1, L−}}/(16 pi G l)", "…e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0"; NRH04 "on the exact NR family: … = {{L+, W1/4},{0, L−}}/(16 pi G l)", "…hence <T_(0)> = 0" | V |
| 20 | `Mainframevariation` | δV_Mp = ½V̄_M^q̄h_pq̄, δV̄_Mq̄ = −½V_M^p h_pq̄ (rotation-free gauge; = SM 6) | NRH04 "delta(V eta V^T) = delta H/2…" | V |
| 21 | `Mainfivebyfive` | connected 5x5 matrix: particular A/M entries + R_s + local contacts | Definition; NRH04 checks the A/M entries. R_s includes hair self-response and other Ward-undetermined terms. | D |
| 22 | `Rcorrelators` | particular A_NR, A_R, M_NR, M_R kernels, with coefficients at x | NRH04 "Rcorrelators:" checks: Ward operators applied to the chosen inverse derivative; source-sign negative control. | V |
| 23 | `Mainworldsheet` | ℒ_GO = (1/2πα′)(β∂̄x⁺ + β̄∂x⁻ + ∂y∂̄y), V_𝒲 = (𝒲/4πα′)∂x⁺∂̄x⁻ | NRH06 "at chi -> 0 the W coupling reduces to…", "with the 1/(2 pi alpha') prefactor this is V_W = …" | V |

## Supplemental Material

| No. | Label | Content | Where verified | Status |
|---:|---|---|---|:-:|
| 1 | `SMpert` | H → H + δH, d → d + δd | — | D |
| 2 | `SMprojectors` | P, P̄ | NRH01 (used throughout); NRH02 "Gamma compatibility: nabla_C P_AB = 0" | D |
| 3 | `SMmixedfluctuation` | h_{pq̄} := δH_MN V^M_p V̄^N_q̄ | — | D |
| 4 | `SMcosetreconstruction` | δH = 2V_(M^p V̄_N)^q̄ h_pq̄ | NRH04 "delta(V eta V^T) = delta H/2…"; the constrained families of NRH04 (H(r,w), H4) | V |
| 5 | — | δV, δV̄ decomposition | NRH04 "delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2 for generic h_{p qbar}" | V |
| 6 | `SMframevariation` | gauge-fixed δV_{Mp} = ½V̄_M^q̄ h_pq̄, δV̄ = −½V h | NRH04 "delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2 for generic h_{p qbar}" | V |
| 7 | — | δL ≃ 2e^{−2d}(−h K + δd T₍₀₎) | NRH02 "RKdef: coefficient matching…" (its consequence) | V° |
| 8 | `SMresponsedef` | K_{pq̄}, T₍₀₎ definitions | — | D |
| 9 | `SMKvariation` | K = −e^{2d}(δL/δH)VV̄ = −½e^{2d}δL/δh | — (rewriting of (7)–(8)) | D |
| 10 | `SMconnectionvariation` | δΓ_KMN in terms of δP and δd (Eq. (2.56) of Park:2025core) | — (quoted from Ref. Park:2025core; not re-derived in this archive) | S |
| 11 | `SMsixprojectors` | six-index projectors 𝒫, 𝒫̄ (Eq. (2.52) of Park:2025core) | — (definition) | D |
| 12 | `GammaDFT` | L_Γ² = e^{−2d}S₍₀₎ − ∂(e^{−2d}B) = e^{−2d}(PP − P̄P̄)(ΓΓ…) | NRH05 "on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)", "on NR (arbitrary W)…" (the toolbox Γ² density is the displayed quadratic form) | V |
| 13 | `variation` | exact variation of L_Γ² | NRH04 "setup: the BTZ Riemannian double vielbein obeys…", "setup: Gamma^L_{LN} = −2 d_N d…", "delta L_Gamma2 = 2e^{-2d}(h S − delta d S_(0)) + d_K(h e^{-2d} A^K + 2 delta d e^{-2d} B^K) exactly, generic tangential h and delta d on BTZ", "the identity is not vacuous…" | V |
| 14 | `defB` | B^K and A^K_{pq̄} | as for SM 13; NRH05 "…B^y = 4 d_y d" | V |
| 15 | `SMrenvariation` | δ(S_Γ² + S_ct) at the cutoff | NRH04 "on shell…", "delta S_ct…" | V |
| 16 | `SMsecondvariation` | δ₂δ₁S_ren at the cutoff: second variation of the boundary expression (15) | NRH04 "SMsecondvariation:" checks local coset coordinates, independent h1, and the nonzero constrained mixed variation of H. The displayed boundary bilinear is a definition. | D/V |
| 17 | `SMmomentumvariation` | δA^K_{pq̄}, δB^K including the frame-variation terms; 𝒜^K_MN | NRH04 "SMmomentumvariation:" differentiates both frame projections of A and all three terms of B on the vacuum; a frozen-frame negative control is nonzero. | V° |
| 18 | `SMsources` | canonical sources 𝒥 = (J₊, J₋, J_h, J_r, J_d) = 2h^{(0)} components, J_r = W₀ | NRH04 "J^{op bop} = … hence the Hessian factor…" (the canonical-source factor) | V° |
| 19 | `SMradialhessian` | five-component response matrix from the radial derivatives (definition, FP scheme) | — (definition) | D |
| 20 | `SMinfinityvielbein` | limiting D = 3 frame | NRH04 "V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty", "V and Vbar are mutually orthogonal"; NRH07 "at hh = 0 and W_0 = 0, lowering the local indices with eta, etabar gives the limiting frame" | V |
| 21 | `SMflatmetrics` | η, η̄ | NRH04 "V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty"; NRH07 "rep: {gamma^p, gamma^q} = 2 eta^{pq} (3d lightcone blocks)" | V |
| 22 | `SMFG` | Fefferman–Graham gauge | imposed in every NRH04 linearization | D |
| 23 | `SMFGcount` | 9 − (3 + 2) = 4 | NRH04 "coset count…" | V |
| 24 | `SMfixedprojection` | fixed-frame projection | NRH04 "exact h^R_{p qbar} matches the displayed closed form", "falloffs h^{(2)} = {2L+, 2L+L-; 2, 2L-} + O(e^{-6y/l})", "exact h^NR = {{e^s sinh chi, (W/2) cosh chi, 0},{0, e^{-s} sinh chi, 0},{0,...}}", "falloffs {2L+ /u, (W0 + W1/u)/2; 0, 2L- /u} with O(u^-3) diagonals" | V |
| 25 | `SMNRlin` | linearized EDFE (nine equations) | NRH04 "solves the displayed system", "reproduces the linearized EDFE…" (re-derived from the nonlinear G_MN) | V |
| 26 | `SMNRradialintegration` | elementary radial identities | NRH04 "D_y{1, e^{−2y/l}} = 0, …" | V |
| 27 | `SMNRsol` | general solution | NRH04 "solves the displayed system", "reproduces the linearized EDFE: G^{(1)}_MN = 0 on the general solution", "SMcontinuitycheck: d_- h^(2)_{op bop} = 0, d_+ h^(2)_{om bom} = 0, d_pm h^(2)_{om bop} = 0", "d_y delta d = -(l/8) d_+ d_- h^(0)_{om bop} in the general solution, so the coefficient is -(l/(16 pi G)) e^{2y/l} d_+ d…" | V |
| 28 | `SMcrosspairs` | source–response pairs | NRH04 "UV momenta…", "the finite piece equals…" (all four pairings) | V |
| 29 | `SMlogfreeconditions` | fixed-dilaton, log-free conditions | NRH04 "the log coefficients are exact multiples of the three r^(0) conditions"; the remaining diagonal combination is the h_pm linear coefficient in the verified general solution SM (27) | V° |
| 30 | `SMrwquadratic` | crossed quadratic action | NRH04 "Gamma^2 quadratic density = −(1/2) e^{2y/l} dy r dy w" | V |
| 31 | `SMrwvariation` | its UV variation | NRH04 "UV momenta give…" | V |
| 32 | `SMSren2` | finite cutoff variation, all four channels | NRH04 "the four-channel coset family obeys H J H = J through O(t^2)" (7 checks) (six checks: coset family, restricted on-shell sector, vacuum momentum, divergent pieces, finite piece, non-triviality) | V |
| 33 | `SMscalarcutoffresponse` | scalar cutoff coefficient | NRH04 "SMscalarcutoffresponse: …", "d_y delta d = …" | V |
| 34 | `SMresponsematrix` | ⟨K_ab̄⟩ = h⁽²⁾_ab̄/(32πGl) | NRH04 "<K_{a bbar}> = h^{(2)}/(32 pi G l) => stress one-points"; normalization from "…finite piece" | V |
| 35 | `SMlinearizedmomenta` | δA^y_{ab̄} = ½∂_y h_{ab̄} + O(e^{−4y/l}), δB^y = 4∂_y δd + O(e^{−4y/l}) | NRH04 "SMlinearizedmomenta:" differentiates the vacuum and both R/NR backgrounds through exp(-2y/l), retaining arbitrary tangential h, delta d and state profiles. | V° |
| 36 | `SMgeneralresponsevariation` | δ⟨K⟩, δ⟨T₍₀₎⟩ with the log-branch coefficients h^{(2L)}, δd^{(2)}, δd^{(2L)} | NRH04 "SMgeneralresponsevariation:" extracts the finite radial-momentum coefficients, including logarithmic branches; finite local counterterms are excluded. | V° |
| 37 | `SMnpointresponse` | connected hierarchy with J^I = −2h^{(0)I} | NRH04 "J^{op bop} = … hence the Hessian factor…" (the canonical-source factor) | D/V° |
| 38 | `SMstresshessian` | Hessian with factor 1/2 | NRH04 "J^{op bop} = -2 h^(0)op bop = 2 h^(0)_{om bom} and J^{om bom} = 2 h^(0)_{op bop}; hence the Hessian …" | V |
| 39 | `SMcontinuitycheck` | ∂₋⟨K⊕⊕̄⟩ = 0, … | NRH04 "SMcontinuitycheck…" | V |
| 40 | `SMtwopointWard` | two-point Ward identities | same check (chirality of the normalizable coefficients) | V° |
| 41 | `SMexactprojectionR` | exact h^R_{pq̄} | NRH04 "exact h^R_{p qbar} matches the displayed closed form" | V |
| 42 | `SMexactprojectionRfalloff` | its falloffs | NRH04 "falloffs …" | V |
| 43 | — | δd_R = −½ln(1 − L₊L₋e^{−4y/l}) | NRH04 "delta d_R = (1/2) L+ L− e^{−4y/l} + O(e^{−8y/l})" | V |
| 44 | `SMonept` | Riemannian diagonal one-point functions | NRH04 "… stress one-points" | V |
| 45 | `SMinvariantstress` | T_{x⁺x̃₋} = 2K⊕⊕̄ | NRH02 "SMinvariantstress…" (two checks) | V |
| 46 | `SMPBHdata` | s₊, r₊ | NRH04 "s_+ = …", "r_+ = …", "the PBH variation carries no (om bop) source deformation" | V |
| 47 | `SMPBHkernel` | Dirichlet-to-Neumann kernel | NRH04 "(l^2/4) d_+^3 [−1/(2 pi x)] = 3 l^2/(4 pi x^4)", "away from coincidence…" | V |
| 48 | `SMRtwopt` | two-point functions, c = 3l/(2G) | NRH04 "(64 pi G l)^{-1} 3 l^2/(4 pi) = (8 pi)^{-2}(c/2)"; NRH02 "c = 3l/2G reproduces…" | V |
| 49 | `SMexactprojectionNR` | exact h^NR_{pq̄} | NRH04 "exact h^NR = …" | V |
| 50 | `SMexactprojectionNRfalloff` | its falloffs | NRH04 "falloffs …" | V |
| 51 | — | δd_NR = ln cosh(χ/2√2) = … | NRH04 "delta d_NR = (1/4) L+ L− e^{−4y/l} + O(e^{−8y/l})" | V |
| 52 | `SMNRonept` | NR one-point functions | NRH04 "NR hair channel…", "h^(2)_{om bop} = 0 … <K_{om bop}> = 0; …", and nonlinearly "on the exact NR family: … {{L+, W1/4},{0, L−}}/(16 pi G l)" | V |
| 53 | `SMNRhairhessian` | type-changing Hessian (kernel undetermined) | NRH04 "…: (1/(64 pi G l)) delta(W_1/2) = (1/(128 pi G l)) delta W_1" (the displayed factor; the kernel is stated to be undetermined) | V° |
| 54 | `SMNRstresskernel` | admissible diagonal kernels vanish in the restricted sector | consequence of SM 27 + 29 (NRH04): the chosen admissible normalizable diagonal kernels vanish in this restricted source sector | V° |
| 55 | `SMgeneralWardoperators` | 𝒟^{NR}_±, 𝒟^{R}_±, ℰ_± (the transformation laws (16)–(17) as operators) | NRH02 "Lhat_xi H - delta_eps H = O(e^{-4y/l})" and NRH03 "NRasympt: Lhat_xi H - delta_(L,W1) H = O(u^-2) componentwise" (4 checks) verify the underlying laws | V° |
| 56 | `SMgeneralparticularkernels` | particular kernels A^s_±, B^R_±, C^{NR}_± from the Ward operators | NRH04 "Rcorrelators:" applies the Ward operators to the chosen inverse derivative. The code uses M_R and M_NR directly, without B_R/C_NR aliases. | V |
| 57 | `SMgeneralpositionkernels` | position-space singular parts of A^{NR}_±, A^{R}_±, C^{NR}_± | NRH04 "Rcorrelators:" verifies all three pole structures and M_R = L_opposite A_R. | V |
| 58 | `SMCPSform` | surface-charge one-form | NRH05 `chargeOneForm` (all three terms), "lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary" (3 checks), "(i) integrability…" | V |
| 59 | — | Q_R | NRH05 "lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''" (2 checks) with "(iii) normalization chain: (16 pi G)^{-1} (4/l) = 1/(4 pi G l), the Letter's charge normalization" | V |
| 60 | `SMRpotentialcomponents` | lim e^{−2d}K̂ components (the minus-sector check uses K̂^{y+} = −K̂^{+y}) | NRH05 "lim e^{-2d} Khat^{-y}[eps+] = …", "…Khat^{y+}[eps-] = … (mirror)" | V |
| 61 | `SMCPSresult` | Θ̂ → 0, k = (4/l)εδL | NRH05 "lim e^{-2d} Thetahat^{+,-,y} = 0", "k^{-y}[eps+] = (4/l) eps+ dL+", "…k^{+y}[eps-…" | V |
| 62 | `SMNRchargefalloffs` | state-dependent falloffs | NRH05 "state-dependent falloffs…" | V |
| 63 | `SMNRchargecancellation` | componentwise W₁ cancellation | NRH05 "W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise" | V |
| 64 | `SMNRcharge` | δQ, Q | NRH05 "(i) k^{-y}[eps+] = delta[(4/l) eps+ L+]", "(i) mirror…" | V |
| 65 | `SMCPSalgebra` | {Q,Q} = Q[[ε,η]], c_charge = 0 (+ footnote) | NRH05 "the same-chirality C-bracket closes up to a closed B-gauge parameter (slot x~+ only)" (4 checks) (four checks), "(ii) {Q[e1+], Q[e2+]} - Q[[e1,e2]] is a total derivative: centerless plus sector" (3 checks), "(iv) the central cocycle is antisymmetric modulo total derivatives" (3 checks), "(v) the adjusted mixed bracket has a leftover (the raw closure fails, as it should)" (2 checks), "(vi) minus-sector C-bracket closes up to a closed B-gauge parameter (slot x~- only)" (2 checks), "footnote: with delta T = eps T' + 2 T eps' - l^2 eps''' the combination L + T/4 obeys the law with -…" | V |
| 66 | `SMgamma2` | Γ² identity | NRH05 "on R…", "on NR (arbitrary W)…" | V |
| 67 | `SMgamma2flux` | B^y = 4∂_y d, e^{−2d}B^y = −(4/l)(e^{2y/l} + μe^{−2y/l}) | NRH05 "on R…", "on NR…" (three checks) | V |
| 68 | `SMmudefinition` | μ dichotomy | NRH05 "on NR: −2 d_y e^{-2d} = … [mu-dichotomy]" | V |
| 69 | `SMgamma2cutoff` | S_ren(Y) | NRH05 "the regulated combination equals…" | V |
| 70 | `SMgamma2value` | S_ren = −8√μ/(16πGl)∫d²x | NRH05 "Y -> Infinity limit gives…" | V |
| 71 | — | Killing horizon e^{4y/l} = L₊L₋ | NRH02 "horizon: e^{-2d} = 0 at e^{4y/l} = L+ L−"; NRH05 "endpoints…" | V |
| 72 | `NRhill` | ψ±, A± = ¼(∂lnL)² − ½∂²lnL = ψ″/ψ | NRH03 "NRhill: …"; A± enters the verified source F | V |
| 73 | `NRradialchange` | q, χ, ∂_y, ∂_q | NRH03 "NRradialchange: d chi/d q = …", "d chi/dy = …" | V |
| 74 | `NRradialoperator` | radial operator identity | NRH03 "NRradialoperator…" | V |
| 75 | `NRchiODE` | ∂²_χ W = F | NRH03 "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F" | V |
| 76 | `NRsource` | the source F | same check (F as displayed) | V |
| 77 | `NRg` | ρ(χ) | NRH03 "NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)" | V |
| 78 | `NRGprofile` | G(χ) normalization | NRH03 "NRGprofile: G'(0) = 0", "…(e^{s chi} − 1 − s chi)…" | V |
| 79 | — | q = √(L₊L₋/2) μ_RG^{−2}, χ, μ_RG = e^{y/l} | — (definitions; used in the next line) | D |
| 80 | — | μ_RG dχ/dμ_RG = −2√2 sinh(χ/√2) | NRH03 "SM: RG rapidity mu d chi/d mu = −2 Sqrt[2] Sinh[chi/Sqrt[2]]" | V |
| 81 | `SMBtransform` | finite B-transformation | NRH03 "SMBtransform, SMWshift…", "Omega_b is O(3,3)" | V |
| 82 | `SMWshift` | Ω_bH(W)Ω_bᵀ = H(W − 2b_{+−}) | NRH03 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| 83 | `SMdyg` | doubled-yet-gauged action | — (defining action; its reductions are verified below) | D |
| 84 | `SMphysicalsectionA` | Â_αμ, D_αx^M on the section | — | D |
| 85 | `SMRfirstorder` | Riemannian first-order form | NRH06 "the auxiliary equations give…", "eliminating the auxiliaries reproduces E_{mu nu}…" | V |
| 86 | `SMceff` | c_eff² = 2F | NRH06 "−det g_par = …", "c_eff^2 = 2F -> 4 Sqrt[L+L−…"; NRH02 "c_eff^2 := 2F…" | V |
| 87 | `SMlongstringE` | winding-string energy of the static probe (SM 5.1) | NRH06 "before : … d_t is Killing…", "before : gamma^{tau a} g_{t nu} d_a X^nu = … = 1", "before : B_{t phi} = …", "before : E = −Int d sigma P_t … reproduces E(y) = …", "−det g_(t,phi) = l^2 (e^{2y/l} − L+L− e^{−2y/l})^2 per winding…", "E(y) = …", "remark: with phi_0 = 0 the area density equals l e^{-2d}…" | V |
| 88 | `SNCtau` | SNC clock forms | NRH03 "SNCtau: tau+ = dx+ − L_− e^{−2y/l} dx− + O(e^{−4y/l})", "unit clock determinant", "W-part of the lower-right block = W(tau+ tau− + tau− tau+)"; NRH07 "x rows = Sqrt[2] tau^pm…" | V |
| 89 | `SMdygconstraints` | τ⁺·∂̄x = 0 = τ⁻·∂x | NRH06 "SM: H^{mu nu} tau^pm_nu = 0", "SM: the dual vectors Y, Ybar exist at every radius" (the kernel and multiplier structure) | V° |
| 90 | `SMdygGO` | exact reduced Lagrangian | NRH06 "symmetric-block route − antisymmetric-clock route = W (tau− . dx)(tau+ . dbar x)", "at chi -> 0 the W coupling reduces to…" | V° |
| 91 | `SMGO` | Gomis–Ooguri limit | NRH06 "at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)"; contraction bookkeeping | V |
| 92 | `SMvertex` | V_W = (1/4πα′)W∂x⁺∂̄x⁻ | NRH06 "with the 1/(2 pi alpha') prefactor this is V_W = …" | V |
| 93 | `SMWisB` | H(W) = Ω_bH(0)Ω_bᵀ, b_{+−} = −W/2 | NRH03 "SMWisB…" | V |
| 94 | `SMdeltaL` | antisymmetric clock coupling | NRH06 "symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)" | V |
| 95 | `SMFTsector` | S_y, S_FT | — (definitions; used in SM 97) | D |
| 96 | `SMWZWlevel` | k = l²/α′ | NRH06 "k = \|Int_{S^3} H\| / (4 pi^2 alpha') = l^2/alpha'" | V |
| 97 | `SMFTweight` | T_y, h_y(a) | NRH06 "the two OPE contributions assemble to h_y(a)/(z−w)^2" | V |
| 98 | `SMBRSTradial` | (∂²_y + (2/l)∂_y)f = 0 | NRH06 "the marginal roots…", "(d_y^2 + (2/l) d_y) f = 0…" | V |
| 99 | `SMBRSTmomentum` | h_y(p_y), roots p_y = 0, 2i/l | NRH06 "h_y(i p_y) = … P_y = p_y − i/l…" | V |
| 100 | `SMBRSTvertex` | U_{W₁} is a weight-(1,1) primary | SM 97–99 plus NRH06 "<x^+ x^−> = 0…" (W₁ has weight (0,0)) | V° |
| 101 | `SMBRSTcentral` | c_y = 1 + 6α′/l², total 3(k+2)/k | NRH06 "c_y/2 = …", "c_{beta gamma} + c_y = …" | V |
| 102 | `SMWgaugeobstruction` | W₀ gauge, e^{−2y/l}W₁ not | NRH06 "Lhat_xi H^infty = D H^infty + H^infty D^T…", "the conditions force…", "(db)_{+−y} = …", "v^y = … so closure forces d_y varpi = 0" | V |
| 103 | `SMBRSTfusion` | self-contraction, h_n, resonances | NRH06 "e^{a y(z)} e^{a y(0)} ~ \|z\|^{-alpha' a^2} = \|z\|^{-4 alpha'/l^2} for a = -2/l (from <y y> = -(alpha'/…" (3 checks) (three checks) | V |
| 104 | `SMupliftblocks` | H₁₀ = H₃ ⊕ H_{S³} ⊕ H_{R⁴}, d₁₀ = d₃ + d_{S³} | NRH07 `d10Assemble` and the symbolic ten-dimensional identities "S_(0)^{(10)} = 0", "(P S Pbar)^{(10)} = 0" for both uplifts (arbitrary chiral L±; one-sided hairy family with arbitrary W₀, W₁) | V |
| 105 | `SMDFTKilling` | generalized Killing equations | — (definition; `GenLieH`/`GenLieD` implement it: NRH07 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm", "Lhat_xi d = 0 (the radial component compensates the divergence)" (2 checks)) | D |
| 106 | `SMtypeIIKS` | type-II Killing-spinor systems | — (definition of the systems solved in NRH07) | D |
| 107 | `SMsusyclosure` | closure on ĥL_X, X = iε̄₂Γε₁ | NRH07 "C = i sigma_2 is the Majorana conjugation…", "X^M = …", "Lhat_X H^infty = 0 and Lhat_X d = 0…", "the bilinear is symmetric…" (three-dimensional calibration) | V° |
| 108 | `SMsemicov` | semi-covariant connection | NRH01 `GammaDFT`; NRH02 "Gamma compatibility: nabla_C P_AB = 0", "Gamma dilaton trace…", "Gamma torsionless…" | V |
| 109 | `SMuplift` | AdS₃×S³×R⁴ with flux | NRH07 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2", "R_{mu nu} = (1/4) H H on the S3 factor", "R(AdS3, Banados) = −6/l^2 and H^2 = −24/l^2…", "R_{mu nu} = (1/4) H H on the AdS3 factor", "S_(0)(S3 block) = +4/l^2", "S_(0)(R4 block) = 0" | V |
| 110 | `SMRDFTKilling` | ordinary-field form of (94) | — (decomposition statement; the six Killing vectors are not displayed) | S |
| 111 | `SMRlocaliso` | k^{(ij)} = s_is_j stabilizers | NRH07 "k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' − (l^2/4) k''' = 0" | V |
| 112 | — | vol₃, H₃ = −(2/l)vol₃, H_{y−+} = +(2/l)√\|g₃\| | NRH07 "flux orientation: H_{y−+} = +(2/l) Sqrt[\|g_3\|…" | V |
| 113 | — | ∇_με± = ±(1/2l)γ_με± | — (the torsionful Killing-spinor equation of the Riemannian branch; its consequence SM 114–115 is checked) | S |
| 114 | `SMRcomponentHill` | first-order pair → Hill equation | NRH07 "the first-order pair … closes into (l^2/2) s'' = L s" | V |
| 115 | `SMRlocalKS` | Hill equations for s± | NRH07 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s", "global AdS3 (L = −1/4)…", "massless BTZ (L = 0): the (u, v) pair shifts by exactly 2 pi v over one circuit (unipotent)", "constant L > 0: monodromy multipliers e^{pm 2 Sqrt[L0] pi} are real and not +-1: (0,0) kernel" | V |
| 116 | — | N_local = 16 = 8 + 8 | NRH07 "counting: (4,4) + (4,4) = the sixteen constant vacuum modes…" | V° |
| 117 | `SMvielbein` | exact non-Riemannian double vielbein | NRH07 "V_M^p V_Np = P_MN and Vbar_M^pbar Vbar_Npbar = Pbar_MN for the exact frame , any hh, sigma, W", "V^M_p Vbar_{M qbar} = 0 and P + Pbar = J", "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)", "at hh = 0 and W_0 = 0, lowering the local indices with eta, etabar gives the limiting frame", "x rows = Sqrt[2] tau^pm (), x~ rows = the dual vectors -Y/Sqrt[2], Ybar/Sqrt[2], and the W entries =…" | V |
| 118 | `SMvielbeincheck` | defining relations | same checks | V |
| 119 | `SMNRlocalstabilizer` | local stabilizer system | NRH07 "eps = c/Sqrt[L] solves…", "exact one-sided identity Lhat_xi H = delta H, with weights (1,2)"; NRH03 "line 2: … delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i" | V |
| 120 | `SMexactiso` | vacuum isometries | NRH07 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm"; "X^M = (omega_+, 0, -(l/2) v'; v^+, 0, -(l/2) v') with v^+ = 2 f1 f2, omega_+ = -2 l^2 f1' f2'" | V |
| 121 | `SMweighteddilaton` | weighted dilaton condition | NRH07 "Lhat_xi d = 0…", "equivalently d_M(e^{-2d} xi^M) = 0" | V |
| 122 | `SMkillingspinor` | vacuum Killing spinor E = (√2f, l∂₊f) | NRH07 "SM: vacuum spin connection = displayed…", "D_{pbar} E = 0…", "gamma^p D_p E = E/(Sqrt[2] l)" | V |
| 123 | `SMcomplexblocks` | σ_i, τ_i, ρ_m | NRH07 (matrices used in SM 124–126 checks); "rep: {gamma^p, gamma^q} = 2 eta^{pq}" | V |
| 124 | `SMgammaten` | ten-dimensional Γ^p̂, Γ₁₁ | NRH07 "{Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32", "Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0" | V |
| 125 | `SMgammaMajorana` | Majorana intertwiner | NRH07 "BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, …" | V |
| 126 | `SMgammabarred` | barred Clifford algebra | NRH07 "{Gammabar, Gammabar} = −2 eta, Gammabar_{pq} = −Gamma_{pq}, same Majorana intertwiner" | V |
| 127 | `SMreducedDirac` | reduced two-component system | NRH07 "reduced system…", "the opposite channel E = (0, g(x+)) has eigenvalue −1/(Sqrt[2] l)", "both channels also satisfy D_{pbar} E = 0" | V |
| 128 | `SMinternalprojectors` | ζ±, κ±, Ξ_{+r}, Ξ′_{−r}; S³ spinor equations | NRH07 "Weyl + S^3-line + zeta_+ leave complex dimension 4 (the Xi_{+r} span)", "below : … integrable (flat connection)" | V° |
| 129 | `SMspinorcountchain` | 32_C → 32_R → 16_R → 4_R | NRH07 "the Majorana condition leaves 32 real components", "…Weyl condition leaves 16…", "…complex dimension 4…" | V |
| 130 | `SMhairyKS` | one-sided Killing-spinor families | NRH07 "the displayed system has rank six…", "the solution space is exactly {e_1, d_+ e_1}…", "the mechanism…" (the displayed jet system; the ten-dimensional evaluation quoted in the SM is the Python/SymPy guard and the co-author notebooks) | V° |
| 131 | `SMcandidateD` | gauge derivative 𝔻 | — (definition; used in NRH08) | D |
| 132 | `SMcandidateaction` | candidate action and its reduction | NRH08 "the P^(0) term contains only D_+…", "each term selects a single independent component of phi_A", "gamma^oplus … rank one" | V |
| 133 | `SMcandidategauge` | internal gauge invariance | NRH08 "delta_Lambda L = 0 for the non-abelian bosonic sector" | V |
| 134 | `SMcandidateWitt` | Witt transformations | NRH08 "delta_v L = d_+(v^+ L) + d_−(v^− L)"; "Lhat_xi H^(0) = 0 <=> d_− v^+ = 0 = d_+ v^−" | V |
| 135 | `SMcandidateWittL` | δ_vL = ∂(vL) | same | V |
| 136 | `SMcandidatefermionic` | chiral fermionic transformations | NRH08 "delta_{eps+} L = …", "delta_{eps−} L = …", "the fermionic parameters carry arbitrary chiral profiles" | V |
| 137 | `SMcandidatefermionicL` | their total-derivative variations | same | V |
| 138 | `SMcandidateextra` | Grassmann-even chiral transformation | NRH08 "delta_zeta L = …" | V |
| 139 | `SMcandidateextraL` | its variation | same | V |
| 140 | `SMcandidatetrivial` | equation-of-motion redundancies | NRH08 "the transformations are built from the fermion equations of motion…" | V |
| 141 | `SMcandidatetrivialL` | δ_αL total derivative | NRH08 "delta_alpha L = d_−(psi^+ delta_alpha psi^+) + d_+(psi^− delta_alpha psi^−) identically (off shell)", "the two-component fermions used here are genuinely Grassmann…" | V |

## Prose statements verified in addition to the numbered displays

* Letter, below Eq. (5): the expansion e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) (NRH03).
* Letter, below Eq. (8): the boundary behaviour of W (exact in the one-sided limits; O(e^{−6y/l}) for constant L±) via the χ-expansion of the two homogeneous modes and the orders of G(χ) and e^{sχ} − 1 − sχ (NRH03).
* Letter, "Non-Riemannian Saddle": the B-shift with 2∂_[+λ_−] = W₀/2 removes W₀ (NRH03, SM 81–82 checks).
* Letter, "Riemannian Saddle": the BTZ horizon e^{2y/l} = √(L₊L₋) and c_eff (NRH02).
* SM 1: the vacuum has no first-order boundary momentum (no vacuum one-point function) (NRH04, SM (32) block).
* SM 3: the endpoints e^{−2d} = 0 for both branches (NRH05).
* SM 5.1: E_{+−} = 0, E_{−+} = −2F, E_{±±} = 2L± (NRH02, NRH06); −det g_∥ = F² − 4L₊L₋ = e^{−4d_R} (NRH06).
* SM 5.1, before (87): the unnumbered Noether-charge display P_μ = (2πα′)⁻¹(−√−γ γ^{τa}g_{μν}∂_aX^ν + B_{μν}X′^ν), E = −∫dσ P_t = (2πα′)⁻¹∫dσ(√−γ − wB_{tφ}), the identity γ^{τa}γ_{aτ} = 1, the Killing property of ∂_t, B_{tφ} = l(e^{2y/l} + L₊L₋e^{−2y/l}), −det g_(t,φ) = l²(e^{2y/l} − L₊L₋e^{−2y/l})², and the reproduction of (87) from the displayed momentum (NRH06, four checks).
* SM 5.2: the Gomis–Ooguri contraction table behind the exact marginality of V_W[W₀] (NRH06).
* SM 6: the symbolic ten-dimensional identities S₍₀₎⁽¹⁰⁾ = 0 = (PSP̄)⁽¹⁰⁾ for both uplifts (NRH07); the displayed vacuum spin-connection components; the Hill-monodromy counts; the complementary-halves bookkeeping.
* SM 6, text below (117): e^{±σ}sinh(χ/2) = uL±[1 + u²Π/3 + …] (NRH03).

## What the suite does not decide

* The SM (13) identity is checked on the exact BTZ background (constant L±) in its own Riemannian double-vielbein frame, for the general tangential fluctuation h_{ab̄}(x⁺,x⁻,y) of the Fefferman–Graham gauge and an arbitrary δd(x⁺,x⁻,y); the radial fluctuation components h_{aȳ}, h_{yb̄}, h_{yȳ} are not switched on.  It is a general identity of the semi-covariant formalism (Ref. Park:2025core), and its consequences (15), (32)–(34), and the nonlinear Eq. (13a) evaluations are verified independently.
* SM (110) and SM (113) are decompositions of DFT statements into ordinary Riemannian ones; the ordinary Killing vectors are not displayed and were not enumerated here.
* SM (130): the full ten-dimensional Killing-spinor evaluation with hair is quoted in the SM from the Python/SymPy guards (`checks/verify_10d_killing_spinor.py`, `verify_hairy_killing_spinor.py`) and the co-author notebooks; the Mathematica suite verifies the displayed reduced jet system and the vacuum problem in three dimensions, and the ten-dimensional *bosonic* equations symbolically (arbitrary chiral L₊ and arbitrary W₀, W₁ on the one-sided family).
* SM (107): the closure is verified in the three-dimensional calibration (real two-component spinors, C = iσ₂); the ten-dimensional bilinear with the full SDFT charge-conjugation convention is not evaluated.
* The fermionic (supercharge) phase-space algebra is not computed anywhere in the archive, exactly as the SM states; the charge algebra verified in NRH05 is the bosonic surface-charge algebra.

* Particular two-point kernels do not fix interior/state/zero-mode completions or local contacts. The matrix in main (21) is an organization of those terms, not a claim that its displayed zeros exhaust the full correlator.
