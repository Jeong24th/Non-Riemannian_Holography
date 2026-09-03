# Equation ledger — every numbered display of the paper, in order

Keyed to the manuscript source with SHA-256
`2E9C2B1DED2D582D44E3A01B55CEDF1045E5E5BCD0EBA55B2464E3E46F41FD41` (2026-09-03, title
"Long-String Holography and Non-Riemannian Hair"), whose isolated build has Letter (1)–(16)
with (4a/4b) and (9a/9b) and SM (1)–(130).  For each display the ledger names the file and the check
label(s) that verify it, or states why no computation applies.

Status codes: **V** = verified by an exact symbolic check; **V°** = verified through its
displayed consequences or at the order stated; **D** = definition/notation (nothing to
compute); **S** = statement whose content is a cited result or an assumption, not an
identity the suite can decide.  Check labels are quoted as printed by the files (`[PASS]`
lines); "→ N" means the label starts with that equation reference.

## Letter

| No. | Label | Content | Where verified | Status |
|---:|---|---|---|:-:|
| 1 | `Rfields` | NS–NS Bañados family (frame ϑ±, metric, B, dilaton) | NRH02 "Eq.(2)…" (built from (1)); the EDFE and horizon checks; NRH06 "SM(74)…E = g − B"; NRH07 "SM(98): R(AdS3, Banados)…" | V |
| 2 | `RDFTfields` | DFT variables H_MN, e^{−2d} of (1) | NRH02 "Eq.(2): H J H = J", "(e^{-2d})^2 = -det g", "e^{-2d} = e^{2y/l}(1 − L+L− e^{−4y/l})"; EDFE "G_MN = 2 l^-2 J_MN" | V |
| 3 | `Rboundary` | asymptotic falloffs | NRH02 "Eq.(3): H − H^infty = O(e^{−2y/l})", "d + y/l = O(e^{−4y/l})" | V |
| 4a | `boundaryH` | H^∞, d^∞ | NRH02 "Eq.(4a) [boundaryH]…" (two checks) | V |
| 4b | `Rboundaryframe` | aligned D = 2 boundary frame | NRH02 "Eq.(4b) [Rboundaryframe]: boundary vielbeins reproduce P^(0) and Pbar^(0)" | V |
| 5 | `Rkilling` | asymptotic-symmetry generator | NRH02 "Eq.(5): Lhat_xi d = O(e^{−4y/l})", "Eq.(5)-(6): Lhat_xi H − delta_eps H = O(e^{−4y/l})"; NRH05 full generators in "SM(49)" | V |
| 6 | `RVirasoro` | δ_εL± with −(l²/4)∂³ε | NRH02 "Eq.(5)-(6)…", "c = 3l/2G reproduces…"; NRH05 "Virasoro cocycle…", "(iii) R bracket − adjoint − central…" | V |
| — | (after 6) | Q = (1/4πGl)∮(L₊ε⁺dx⁺ − L₋ε⁻dx⁻) | NRH05 "SM(49)…" (components) and "(iii) normalization chain (16 pi G)^{-1}(4/l) = 1/(4 pi G l)" | V |
| 7 | `RGKPW` | GKPW relation Z_DFT = Z_CFT | — (defining relation) | D |
| 8 | `Rrenvariation` | δS_ren on shell | NRH04 "SM(11)-(12): delta L_Gamma2 = …" (the exact identity), "SM(11)->(13)…" (on-shell reduction), "SM(13): delta S_ct…"; NRH02 "Eq.(9a) [RKdef]: coefficient matching…" | V |
| — | (after 8) | A^K_{pq̄}, B^K definitions | NRH05 "SM(56) on R/NR: B^y = 4 d_y d"; A^K enters the SM (11) check | V |
| 9a | `RKdef` | ⟨K⟩, ⟨T₍₀₎⟩ from the rescaled momenta | NRH02 "Eq.(9a) [RKdef]: coefficient matching…"; **nonlinear evaluation on both exact saddles**: NRH02 "Eq.(9a) on the exact family: −(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L−},{1, L−}}/(16 pi G l)", "…e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0"; NRH04 "Eq.(9a) on the exact NR family: … = {{L+, W1/4},{0, L−}}/(16 pi G l)", "…hence <T_(0)> = 0"; NRH04 "SM(27)…" | V |
| 9b | `RDFTconservation` | T^CFT_AB and ∇^A T_AB = 0 | NRH02 "Eq.(10): div T = {…} exactly"; "SM(37) [SMinvariantstress]…" | V |
| 10 | `Rcontinuity` | Ward identities | NRH02 "Eq.(10): div T…", "Eq.(10): no local condition on K_{op bom}" | V |
| 11 | `Rcorrelators` | ⟨K⊕⊕̄⟩ = L₊/(16πGl), ⟨KK⟩_c = 3l/(256π²G(x⁺)⁴) | NRH02 "Eq.(11): (64 pi G l)^{-1}(3 l^2/(4 pi)) = (8 pi)^{-2}(c/2)"; NRH04 "SM(28)…stress one-points", "SM(38)…", "SM(39)…", "SM(40)…" | V |
| 12 | `NRvariables` | (Π, q, e^σ, χ) and e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) | NRH03 "Eq.(12): e^{sigma} sinh(chi/2) = …", "SM(62) [NRradialchange]: d chi/d q…" | V |
| 13 | `NRHcompact` | the everywhere non-Riemannian H_MN | NRH03 "Eq.(13): H J H = J…", "type (1,1) at every radius"; EDFE checks; NRH07 "SM(106) -> Eq.(13)…" | V |
| 14 | `NRdilaton` | e^{−2d} = e^{2y/l}(1 − q²), d = −y/l + ln cosh(χ/2√2) | NRH03 "Eq.(14): e^{-2d} e^{-2y/l} = 1 − q^2…"; "SM(61)-(63): d chi/dy…" | V |
| 15 | `NRWgeneral` | exact hair profile W | NRH03 "EDFE scalar: S_(0) = −4/l^2 for ARBITRARY W", "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F", "G_MN = 2 l^-2 J_MN on the ODE shell", "Eq.(15) [NRWgeneral] solves d^2W/dchi^2 = F", "W_0 and W_1 multiply the two homogeneous modes", "near the boundary…", "falloff bookkeeping…" | V |
| 16 | `NRasympt` | δ_εL±, δ_εW₁ at W₀ = 0 | NRH03 "Eq.(16) [NRasympt]…" (four checks) | V |

## Supplemental Material

| No. | Label | Content | Where verified | Status |
|---:|---|---|---|:-:|
| 1 | `SMpert` | H → H + δH, d → d + δd | — | D |
| 2 | `SMprojectors` | P, P̄ | NRH01 (used throughout); NRH02 "Gamma compatibility: nabla_C P_AB = 0" | D |
| 3 | — | h_{pq̄} := δH_MN V^M_p V̄^N_q̄ | — | D |
| 4 | `SMcosetreconstruction` | δH = 2V_(M^p V̄_N)^q̄ h_pq̄ | NRH04 "SM(5)-(6): delta(V eta V^T) = delta H/2…"; the constrained families of NRH04 (H(r,w), H4) | V |
| 5 | — | δV, δV̄ decomposition | NRH04 "SM(5)-(6)…" | V |
| 6 | — | gauge-fixed δV_{Mp} = ½V̄_M^q̄ h_pq̄, δV̄ = −½V h | NRH04 "SM(5)-(6)…" | V |
| 7 | — | δL ≃ 2e^{−2d}(−h K + δd T₍₀₎) | NRH02 "Eq.(9a) [RKdef]: coefficient matching…" (its consequence) | V° |
| 8 | `SMresponsedef` | K_{pq̄}, T₍₀₎ definitions | — | D |
| 9 | `SMKvariation` | K = −e^{2d}(δL/δH)VV̄ = −½e^{2d}δL/δh | — (rewriting of (7)–(8)) | D |
| 10 | `GammaDFT` | L_Γ² = e^{−2d}S₍₀₎ − ∂(e^{−2d}B) = e^{−2d}(PP − P̄P̄)(ΓΓ…) | NRH05 "SM(55) on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)", "SM(55) on NR (arbitrary W)…" (the toolbox Γ² density is the displayed quadratic form) | V |
| 11 | `variation` | exact variation of L_Γ² | NRH04 "SM(11) setup: the BTZ Riemannian double vielbein obeys…", "SM(12) setup: Gamma^L_{LN} = −2 d_N d…", "SM(11)-(12): delta L_Gamma2 = 2e^{-2d}(h S − delta d S_(0)) + d_K(h e^{-2d} A^K + 2 delta d e^{-2d} B^K) exactly, generic tangential h and delta d on BTZ", "SM(11): the identity is not vacuous…" | V |
| 12 | `defB` | B^K and A^K_{pq̄} | as for SM 11; NRH05 "SM(56)…B^y = 4 d_y d" | V |
| 13 | `SMrenvariation` | δ(S_Γ² + S_ct) at the cutoff | NRH04 "SM(11)->(13): on shell…", "SM(13): delta S_ct…" | V |
| 14 | `SMinfinityvielbein` | limiting D = 3 frame | NRH04 "SM(14)-(15): V eta V^T = P^infty…", "V and Vbar are mutually orthogonal"; NRH07 "SM(106) -> SM(14)…" | V |
| 15 | `SMflatmetrics` | η, η̄ | NRH04 "SM(14)-(15)…"; NRH07 "SM(112) rep: {gamma^p, gamma^q} = 2 eta^{pq}" | V |
| 16 | `SMFG` | Fefferman–Graham gauge | imposed in every NRH04 linearization | D |
| 17 | `SMFGcount` | 9 − (3 + 2) = 4 | NRH04 "SM(17): coset count…" | V |
| 18 | `SMfixedprojection` | fixed-frame projection | NRH04 "SM(33)…", "SM(34)…", "SM(41)…", "SM(42)…" | V |
| 19 | `SMNRlin` | linearized EDFE (nine equations) | NRH04 "SM(21) solves the displayed system SM(19)", "SM(19) reproduces the linearized EDFE…" (re-derived from the nonlinear G_MN) | V |
| 20 | `SMNRradialintegration` | elementary radial identities | NRH04 "SM(20): D_y{1, e^{−2y/l}} = 0, …" | V |
| 21 | `SMNRsol` | general solution | NRH04 "SM(21) solves…", "SM(19) reproduces…", "SM(31)-(32)…", "SM(27): d_y delta d = …" | V |
| 22 | `SMcrosspairs` | source–response pairs | NRH04 "SM(25): UV momenta…", "SM(26): the finite piece equals…" (all four pairings) | V |
| 23 | `SMlogfreeconditions` | fixed-dilaton, log-free conditions | NRH04 "SM(23): the log coefficients are exact multiples of the three r^(0) conditions" | V |
| 24 | `SMrwquadratic` | crossed quadratic action | NRH04 "SM(24): Gamma^2 quadratic density = −(1/2) e^{2y/l} dy r dy w" | V |
| 25 | `SMrwvariation` | its UV variation | NRH04 "SM(25): UV momenta give…" | V |
| 26 | `SMSren2` | finite cutoff variation, all four channels | NRH04 "SM(26): …" (six checks: coset family, restricted on-shell sector, vacuum momentum, divergent pieces, finite piece, non-triviality) | V |
| 27 | `SMscalarcutoffresponse` | scalar cutoff coefficient | NRH04 "SM(27) [SMscalarcutoffresponse]: …", "SM(27): d_y delta d = …" | V |
| 28 | `SMresponsematrix` | ⟨K_ab̄⟩ = h⁽²⁾_ab̄/(32πGl) | NRH04 "SM(28): <K_{a bbar}> = h^{(2)}/(32 pi G l) => stress one-points"; normalization from "SM(26)…finite piece" | V |
| 29 | `SMnpointresponse` | connected hierarchy with J^I = −2h^{(0)I} | NRH04 "SM(28)-(30): J^{op bop} = … hence the Hessian factor…" (the canonical-source factor) | D/V° |
| 30 | `SMstresshessian` | Hessian with factor 1/2 | NRH04 "SM(28)-(30)…" | V |
| 31 | `SMcontinuitycheck` | ∂₋⟨K⊕⊕̄⟩ = 0, … | NRH04 "SM(31)-(32) [SMcontinuitycheck]…" | V |
| 32 | `SMtwopointWard` | two-point Ward identities | same check (chirality of the normalizable coefficients) | V° |
| 33 | `SMexactprojectionR` | exact h^R_{pq̄} | NRH04 "SM(33): exact h^R_{p qbar} matches the displayed closed form" | V |
| 34 | `SMexactprojectionRfalloff` | its falloffs | NRH04 "SM(34): falloffs …" | V |
| 35 | — | δd_R = −½ln(1 − L₊L₋e^{−4y/l}) | NRH04 "delta d_R = (1/2) L+ L− e^{−4y/l} + O(e^{−8y/l})" | V |
| 36 | `SMonept` | Riemannian diagonal one-point functions | NRH04 "SM(28): … stress one-points" | V |
| 37 | `SMinvariantstress` | T_{x⁺x̃₋} = 2K⊕⊕̄ | NRH02 "SM(37) [SMinvariantstress]…" (two checks) | V |
| 38 | `SMPBHdata` | s₊, r₊ | NRH04 "SM(38): s_+ = …", "SM(38): r_+ = …", "SM(38): the PBH variation carries no (om bop) source deformation" | V |
| 39 | `SMPBHkernel` | Dirichlet-to-Neumann kernel | NRH04 "SM(39): (l^2/4) d_+^3 [−1/(2 pi x)] = 3 l^2/(4 pi x^4)", "away from coincidence…" | V |
| 40 | `SMRtwopt` | two-point functions, c = 3l/(2G) | NRH04 "SM(40): (64 pi G l)^{-1} 3 l^2/(4 pi) = (8 pi)^{-2}(c/2)"; NRH02 "c = 3l/2G reproduces…" | V |
| 41 | `SMexactprojectionNR` | exact h^NR_{pq̄} | NRH04 "SM(41): exact h^NR = …" | V |
| 42 | `SMexactprojectionNRfalloff` | its falloffs | NRH04 "SM(42): falloffs …" | V |
| 43 | — | δd_NR = ln cosh(χ/2√2) = … | NRH04 "delta d_NR = (1/4) L+ L− e^{−4y/l} + O(e^{−8y/l})" | V |
| 44 | `SMNRonept` | NR one-point functions | NRH04 "SM(44): NR hair channel…", "SM(44): h^(2)_{om bop} = 0 … <K_{om bop}> = 0; SM(45)…", and nonlinearly "Eq.(9a) on the exact NR family: … {{L+, W1/4},{0, L−}}/(16 pi G l)" | V |
| 45 | `SMNRhairhessian` | type-changing Hessian (kernel undetermined) | NRH04 "…SM(45): (1/(64 pi G l)) delta(W_1/2) = (1/(128 pi G l)) delta W_1" (the displayed factor; the kernel is stated to be undetermined) | V° |
| 46 | `SMNRstresskernel` | admissible diagonal kernels vanish in the restricted sector | consequence of SM 21 + 23 (NRH04): the normalizable diagonal coefficients are chiral functions independent of the constant zero-mode source | V° |
| 47 | `SMCPSform` | surface-charge one-form | NRH05 `chargeOneForm` (all three terms), "SM(50)…", "(i) integrability…" | V |
| 48 | — | Q_R | NRH05 "SM(49)…" with "(iii) normalization chain…" | V |
| 49 | `SMRpotentialcomponents` | lim e^{−2d}K̂ components | NRH05 "SM(49): lim e^{-2d} Khat^{-y}[eps+] = …", "…Khat^{y+}[eps-] = … (mirror)" | V |
| 50 | `SMCPSresult` | Θ̂ → 0, k = (4/l)εδL | NRH05 "SM(50): lim e^{-2d} Thetahat^{+,-,y} = 0", "SM(50): k^{-y}[eps+] = (4/l) eps+ dL+", "…k^{+y}[eps-]…" | V |
| 51 | `SMNRchargefalloffs` | state-dependent falloffs | NRH05 "SM(51): state-dependent falloffs…" | V |
| 52 | `SMNRchargecancellation` | componentwise W₁ cancellation | NRH05 "SM(52): W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise" | V |
| 53 | `SMNRcharge` | δQ, Q | NRH05 "(i) k^{-y}[eps+] = delta[(4/l) eps+ L+]", "(i) mirror…" | V |
| 54 | `SMCPSalgebra` | {Q,Q} = Q[[ε,η]], c_charge = 0 (+ footnote) | NRH05 "SM(54): …" (four checks), "(ii)…", "(iv)…", "(v)…", "(vi)…", "SM(54) footnote: …" | V |
| 55 | `SMgamma2` | Γ² identity | NRH05 "SM(55) on R…", "SM(55) on NR (arbitrary W)…" | V |
| 56 | `SMgamma2flux` | B^y = 4∂_y d, e^{−2d}B^y = −(4/l)(e^{2y/l} + μe^{−2y/l}) | NRH05 "SM(56) on R…", "SM(56) on NR…" (three checks) | V |
| 57 | `SMmudefinition` | μ dichotomy | NRH05 "SM(56) on NR: −2 d_y e^{-2d} = … [mu-dichotomy]" | V |
| 58 | `SMgamma2cutoff` | S_ren(Y) | NRH05 "SM(58): the regulated combination equals…" | V |
| 59 | `SMgamma2value` | S_ren = −8√μ/(16πGl)∫d²x | NRH05 "SM(59): Y -> Infinity limit gives…" | V |
| 60 | — | Killing horizon e^{4y/l} = L₊L₋ | NRH02 "horizon: e^{-2d} = 0 at e^{4y/l} = L+ L−"; NRH05 "endpoints…" | V |
| 61 | `NRhill` | ψ±, A± = ¼(∂lnL)² − ½∂²lnL = ψ″/ψ | NRH03 "SM(61) [NRhill]: …"; A± enters the verified source F | V |
| 62 | `NRradialchange` | q, χ, ∂_y, ∂_q | NRH03 "SM(62) [NRradialchange]: d chi/d q = …", "SM(61)-(63): d chi/dy = …" | V |
| 63 | `NRradialoperator` | radial operator identity | NRH03 "SM(63) [NRradialoperator]…" | V |
| 64 | `NRchiODE` | ∂²_χ W = F | NRH03 "EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F" | V |
| 65 | `NRsource` | the source F | same check (F as displayed) | V |
| 66 | `NRg` | ρ(χ) | NRH03 "SM(66)-(67) [NRg, NRGprofile]: d^2 G/d chi^2 = rho(chi)" | V |
| 67 | `NRGprofile` | G(χ) normalization | NRH03 "SM(67) [NRGprofile]: G'(0) = 0", "…(e^{s chi} − 1 − s chi)…" | V |
| 68 | — | q = √(L₊L₋/2) μ_RG^{−2}, χ, μ_RG = e^{y/l} | — (definitions; used in the next line) | D |
| 69 | — | μ_RG dχ/dμ_RG = −2√2 sinh(χ/√2) | NRH03 "SM: RG rapidity mu d chi/d mu = −2 Sqrt[2] Sinh[chi/Sqrt[2]]" | V |
| 70 | `SMBtransform` | finite B-transformation | NRH03 "SM(70)-(71) [SMBtransform, SMWshift]…", "Omega_b is O(3,3)" | V |
| 71 | `SMWshift` | Ω_bH(W)Ω_bᵀ = H(W − 2b_{+−}) | NRH03 "SM(70)-(71)…" | V |
| 72 | `SMdyg` | doubled-yet-gauged action | — (defining action; its reductions are verified below) | D |
| 73 | `SMphysicalsectionA` | Â_αμ, D_αx^M on the section | — | D |
| 74 | `SMRfirstorder` | Riemannian first-order form | NRH06 "SM(74): the auxiliary equations give…", "SM(74): eliminating the auxiliaries reproduces E_{mu nu}…" | V |
| 75 | `SMceff` | c_eff² = 2F | NRH06 "SM(75): −det g_par = …", "SM(75): c_eff^2 = 2F -> 4 Sqrt[L+L−]…"; NRH02 "c_eff^2 := 2F…" | V |
| 76 | `SMlongstringE` | winding-string energy of the static probe (SM 5.1) | NRH06 "before SM(76): … d_t is Killing…", "before SM(76): gamma^{tau a} g_{t nu} d_a X^nu = … = 1", "before SM(76): B_{t phi} = …", "before SM(76) -> SM(76): E = −Int d sigma P_t … reproduces E(y) = …", "SM(76): −det g_(t,phi) = l^2 (e^{2y/l} − L+L− e^{−2y/l})^2 per winding…", "SM(76): E(y) = …", "SM(76) remark: with phi_0 = 0 the area density equals l e^{-2d}…" | V |
| 77 | `SNCtau` | SNC clock forms | NRH03 "SM(77) [SNCtau]: tau+ = dx+ − L_− e^{−2y/l} dx− + O(e^{−4y/l})", "unit clock determinant", "W-part of the lower-right block = W(tau+ tau− + tau− tau+)"; NRH07 "SM(106): x rows = Sqrt[2] tau^pm…" | V |
| 78 | `SMdygconstraints` | τ⁺·∂̄x = 0 = τ⁻·∂x | NRH06 "SM: H^{mu nu} tau^pm_nu = 0", "SM: the dual vectors Y, Ybar exist at every radius" (the kernel and multiplier structure) | V° |
| 79 | `SMdygGO` | exact reduced Lagrangian | NRH06 "SM(83): symmetric-block route − antisymmetric-clock route = W (tau− . dx)(tau+ . dbar x)", "SM(80)-(81): at chi -> 0 the W coupling reduces to…" | V° |
| 80 | `SMGO` | Gomis–Ooguri limit | NRH06 "SM(80)-(81)…"; contraction bookkeeping | V |
| 81 | `SMvertex` | V_W = (1/4πα′)W∂x⁺∂̄x⁻ | NRH06 "SM(81): with the 1/(2 pi alpha') prefactor this is V_W = …" | V |
| 82 | `SMWisB` | H(W) = Ω_bH(0)Ω_bᵀ, b_{+−} = −W/2 | NRH03 "SM(82) [SMWisB]…" | V |
| 83 | `SMdeltaL` | antisymmetric clock coupling | NRH06 "SM(83): …" | V |
| 84 | `SMFTsector` | S_y, S_FT | — (definitions; used in SM 86) | D |
| 85 | `SMWZWlevel` | k = l²/α′ | NRH06 "SM(85): k = |Int_{S^3} H| / (4 pi^2 alpha') = l^2/alpha'" | V |
| 86 | `SMFTweight` | T_y, h_y(a) | NRH06 "SM(86): the two OPE contributions assemble to h_y(a)/(z−w)^2" | V |
| 87 | `SMBRSTradial` | (∂²_y + (2/l)∂_y)f = 0 | NRH06 "SM(87): the marginal roots…", "SM(87): (d_y^2 + (2/l) d_y) f = 0…" | V |
| 88 | `SMBRSTmomentum` | h_y(k), roots k = 0, 2i/l | NRH06 "SM(88): h_y(i k) = … P = k − i/l…" | V |
| 89 | `SMBRSTvertex` | U_{W₁} is a (1,1) primary | SM 86–88 plus NRH06 "<x^+ x^−> = 0…" (W₁ has weight (0,0)) | V° |
| 90 | `SMBRSTcentral` | c_y = 1 + 6α′/l², total 3(k+2)/k | NRH06 "SM(90): c_y/2 = …", "SM(90): c_{beta gamma} + c_y = …" | V |
| 91 | `SMWgaugeobstruction` | W₀ gauge, e^{−2y/l}W₁ not | NRH06 "SM(91): Lhat_xi H^infty = D H^infty + H^infty D^T…", "SM(91): the conditions force…", "SM(91): (db)_{+−y} = …", "SM(91): v^y = … so closure forces d_y w = 0" | V |
| 92 | `SMBRSTfusion` | self-contraction, h_n, resonances | NRH06 "SM(92): …" (three checks) | V |
| 93 | `SMupliftblocks` | H₁₀ = H₃ ⊕ H_{S³} ⊕ H_{R⁴}, d₁₀ = d₃ + d_{S³} | NRH07 `d10Assemble` and the two ten-dimensional probes "S_(0)^{(10)} = 0", "(P S Pbar)^{(10)} = 0" | V |
| 94 | `SMDFTKilling` | generalized Killing equations | — (definition; `GenLieH`/`GenLieD` implement it: NRH07 "SM(109)…", "SM(110)…") | D |
| 95 | `SMtypeIIKS` | type-II Killing-spinor systems | — (definition of the systems solved in NRH07) | D |
| 96 | `SMsusyclosure` | closure on ĥL_X, X = iε̄₂Γε₁ | NRH07 "C = i sigma_2 is the Majorana conjugation…", "SM(96)->SM(109): X^M = …", "SM(96): Lhat_X H^infty = 0 and Lhat_X d = 0…", "SM(96): the bilinear is symmetric…" (three-dimensional calibration) | V° |
| 97 | `SMsemicov` | semi-covariant connection | NRH01 `GammaDFT`; NRH02 "Gamma compatibility: nabla_C P_AB = 0", "Gamma dilaton trace…", "Gamma torsionless…" | V |
| 98 | `SMuplift` | AdS₃×S³×R⁴ with flux | NRH07 "SM(98): R(S3) = +6/l^2 and H^2(S3) = +24/l^2", "SM(98): R_{mu nu} = (1/4) H H on the S3 factor", "SM(98): R(AdS3, Banados) = −6/l^2 and H^2 = −24/l^2…", "SM(98): R_{mu nu} = (1/4) H H on the AdS3 factor", "S_(0)(S3 block) = +4/l^2", "S_(0)(R4 block) = 0" | V |
| 99 | `SMRDFTKilling` | ordinary-field form of (94) | — (decomposition statement; the six Killing vectors are not displayed) | S |
| 100 | `SMRlocaliso` | k^{(ij)} = s_is_j stabilizers | NRH07 "SM(100): k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' − (l^2/4) k''' = 0" | V |
| 101 | — | vol₃, H₃ = −(2/l)vol₃, H_{y−+} = +(2/l)√|g₃| | NRH07 "orientation before SM(103): H_{y−+} = +(2/l) Sqrt[|g_3|]…" | V |
| 102 | — | ∇_με± = ±(1/2l)γ_με± | — (the torsionful Killing-spinor equation of the Riemannian branch; its consequence SM 103–104 is checked) | S |
| 103 | `SMRcomponentHill` | first-order pair → Hill equation | NRH07 "SM(103)-(104): the first-order pair … closes into (l^2/2) s'' = L s" | V |
| 104 | `SMRlocalKS` | Hill equations for s± | NRH07 "SM(103)-(104)…", "global AdS3 (L = −1/4)…", "massless BTZ (L = 0)…", "constant L > 0…" | V |
| 105 | — | N_local = 16 = 8 + 8 | NRH07 "counting: (4,4) + (4,4) = the sixteen constant vacuum modes…" | V° |
| 106 | `SMvielbein` | exact non-Riemannian double vielbein | NRH07 "SM(107): V_M^p V_Np = P_MN and …", "SM(107): V^M_p Vbar_{M qbar} = 0 and P + Pbar = J", "SM(106) -> Eq.(13)…", "SM(106) -> SM(14)…", "SM(106): x rows = Sqrt[2] tau^pm…" | V |
| 107 | `SMvielbeincheck` | defining relations | same checks | V |
| 108 | `SMNRlocalstabilizer` | local stabilizer system | NRH07 "SM(108): eps = c/Sqrt[L] solves…", "SM(108): exact one-sided identity Lhat_xi H = delta H, with weights (1,2)"; NRH03 "SM(108) line 2: … delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i" | V |
| 109 | `SMexactiso` | vacuum isometries | NRH07 "SM(109): Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm"; "SM(96)->SM(109)…" | V |
| 110 | `SMweighteddilaton` | weighted dilaton condition | NRH07 "SM(110): Lhat_xi d = 0…", "SM(110): equivalently d_M(e^{-2d} xi^M) = 0" | V |
| 111 | `SMkillingspinor` | vacuum Killing spinor E = (√2f, l∂₊f) | NRH07 "SM: vacuum spin connection = displayed…", "SM(111): D_{pbar} E = 0…", "SM(111): gamma^p D_p E = E/(Sqrt[2] l)" | V |
| 112 | `SMcomplexblocks` | σ_i, τ_i, ρ_m | NRH07 (matrices used in SM 113–115 checks); "SM(112) rep: {gamma^p, gamma^q} = 2 eta^{pq}" | V |
| 113 | `SMgammaten` | ten-dimensional Γ^p̂, Γ₁₁ | NRH07 "SM(113): {Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32", "SM(113): Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0" | V |
| 114 | `SMgammaMajorana` | Majorana intertwiner | NRH07 "SM(114): BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, …" | V |
| 115 | `SMgammabarred` | barred Clifford algebra | NRH07 "SM(115): {Gammabar, Gammabar} = −2 eta, Gammabar_{pq} = −Gamma_{pq}, same Majorana intertwiner" | V |
| 116 | `SMreducedDirac` | reduced two-component system | NRH07 "SM(116): reduced system…", "SM(116): the opposite channel E = (0, g(x+)) has eigenvalue −1/(Sqrt[2] l)", "SM(116): both channels also satisfy D_{pbar} E = 0" | V |
| 117 | `SMinternalprojectors` | ζ±, κ±, Ξ_{+r}, Ξ′_{−r}; S³ spinor equations | NRH07 "SM(118): Weyl + S^3-line + zeta_+ leave complex dimension 4 (the Xi_{+r} span)", "below SM(116): … integrable (flat connection)" | V° |
| 118 | `SMspinorcountchain` | 32_C → 32_R → 16_R → 4_R | NRH07 "SM(118): the Majorana condition leaves 32 real components", "…Weyl condition leaves 16…", "…complex dimension 4…" | V |
| 119 | `SMhairyKS` | one-sided Killing-spinor families | NRH07 "SM(119): the displayed system has rank six…", "SM(119): the solution space is exactly {e_1, d_+ e_1}…", "SM(119): the mechanism…" (the displayed jet system; the ten-dimensional evaluation quoted in the SM is the Python/SymPy guard and the co-author notebooks) | V° |
| 120 | `SMcandidateD` | gauge derivative 𝔻 | — (definition; used in NRH08) | D |
| 121 | `SMcandidateaction` | candidate action and its reduction | NRH08 "SM(121): the P^(0) term contains only D_+…", "SM(121): each term selects a single independent component of phi_A", "SM(121): gamma^oplus … rank one" | V |
| 122 | `SMcandidategauge` | internal gauge invariance | NRH08 "SM(122): delta_Lambda L = 0 for the non-abelian bosonic sector" | V |
| 123 | `SMcandidateWitt` | Witt transformations | NRH08 "SM(124): delta_v L = d_+(v^+ L) + d_−(v^− L)"; "Lhat_xi H^(0) = 0 <=> d_− v^+ = 0 = d_+ v^−" | V |
| 124 | `SMcandidateWittL` | δ_vL = ∂(vL) | same | V |
| 125 | `SMcandidatefermionic` | chiral fermionic transformations | NRH08 "SM(126): delta_{eps+} L = …", "SM(126): delta_{eps−} L = …", "the fermionic parameters carry arbitrary chiral profiles" | V |
| 126 | `SMcandidatefermionicL` | their total-derivative variations | same | V |
| 127 | `SMcandidateextra` | Grassmann-even chiral transformation | NRH08 "SM(128): delta_zeta L = …" | V |
| 128 | `SMcandidateextraL` | its variation | same | V |
| 129 | `SMcandidatetrivial` | equation-of-motion redundancies | NRH08 "SM(129): the transformations are built from the fermion equations of motion…" | V |
| 130 | `SMcandidatetrivialL` | δ_αL total derivative | NRH08 "SM(130): delta_alpha L = d_−(psi^+ delta_alpha psi^+) + d_+(psi^− delta_alpha psi^−) identically (off shell)", "the two-component fermions used here are genuinely Grassmann…" | V |

## Prose statements verified in addition to the numbered displays

* Letter, below Eq. (12): the expansion e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) (NRH03).
* Letter, below Eq. (15): the boundary behaviour of W (exact in the one-sided limits; O(e^{−6y/l}) for constant L±) via the χ-expansion of the two homogeneous modes and the orders of G(χ) and e^{sχ} − 1 − sχ (NRH03).
* Letter, "Radial hair": the B-shift with 2∂_[+λ_−] = W₀/2 removes W₀ (NRH03, SM 70–71 checks).
* Letter, "Riemannian Reference": the BTZ horizon e^{2y/l} = √(L₊L₋) and c_eff (NRH02).
* SM 1: the vacuum has no first-order boundary momentum (no vacuum one-point function) (NRH04, SM (26) block).
* SM 3: the endpoints e^{−2d} = 0 for both branches (NRH05).
* SM 5.1: E_{+−} = 0, E_{−+} = −2F, E_{±±} = 2L± (NRH02, NRH06); −det g_∥ = F² − 4L₊L₋ = e^{−4d_R} (NRH06).
* SM 5.1, before (76): the unnumbered Noether-charge display P_μ = (2πα′)⁻¹(−√−γ γ^{τa}g_{μν}∂_aX^ν + B_{μν}X′^ν), E = −∫dσ P_t = (2πα′)⁻¹∫dσ(√−γ − wB_{tφ}), the identity γ^{τa}γ_{aτ} = 1, the Killing property of ∂_t, B_{tφ} = l(e^{2y/l} + L₊L₋e^{−2y/l}), −det g_(t,φ) = l²(e^{2y/l} − L₊L₋e^{−2y/l})², and the reproduction of (76) from the displayed momentum (NRH06, four checks).
* SM 5.2: the Gomis–Ooguri contraction table behind the exact marginality of V_W[W₀] (NRH06).
* SM 6: the two ten-dimensional exact probes of S₍₀₎⁽¹⁰⁾ = 0 = (PSP̄)⁽¹⁰⁾ for both uplifts (NRH07); the displayed vacuum spin-connection components; the Hill-monodromy counts; the complementary-halves bookkeeping.
* SM 6, text below (106): e^{±σ}sinh(χ/2) = uL±[1 + u²Π/3 + …] (NRH03).

## What the suite does not decide

* The SM (11) identity is checked on the exact BTZ background (constant L±) in its own Riemannian double-vielbein frame, for the general tangential fluctuation h_{ab̄}(x⁺,x⁻,y) of the Fefferman–Graham gauge and an arbitrary δd(x⁺,x⁻,y); the radial fluctuation components h_{aȳ}, h_{yb̄}, h_{yȳ} are not switched on.  It is a general identity of the semi-covariant formalism (Ref. Park:2025core), and its consequences (13), (26)–(28), and the nonlinear Eq. (9a) evaluations are verified independently.
* SM (99) and SM (102) are decompositions of DFT statements into ordinary Riemannian ones; the ordinary Killing vectors are not displayed and were not enumerated here.
* SM (119): the full ten-dimensional Killing-spinor evaluation with hair is quoted in the SM from the Python/SymPy guards (`checks/verify_10d_killing_spinor.py`, `verify_hairy_killing_spinor.py`) and the co-author notebooks; the Mathematica suite verifies the displayed reduced jet system and the vacuum problem in three dimensions, and the ten-dimensional *bosonic* equations at exact rational probes.
* SM (96): the closure is verified in the three-dimensional calibration (real two-component spinors, C = iσ₂); the ten-dimensional bilinear with the full SDFT charge-conjugation convention is not evaluated.
* The fermionic (supercharge) phase-space algebra is not computed anywhere in the archive, exactly as the SM states; the charge algebra verified in NRH05 is the bosonic surface-charge algebra.
