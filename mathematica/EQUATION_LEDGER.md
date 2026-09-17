# Equation ledger

Current source SHA-256: `09B04C8BCF63EBDB4A879E1C893D28F43455756535340161D1138123B55DDA84` (2026-09-17).

Letter (1)-(22); Supplemental Material (1)-(205); 231 numbered displays.
Status counts: V 184, V° 14, D 25, S 5, N 3.

Every numbered display of the current source is listed below in manuscript order. Coverage
refers to the public Mathematica suite in this directory; the quoted check names are the
descriptions printed by the files at run time (the branch iterator is expanded to R or NR), so
each entry can be located with a text search in the `.wl` sources or in the run log.

**V**: symbolic check of the stated relation, for arbitrary chiral L± and, on the NR branch,
arbitrary W₁(x⁺,x⁻), unless the check name says otherwise; **V°**: only a stated consequence,
background value or asymptotic order is checked; **D**: definition or notation; **S**: cited or
stated without a machine check; **N**: no current public check. The execution record is in
[REFERENCE_RUN.md](REFERENCE_RUN.md).

The methods follow the manuscript. The near-boundary solutions are obtained by solving the
coupled hierarchy of the linearized EDFE order by order (SM 3.4), with the undetermined
responses R±, H_s and the zero mode c_s carried as free data. The two-point functions are
computed from the ordered second variation of the on-shell action with two independent bulk
solutions inserted in the two slots (SM 3.8); they are never obtained by differentiating
one-point functions. The counterterm cancellation is tested modulo total tangential
derivatives with the stress constraints imposed in both slots. The nonlocal kernels follow
from the Ward operators and the manuscript's Lorentzian inverse-derivative convention.

Not covered (**N**/**S**): the general variation formulas of SM 3.5 (the linearized EDFE are
obtained by direct linearization of the exact curvature instead), the symmetry statement of the
action-response matrix, and the ten-dimensional Killing-spinor statements marked **S**. A
passing reduced jet system is not a construction of nonzero supercharges.


## Letter

| No. | Label | Content | Public coverage | Status |
|---:|---|---|---|:---:|
| 1 | `Rfields` | NS–NS Bañados family (frame ϑ±, metric, B, dilaton) | NRH02 "Rfields: dy^2 - 2 theta+ theta- reproduces the displayed metric components", "Rfields: B_{-+} = e^{2y/l} + e^{-2y/l} L+ L- (B = (...) dx^- ^ dx^+)" | V |
| 2 | `RDFTfields` | DFT variables H_MN and e^{-2d} for the Riemannian family | NRH02 "RDFTfields: H J H = J (O(3,3) constraint)", "RDFTfields: e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l}) = sqrt(-g) at phi_0 = 0" | V |
| 3 | `Rboundary` | asymptotic falloffs | NRH02 "Rboundary: H - H^infty = O(e^{-2y/l}) and d + y/l = O(e^{-4y/l})" | V |
| 4a | `boundaryH` | H^∞, d^∞ | NRH02 "boundaryH: H^infty is the displayed 6x6 matrix", "boundaryH: the induced boundary H^(0) is of type (1,1): vanishing upper-left block" | V |
| 4b | `Rboundaryframe` | aligned D = 2 boundary frame | NRH02 "Rboundaryframe: the aligned D = 2 double vielbein reproduces P^(0), Pbar^(0)" | V |
| 5 | `NRvariables` | (Π, q, e^σ, χ) and e^{±σ}sinh(χ/2) = e^{−2y/l}L± + O(e^{−6y/l}) | NRH02 "NRvariables: chi = 2 Sqrt[2] arctanh q and e^{pm sigma} sinh(chi/2) = e^{-2y/l} L_pm + O(e^{-6y/l})" | V |
| 6 | `NRHcompact` | the everywhere non-Riemannian H_MN | NRH02 "NRHcompact: type (1,1) at every radius (upper-left block diag(0,0,1))", "NRHcompact: H J H = J for the exact non-Riemannian matrix (generic W)" | V |
| 7 | `NRdilaton` | e^{−2d} = e^{2y/l}(1 − q²), d = −y/l + ln cosh(χ/2√2) | NRH02 "NRdilaton: e^{-2d} = e^{2y/l}(1 - q^2) with q = tanh(chi/(2 Sqrt[2]))" | V |
| 8 | `NRWgeneral` | exact hair profile W | NRH02 "NRWgeneral: W_0 and W_1 multiply the two homogeneous modes {1, chi/(2 Sqrt[Pi])}", "NRWgeneral: solves d^2W/dchi^2 = F with G'' = rho"; NRH04 "NRWgeneral: twice integrating NRchiODE gives the exact solution (with G'' = rho)"; NRH05 "SMbackgroundexpansion: W_NR = W_1 u - (l^2/4) L+' L-' u^2 + O(u^3) from the exact NRWgeneral" | V |
| 9 | `MainGprofile` | radial profile G(χ) = 4√2∫₀^χ[sinh t/sinh(t/√2) − √2]dt, G(0) = G′(0) = 0 | NRH02 "MainGprofile: G(0) = G'(0) = 0 and G'' = rho: G'(chi) = (2/3) chi^2 + O(chi^4), hence G = (2/9) chi^3 + O(c…" | V |
| 10 | `RGKPW` | GKPW relation Z_DFT = Z_CFT | — (defining relation) | D |
| 11 | `Rrenvariation` | On-shell boundary variation of the renormalized action | NRH06 "SMrenvariation, SMvolumecounterterm: on shell the bulk term cancels delta(-2 Lambda e^{-2d}) and the counte…" | V° |
| 12 | `Mainmomenta` | Generalized-metric and scalar radial momenta A and B | NRH06 "SMmetricradialdictionary, SMscalarradialdictionary: coefficient matching -2K = (16 pi G)^{-1} A^y and 2 T_(…", "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" (+1 more) | V° |
| 13a | `RKdef` | ⟨K⟩, ⟨T₍₀₎⟩ from the rescaled momenta | NRH02 "RKdef: -2K = (16 pi G)^{-1} A^y and 2 T_(0) = (16 pi G)^{-1} 2 (B^y + 4/l) fix the coefficients -1/(32 pi G…" | V |
| 13b | `RDFTconservation` | T^CFT_AB and ∇^A T_AB = 0 | NRH02 "RDFTconservation, Rcontinuity: div T = {d_- K_mp, -d_+ K_mp, -2 Ward_+, -2 Ward_-}" | V |
| 14 | `Rcontinuity` | Ward identities | NRH02 "RDFTconservation, Rcontinuity: div T = {d_- K_mp, -d_+ K_mp, -2 Ward_+, -2 Ward_-}", "Rcontinuity: no local condition on K_{op bom}; K_{om bop} must be constant" | V |
| 15 | `Rkilling` | common R/NR asymptotic-symmetry generator (ξ^±_R with the l²e^{−2y/l}∂²ε/4 tail, ξ^±_NR = ε^±; common radial and dual components) | NRH02 "Rkilling, RVirasoro: Lhat_xi H - delta_eps H = O(e^{-4y/l}) with delta L = eps L' + 2 L eps' - (l^2/4) eps'''", "Rkilling: Lhat_xi d = O(e^{-4y/l}) (the radial component compensates the weight)" (+1 more) | V |
| 16 | `RVirasoro` | δ_εL± with −(l²/4)∂³ε | NRH02 "Rkilling, RVirasoro: Lhat_xi H - delta_eps H = O(e^{-4y/l}) with delta L = eps L' + 2 L eps' - (l^2/4) eps'''", "charges: with T = 8 pi K = L+/(2 G l), RVirasoro is delta T = eps T' + 2 T eps' - (c/12) eps''' with c = 3l…" | V |
| 17 | `NRasympt` | δ_εL±, δ_εW₁ at W₀ = 0 | NRH02 "Rkilling (NR form), NRasympt: Lhat_xi H - delta_(L, W1) H = O(u^-2) componentwise", "NRasympt: delta psi is equivalent to delta L = eps L' + 2 L eps' (no central term)" (+1 more) | V |
| 18 | `Mainonepoints` | R/NR one-point matrices and vanishing scalar responses | NRH02 "Mainonepoints on R: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) (exact fam…", "Mainonepoints on R: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0" (+2 more) | V |
| 19 | `Mainframevariation` | Frame variations in the gauge with no local Lorentz rotation | NRH02 "Mainframevariation: delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2" | V |
| 20 | `Mainfivebyfive` | connected 5x5 matrix: particular A/M entries + R_s + local contacts | NRH02 "Mainfivebyfive: M_+^R = L_- A_+^R" | V |
| 21 | `Rcorrelators` | particular A_NR, A_R, M_NR, M_R kernels, with coefficients at x | NRH02 "Rcorrelators: A_+^NR = [L+'/D+ - 2 L+/D+^2]/(128 pi^2 G l)", "Rcorrelators: A_+^R = A_+^NR + 3 l/(256 pi^2 G D+^4)" (+1 more); NRH06 "SMgeneralpositionkernels, Rcorrelators: A_+^NR = [L+'/D+ - 2L+/D+^2]/(128 pi^2 G l)", "SMgeneralpositionkernels, Rcorrelators: A_+^R = A_+^NR + 3l/(256 pi^2 G D+^4)" (+1 more) | V |
| 22 | `Mainworldsheet` | ℒ_GO = (1/2πα′)(β∂̄x⁺ + β̄∂x⁻ + ∂y∂̄y), V_𝒲 = (𝒲/4πα′)∂x⁺∂̄x⁻ | NRH02 "Mainworldsheet: eliminating beta, betabar from the first-order action reproduces E_{mu nu} dx dbar x (E = g…", "Mainworldsheet: V_W = (1/(4 pi alpha')) W dx+ dbar x- with the 1/(2 pi alpha') prefactor of the reduced action" | V |


## Supplemental Material

| No. | Label | Content | Public coverage | Status |
|---:|---|---|---|:---:|
| SM 1 | `SMgamma2` | Γ² identity | NRH03 "SMgamma2 on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M) for arbitrary chiral L_pm", "SMgamma2 on R: the section fluxes B^{x pm} are built from the W-independent blocks and are total x-derivati…" (+2 more) | V |
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
| SM 12 | `NRg` | ρ(χ) | NRH04 "NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)" | V |
| SM 13 | `NRGprofile` | G(χ) normalization | NRH04 "NRGprofile: G(0) = G'(0) = 0: G'(chi) = (2/3) chi^2 + O(chi^4), hence G = (2/9) chi^3 + O(chi^5)", "NRGprofile: d^2/dchi^2 (e^{s chi} - 1 - s chi) = e^{s chi}" (+1 more) | V |
| SM 14 | — | constant-L radial variables q = √(L₊L₋/2) μ_RG^{−2}, χ = 2√2 arctanh q, μ_RG = e^{y/l} | definition; used by the NRH04 check of SM15 | D |
| SM 15 | — | μ_RG dχ/dμ_RG = −2√2 sinh(χ/√2) | NRH04 "constant L: q = Sqrt[L+L-/2] mu_RG^{-2}; mu d chi/d mu = -2 Sqrt[2] sinh(chi/Sqrt[2])" | V |
| SM 16 | `SMBtransform` | finite B-transformation | NRH04 "SMBtransform: Omega_b is O(3,3)", "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| SM 17 | `SMWshift` | Ω_bH(W)Ω_bᵀ = H(W − 2b_{+−}) | NRH04 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched"; NRH05 "SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly" | V |
| SM 18 | `SMprojectors` | P, P̄ | NRH05 "SMprojectors: P + Pbar = J, P J P = P, Pbar J Pbar = Pbar, P J Pbar = 0 on the limiting metric" | V |
| SM 19 | `SMmixedfluctuation` | h_{pq̄} := δH_MN V^M_p V̄^N_q̄ | NRH05 "SMmixedfluctuation: projecting the reconstruction back gives h_{p qbar}" | V |
| SM 20 | `SMcosetreconstruction` | δH = 2V_(M^p V̄_N)^q̄ h_pq̄ | NRH05 "SMcosetreconstruction: P delta H P = 0 = Pbar delta H Pbar for the reconstructed mixed variation" | V |
| SM 21 | `SMframevariation` | gauge-fixed δV_{Mp} = ½V̄_M^q̄ h_pq̄, δV̄ = −½V h | NRH05 "SMframevariation: delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2" | V |
| SM 22 | `SMresponsedef` | K_{pq̄}, T₍₀₎ definitions | — | D |
| SM 23 | `SMbackgroundconnection` | Torsionless DFT connection at D=3 | NRH05 "SMbackgroundconnection on R: nabla_C P_AB = 0 through z^2", "SMbackgroundconnection on NR: nabla_C P_AB = 0 through z^2" (+4 more) | V |
| SM 24 | `SMconnectionvariation` | δΓ_KMN in terms of δP and δd (Eq. (2.56) of Park:2025core) | — (quoted from Ref. Park:2025core; not re-derived in this archive) | S |
| SM 25 | `SMsixprojectors` | six-index projectors 𝒫, 𝒫̄ (Eq. (2.52) of Park:2025core) | — (definition) | D |
| SM 26 | `SMflatmetrics` | η, η̄ | NRH05 "SMflatmetrics, SMinfinityvielbein: V eta V^T = P^infty, Vbar etabar Vbar^T = Pbar^infty, V^T J Vbar = 0" | V |
| SM 27 | `SMinfinityvielbein` | limiting D = 3 frame | NRH05 "SMflatmetrics, SMinfinityvielbein: V eta V^T = P^infty, Vbar etabar Vbar^T = Pbar^infty, V^T J Vbar = 0", "both saddle frames tend to the limiting frame SMinfinityvielbein at u = 0" | V |
| SM 28 | `SMFG` | Fefferman–Graham gauge | imposed in every NRH04 linearization | D |
| SM 29 | `SMsources` | Four crossed lower-index source couplings | NRH05 "SMsources: -h^{(0)a bbar} K_{a bbar} = h_mm K_pp + h_pp K_mm + h_mp K_pm + h_pm K_mp" | V |
| SM 30 | `SMW0completion` | Transport of fields and frames by the closed W0 B shift | NRH05 "SMW0completion: b_{+-} = -W_0/2 gives H_s[W_0] with W -> W + W_0, and Omega_b is O(3,3)", "SMW0completion: the transported frames V_s[W_0] = Omega_b V_s[0] reconstruct P, Pbar of H_s[W_0]" (+1 more) | V |
| SM 31 | `SMW0variation` | Product rule when the W0 source varies | NRH05 "SMW0variation: varying W_0 gives h^{(0)}_{op bom} = delta W_0/2 at the reference boundary" | V |
| SM 32 | `SMbackgroundexpansion` | General R/NR dilatons and derivative-dependent W expansion | NRH05 "SMbackgroundexpansion: d_R = -y/l - (1/2) ln(1 - Pi u^2) and d_NR = -y/l - (1/2) ln(1 - Pi u^2/2)", "SMbackgroundexpansion: e^{-2 d_s} = u^{-1}[1 + O(u^2)] on both saddles" (+1 more) | V |
| SM 33 | `SMexactboxEDFE` | Exact linearized scalar and mixed EDFE via the universal box | NRH05 "SMexactnormalform: M = diag(-1/4,-1/4,-1/4,-1/4,1) multiplies d_y^2 (h_pp,h_pm,h_mp,h_mm,delta d) in (E_pp,…", "SMEzerocomponents: the ten displayed leading operators E^{(0)} equal the u^0 part of the linearized EDFE on…" | V° |
| SM 34 | `SMboxdefinition` | Universal box and section-condition identities | — | D |
| SM 35 | `SMmixedboxdefinition` | Full mixed-tensor box with curvature and connection terms | — | D |
| SM 36 | — | □δd = H^{MN}∇_M∇_N δd = e^{2d}∂_M(e^{−2d}H^{MN}∂_N δd) (dilaton Laplacian identity) | identity used in the SM 3.5 variation formulas; not separately checked | D |
| SM 37 | `SMboxcurvature` | Connection field strength and its Ricci contraction | Definition. | D |
| SM 38 | `SMexactRdata` | Exact Riemannian fields and lower-index saddle frames | NRH05 "SMexactRdata: e eta e^T = g and ebar etabar ebar^T = -g; V, Vbar reconstruct P, Pbar and are orthogonal" | V |
| SM 39 | `SMexactNRdata` | Exact NR substitution with full derivative-dependent W | NRH05 "SMexactNRdata: the lowered SMvielbein frames reconstruct P, Pbar of the exact hairy metric and are orthogonal" | V |
| SM 40 | `SMexactcomponentrecipe` | Definition of the ten component equations | NRH05 "SMexactnormalform: M = diag(-1/4,-1/4,-1/4,-1/4,1) multiplies d_y^2 (h_pp,h_pm,h_mp,h_mm,delta d) in (E_pp,…", "SMEzerocomponents: the ten displayed leading operators E^{(0)} equal the u^0 part of the linearized EDFE on…" | V° |
| SM 41 | — | projected linearized EDFE components E_{pq̄} = V^M_p V̄^N_q̄ [ε](P_ε S_ε P̄_ε)_{MN} and E₀ = [ε](S₍₀₎,ε + 4/l²) | definition implemented by `LinearizedEDFEComponents` (NRH01) and used throughout NRH05 | D |
| SM 42 | `SMexactfieldorder` | Five-field and five-evolution-equation order | Definition. | D |
| SM 43 | `SMexactnormalform` | Radial normal form and its principal matrix | NRH05 "SMexactnormalform: M = diag(-1/4,-1/4,-1/4,-1/4,1) multiplies d_y^2 (h_pp,h_pm,h_mp,h_mm,delta d) in (E_pp,…" | V |
| SM 44 | `SMexactcoefficientpolynomial` | Polynomial extraction of exact operator coefficients | — | D |
| SM 45 | `SMexactTaylor` | Taylor coefficients at a regular finite slice | Definition. | D |
| SM 46 | `SMexactTaylorSolution` | Constrained local radial Taylor recursion | NRH05 "SMexactTaylorSolution: the displayed recursion makes the t^0 and t^1 coefficients of the normal form vanish…" | V |
| SM 47 | `SMexactCauchyconstraints` | Five radial constraints | NRH05 "SMexactCauchyconstraints: E_{p ybar}, E_{y qbar} and E_0 - 4 E_{y ybar} contain no second radial derivative…" | V |
| SM 48 | `SMexactconstraintpropagation` | Linearized Bianchi identity and constraint propagation | NRH05 "SMexactconstraintpropagation: delta G_{y~ N}: tangential = sqrt2(Vbar_N^a E_{y a} + V_N^a E_{a y}), N=y giv…" | V |
| SM 49 | `SMoperatorseries` | Operator expansion through exp(-2y/l) | Definition of the asymptotic expansion. | D |
| SM 50 | `SMradialansatz` | Independent sources and radial logarithmic coefficients | NRH05 "SMNRradialintegration: (d_y^2 + 2 d_y/l) annihilates 1 and u, gives 2/l on y and -2u/l on u y", "u^1 order on R: with the type and stress constraints imposed every remaining order-u equation is satisfied …" (+10 more) | V |
| SM 51 | `SMNRradialintegration` | Radial derivatives of u times a polynomial | NRH05 "SMNRradialintegration: (d_y^2 + 2 d_y/l) annihilates 1 and u, gives 2/l on y and -2u/l on u y" | V |
| SM 52 | `SMtotalorderhierarchy` | Total-order operator/field hierarchy | NRH05 "SMtotalorderhierarchy: E_I[u^n f] = u^n E_I(d_y - 2n/l, d_+, d_-) f for the ten E^{(0)}_I and all five fiel…" | V |
| SM 53 | `SMcoupledorders` | Combined equations at orders u0 and u1 | NRH05 "u^1 order on R: with the type and stress constraints imposed every remaining order-u equation is satisfied …", "u^1 order on NR: with the type and stress constraints imposed every remaining order-u equation is satisfied…" (+9 more) | V |
| SM 54 | `SMmixingexample` | R type-changing coefficient example | NRH05 "E^{(1),NR}_{om bop} = 0 while Delta E^{(1)}_{om bop} = -2 d_y delta d/l (the branch dependence of SMlogexam…" | V |
| SM 55 | `SMEzerocomponents` | the nine tensor components and the scalar of the common leading operator E^(0) on the radial-gauge fluctuation | NRH05 "SMEzerocomponents: the ten displayed leading operators E^{(0)} equal the u^0 part of the linearized EDFE on…" | V |
| SM 56 | `SMEoneNRcomponents` | the components of the first-order NR operator E^(1),NR for arbitrary L±(x^±) and W₁(x^+,x^−) | NRH05 "SMEoneNRcomponents: the ten displayed E^{(1),NR} equal the u^1 part of the linearized EDFE on the general N…" | V |
| SM 57 | `SMRoperatorrule` | Relation between the displayed R and NR first-order operators | NRH05 "SMRoperatorrule, SMDeltaEonecomponents: E^{(1),R} = E^{(1),NR}\|_{W1 -> 4 L+ L-} + Delta E^{(1)} with the te…" | V |
| SM 58 | `SMDeltaEonecomponents` | ΔE^(1) = E^(1),R − E^(1),NR|_{W₁→4L₊L₋}: the Riemannian-only first-order terms | NRH05 "SMRoperatorrule, SMDeltaEonecomponents: E^{(1),R} = E^{(1),NR}\|_{W1 -> 4 L+ L-} + Delta E^{(1)} with the te…" | V |
| SM 59 | `SMleadingdilaton` | Leading dilaton slope and type-changing logarithms | NRH05 "SMleadingdilaton, SMleadingstress, SMNRsol: the solved u^0 coefficients equal the displayed ones" | V |
| SM 60 | `SMleadingstress` | Leading diagonal logarithmic slopes | NRH05 "SMleadingdilaton, SMleadingstress, SMNRsol: the solved u^0 coefficients equal the displayed ones" | V |
| SM 61 | `SMNRsol` | Leading hair-channel polynomial coefficients | NRH05 "SMleadingdilaton, SMleadingstress, SMNRsol: the solved u^0 coefficients equal the displayed ones" | V |
| SM 62 | `SMlogexample` | General-background logarithmic coefficient or constraint | NRH05 "SMlogfreeconditions: dd^{(1)} = h^{(1)} = h^{(1b)} = 0 <=> d_+ d_- r = d_+^2 r = d_-^2 r = 0, delta d^{(0)}…", "E^{(1),NR}_{om bop} = 0 while Delta E^{(1)}_{om bop} = -2 d_y delta d/l (the branch dependence of SMlogexam…" | V |
| SM 63 | `SMRh2L_omopb` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 64 | `SMRdd2L` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 65 | `SMRh2L_opopb` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 66 | `SMRh2L_omomb` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 67 | `SMRdd2` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 68 | `SMRh2L_opomb` | General-background logarithmic coefficient or constraint | NRH05 "SMRh2L_*, SMRdd2L, SMRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 69 | `SMNRh2L_omopb` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 70 | `SMNRdd2L` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 71 | `SMNRdd2` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 72 | `SMNRh2L_opopb` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 73 | `SMNRh2L_omomb` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 74 | `SMNRh2L_opomb` | General-background logarithmic coefficient or constraint | NRH05 "SMNRh2L_*, SMNRdd2L, SMNRdd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})" | V |
| SM 75 | `SMtypeconstraint` | Integrated type-changing constraints with source zero modes | NRH05 "SMtypeconstraint on R: they are solved exactly by h^{(2)}_{om bop} = c_s + (…)", "SMtypeconstraint on NR: they are solved exactly by h^{(2)}_{om bop} = c_s + (…)" (+3 more) | V |
| SM 76 | `SMRconstraint0` | General-background logarithmic coefficient or constraint | NRH05 "SMRconstraint0, SMRconstraint1: d_- h^{(2)}_{op bop} = F_+^R and d_+ h^{(2)}_{om bom} = F_-^R with the disp…" | V |
| SM 77 | `SMRconstraint1` | General-background logarithmic coefficient or constraint | NRH05 "SMRconstraint0, SMRconstraint1: d_- h^{(2)}_{op bop} = F_+^R and d_+ h^{(2)}_{om bom} = F_-^R with the disp…" | V |
| SM 78 | `SMNRconstraint0` | General-background logarithmic coefficient or constraint | NRH05 "SMNRconstraint0, SMNRconstraint1: d_- h^{(2)}_{op bop} = F_+^NR and d_+ h^{(2)}_{om bom} = F_-^NR with the …" | V |
| SM 79 | `SMNRconstraint1` | General-background logarithmic coefficient or constraint | NRH05 "SMNRconstraint0, SMNRconstraint1: d_- h^{(2)}_{op bop} = F_+^NR and d_+ h^{(2)}_{om bom} = F_-^NR with the …" | V |
| SM 80 | `SMgeneralradialsolution` | Remaining integration functions and inverse derivatives | NRH05 "u^1 order on R: with the type and stress constraints imposed every remaining order-u equation is satisfied …", "u^1 order on NR: with the type and stress constraints imposed every remaining order-u equation is satisfied…" (+2 more) | V |
| SM 81 | `SMlogfreeconditions` | Restricted fixed-dilaton ensemble on the common vacuum | NRH05 "SMlogfreeconditions: dd^{(1)} = h^{(1)} = h^{(1b)} = 0 <=> d_+ d_- r = d_+^2 r = d_-^2 r = 0, delta d^{(0)}…" | V |
| SM 82 | `SMfullscalarvariation` | general variation of the scalar curvature: δS₍₀₎ = 2□δd − ∇_MΩ^M + 2(PSP̄)^{MN}V_M^p V̄_N^q̄ h_{pq̄} | — | N |
| SM 83 | `SMfullcurvaturevariation` | general variation of the tensor curvature δ(PSP̄)_{MN} (□ term, ∇Ω term and the two curvature-times-h terms) | — | N |
| SM 84 | `SMfullOmega` | Ω_M = e^{2d}∇_N δ(e^{−2d}H^N_M) = (V^{Np}V̄_M^q̄ + V_M^p V̄^{Nq̄})D_N h_{pq̄} − 2H_M^N ∂_N δd | — | N |
| SM 85 | `GammaDFT` | L_Γ² = e^{−2d}S₍₀₎ − ∂(e^{−2d}B) = e^{−2d}(PP − P̄P̄)(ΓΓ…) | NRH03 "SMgamma2 on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M) for arbitrary chiral L_pm", "SMgamma2 on R: the section fluxes B^{x pm} are built from the W-independent blocks and are total x-derivati…" (+2 more) | V |
| SM 86 | `variation` | Variation of the Gamma-squared Lagrangian | NRH06 "variation, SMBdefinition, SMAdefinition: delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) +…" | V |
| SM 87a | `SMBdefinition` | B^K := 2(P^{KM}P^{LN} − P̄^{KM}P̄^{LN})Γ_{LMN} = 4H^{KL}∂_L d − ∂_L H^{KL} | NRH06 "variation, SMBdefinition, SMAdefinition: delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) +…", "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" (+1 more) | V |
| SM 87b | `SMAdefinition` | A^K_{pq̄} := V^K_p Γ^L_{Lq̄} + V̄^K_q̄ Γ^L_{Lp} − Γ_p^K_q̄ − Γ_q̄^K_p | NRH06 "variation, SMBdefinition, SMAdefinition: delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) +…", "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" (+3 more) | V |
| SM 88 | `SMrenvariation` | On-shell variation of the action and leading counterterm at the cutoff | NRH06 "SMrenvariation, SMvolumecounterterm: on shell the bulk term cancels delta(-2 Lambda e^{-2d}) and the counte…" | V |
| SM 89 | `SMvolumecounterterm` | volume counterterm S_ct = −(1/16πG)(4/l)∫_{y=Y} d²x e^{−2d} | NRH06 "SMrenvariation, SMvolumecounterterm: on shell the bulk term cancels delta(-2 Lambda e^{-2d}) and the counte…" | V |
| SM 90a | `SMmetricradialdictionary` | dictionary ⟨K_{ab̄}⟩_s = −(1/32πG) FP_{Y→∞} e^{2Y/l} A^y_{ab̄} | NRH02 "Mainonepoints on R: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) (exact fam…", "Mainonepoints on R: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0" (+2 more); NRH06 "SMmetricradialdictionary, SMscalarradialdictionary: coefficient matching -2K = (16 pi G)^{-1} A^y and 2 T_(…", "SMonept/SMNRonept on R: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0" (+1 more) | V |
| SM 90b | `SMscalarradialdictionary` | dictionary ⟨T₍₀₎⟩_s = (1/16πG) FP_{Y→∞} e^{2Y/l}(B^y + 4/l) | NRH02 "Mainonepoints on R: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) (exact fam…", "Mainonepoints on R: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0" (+2 more); NRH06 "SMmetricradialdictionary, SMscalarradialdictionary: coefficient matching -2K = (16 pi G)^{-1} A^y and 2 T_(…", "SMonept/SMNRonept on R: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0" (+1 more) | V |
| SM 91 | `SMsecondvariation` | Independent tangent variations, frame and measure terms | NRH06 "SMsecondvariation: mixed coset constraint at second order", "SMsecondvariation: h_1 projection and delta_2 h_1^{p qbar} = 0" | V |
| SM 92 | `SMmomentumvariation` | δA^K_{pq̄}, δB^K including the frame-variation terms; 𝒜^K_MN | NRH06 "SMframemomentumcheck on R: the two frame terms give -(1/l) h_{a bbar} + O(u^2)", "SMframemomentumcheck on NR: the two frame terms give -(1/l) h_{a bbar} + O(u^2)" (+4 more) | V° |
| SM 93 | `SMBmomentumvariation` | δB^K = 8V^{(Kp}V̄^{L)q̄}h_{pq̄}∂_L d + 4H^{KL}∂_L δd − 2∂_L(V^{(Kp}V̄^{L)q̄}h_{pq̄}) | NRH06 "SMlinearizedmomenta on R: delta A^y_{a bbar} = (1/2) d_y h_{a bbar} + O(u^2), delta B^y = 4 d_y delta d + O…", "SMlinearizedmomenta on NR: delta A^y_{a bbar} = (1/2) d_y h_{a bbar} + O(u^2), delta B^y = 4 d_y delta d + …" | V° |
| SM 94 | `SMcurvedmomentumdefinition` | curved-index momentum 𝒜^K_{MN} := δ^K_M Γ^L_{LN} + δ^K_N Γ^L_{LM} − Γ_M^K_N − Γ_N^K_M | — | D |
| SM 95 | `SMprojectedconnectionmomentum` | V^M_p V̄^N_q̄ δ𝒜^K_{MN} = ½H^{KR}D_R h_{pq̄} − V^{Kr}D_p h_{rq̄} + V̄^{Kr̄}D_q̄ h_{pr̄} − 2V^K_p ∂_q̄ δd − 2V̄^K_q̄ ∂_p δd | NRH06 "SMconnectionmomentumcheck on R: V^M_a Vbar^N_bbar delta A^y_MN = (1/2) d_y h_{a bbar} + h_{a bbar}/l + O(u^2)", "SMconnectionmomentumcheck on NR: V^M_a Vbar^N_bbar delta A^y_MN = (1/2) d_y h_{a bbar} + h_{a bbar}/l + O(u^2)" | V° |
| SM 96 | `SMbackgroundmomenta` | Leading R/NR mixed momenta and scalar momentum | NRH06 "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)", "SMbackgroundmomenta on NR: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" | V |
| SM 97 | `SMsamechiralitymomenta` | Same-chirality tangential momenta | NRH06 "SMsamechiralitymomenta on R: tangential A^y_{ab} = eta/l, A^y_{abar bbar} = -etabar/l + O(u^2)", "SMsamechiralitymomenta on NR: tangential A^y_{ab} = eta/l, A^y_{abar bbar} = -etabar/l + O(u^2)" | V |
| SM 98 | `SMonept` | Riemannian one-point matrix | NRH06 "SMonept/SMNRonept on R: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0", "SMonept/SMNRonept on NR: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0" | V |
| SM 99 | `SMNRonept` | Non-Riemannian one-point matrix | NRH06 "SMonept/SMNRonept on R: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0", "SMonept/SMNRonept on NR: <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0" | V |
| SM 100 | `SMstatefluctuations` | Source-free tangent variations of the two saddle families | NRH06 "SMstatefluctuations on R: h_{a bbar} = u {{2 dL+, 2 d(L+L-)},{0, 2 dL-}} + O(u^2), delta d = O(u^2)", "SMstatefluctuations on NR: h_{a bbar} = u {{2 dL+, dW1/2},{0, 2 dL-}} + O(u^2), delta d = O(u^2)" | V |
| SM 101 | `SMderivativectNR` | quadratic derivative counterterm S^(2)_{ct,NR} on cutoff fields (e^{2Y/l}, Y and Y² terms) | NRH06 "SMderivativect on NR: delta_2 delta_1 S_ct^(2) cancels every e^{2Y/l} poly(Y) and Y^k divergence of the ord…" | V |
| SM 102 | `SMderivativectR` | S^(2)_{ct,R} = S^(2)_{ct,NR} − (1/32πG)∫[l²Y δd ∂₊²∂₋²h_{⊖⊕̄} + (l³Y²/16)h_{⊖⊕̄}∂₊³∂₋³h_{⊖⊕̄}] | NRH06 "SMderivativect on R: delta_2 delta_1 S_ct^(2) cancels every e^{2Y/l} poly(Y) and Y^k divergence of the orde…" | V |
| SM 103 | `SMctlocality` | stress constraints in second-derivative form ∂₋²h^(2)_{⊕⊕̄} = ∂₋F₊^s, ∂₊²h^(2)_{⊖⊖̄} = ∂₊F₋^s (locality of the divergent response terms) | NRH06 "SMctlocality on R: with the stress constraints the divergent rows of the bare bilinear are local in both sl…", "SMctlocality on NR: with the stress constraints the divergent rows of the bare bilinear are local in both s…" | V |
| SM 104 | `SMctadjoint` | (L₊∂₊ + ½∂₊L₊)† = −(L₊∂₊ + ½∂₊L₊), [∂₋, L₊] = 0 | NRH06 "SMctadjoint: (L+ d+ + (1/2) d+ L+)^dagger = -(L+ d+ + (1/2) d+ L+) on test functions" | V |
| SM 105 | `SMdirectsources` | source abbreviations j = (a,b,r,c,v) and (ε_s, J_s) = (1, L₊L₋)_R, (0, W₁/4)_NR | definition; the slot structure of NRH06 | D |
| SM 106 | `SMdirectbackground` | background terms e^{−2d_s} = u^{−1}[1+O(u²)], A^y_{ab̄} = −(2u/l)((L₊, J_s),(ε_s, L₋)) + O(u²), B^y + 4/l = O(u²) | NRH06 "SMbackgroundmomenta on R: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)", "SMbackgroundmomenta on NR: A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)" (+2 more) | V° |
| SM 107 | `SMframemomentumcheck` | Contribution of both moving-frame projections | NRH06 "SMframemomentumcheck on R: the two frame terms give -(1/l) h_{a bbar} + O(u^2)", "SMframemomentumcheck on NR: the two frame terms give -(1/l) h_{a bbar} + O(u^2)" | V |
| SM 108 | `SMconnectionmomentumcheck` | Projected connection contribution through u | NRH06 "SMconnectionmomentumcheck on R: V^M_a Vbar^N_bbar delta A^y_MN = (1/2) d_y h_{a bbar} + h_{a bbar}/l + O(u^2)", "SMconnectionmomentumcheck on NR: V^M_a Vbar^N_bbar delta A^y_MN = (1/2) d_y h_{a bbar} + h_{a bbar}/l + O(u^2)" | V |
| SM 109 | `SMlinearizedmomenta` | Off-shell tangential momentum variations through u | NRH06 "SMlinearizedmomenta on R: delta A^y_{a bbar} = (1/2) d_y h_{a bbar} + O(u^2), delta B^y = 4 d_y delta d + O…", "SMlinearizedmomenta on NR: delta A^y_{a bbar} = (1/2) d_y h_{a bbar} + O(u^2), delta B^y = 4 d_y delta d + …" | V |
| SM 110 | `SMdirectcutoffbilinear` | ordered cutoff bilinear δ₂δ₁S_ren|_Y (u^{−1}[½h₁^{ab̄}∂_y h₂ + 8δ₁d∂_yδ₂d] − (4/l)δ₂d(L₊a₁ + L₋b₁ + J_s r₁ + ε_s c₁)) + δ₂δ₁S^(2)_ct | NRH06 "SMdirectcutoffbilinear on R: the general second variation with the computed momenta equals the displayed or…", "SMdirectcutoffbilinear on NR: the general second variation with the computed momenta equals the displayed o…" | V |
| SM 111 | `SMfinitepartexplicit` | Finite coefficients after radial differentiation | NRH06 "SMfinitepartexplicit on R: u^{-1} delta A^y_{a bbar} = u^{-1}(h^(1)/2 + y h^(1b)) - h^(2)/l + h^(2L)/2 - (y…", "SMfinitepartexplicit on NR: u^{-1} delta A^y_{a bbar} = u^{-1}(h^(1)/2 + y h^(1b)) - h^(2)/l + h^(2L)/2 - (…" | V |
| SM 112 | `SMfullfinitehessian` | finite part FP δ₂δ₁S_ren in terms of the expansion coefficients h^(0), h^(1), h^(2), h^(2L), δd^(k) | NRH06 "SMfullfinitehessian on R: the Y^0 coefficient of the ordered bilinear equals the displayed finite integrand", "SMfullfinitehessian on NR: the Y^0 coefficient of the ordered bilinear equals the displayed finite integrand" | V |
| SM 113 | `SMcontactconstraintreduction` | integration by parts of the slot-1 response terms using the stress constraints: (l/4)∫[h^(2)_{1,⊖⊖̄}∂₊²r₂ + h^(2)_{1,⊕⊕̄}∂₋²r₂] = (l/4)∫r₂[∂₊F₋ + ∂₋F₊] | NRH06 "SMcontactconstraintreduction on R: after integrating the slot-1 response terms by parts the finite part is …", "SMcontactconstraintreduction on NR: after integrating the slot-1 response terms by parts the finite part is…" | V |
| SM 114 | `SMdirectstressrows` | rows q_{s,1}, q_{s,2} (stress channels) of the finite bilinear | NRH06 "SMdirectstressrows, SMdirectremainingrows on R: the Euler-Lagrange derivatives of the finite bilinear with …", "SMdirectstressrows, SMdirectremainingrows on NR: the Euler-Lagrange derivatives of the finite bilinear with…" | V |
| SM 115 | `SMdirectremainingrows` | rows q_{s,3}, q_{s,4}, q_{s,5} (type-changing, zero-mode and dilaton channels) | NRH06 "SMdirectstressrows, SMdirectremainingrows on R: the Euler-Lagrange derivatives of the finite bilinear with …", "SMdirectstressrows, SMdirectremainingrows on NR: the Euler-Lagrange derivatives of the finite bilinear with…" | V |
| SM 116 | `SMdirectkernel` | δ₂δ₁S_ren = (1/16πG)∫Σ_I j_{1,I} q_{s,I}[j₂] + δ₂δ₁S_fin, kernels Q_{s,IJ}, (G_s)_{IJ}|action response = Q_{s,IJ}/(64πG) + ¼δ²S_fin/δjδj | NRH06 "SMdirectkernel, SMRscalarcontactrow: (G_s)_{5J}\|action response = q_5/(4 * 16 pi G) reproduces the displaye…", "SMdirectkernel, SMNRscalarcontactrow: (G_s)_{5J}\|action response = q_5/(4 * 16 pi G) reproduces the display…" | V |
| SM 117 | `SMRscalarcontactrow` | Riemannian scalar contact row (G_R)_{5J} = (l/64πG)(−∂₊², −∂₋², (l²/4)∂₊²∂₋² − L₊∂₋² − L₋∂₊², 0, 4∂₊∂₋)δ² | NRH06 "SMdirectkernel, SMRscalarcontactrow: (G_s)_{5J}\|action response = q_5/(4 * 16 pi G) reproduces the displaye…" | V |
| SM 118 | `SMNRscalarcontactrow` | non-Riemannian scalar contact row (G_NR)_{5J} = −(l/64πG)(0, 0, L₊∂₋² + L₋∂₊², 0, 0)δ² | NRH06 "SMdirectkernel, SMNRscalarcontactrow: (G_s)_{5J}\|action response = q_5/(4 * 16 pi G) reproduces the display…" | V |
| SM 119 | `SMLorentziandeltaconvention` | Lorentzian delta convention δ²(x−x′) = (1/2πi)∂₊∂₋ ln(−Δ₊Δ₋ + iε) | NRH06 "SMLorentziandeltaconvention: the box integral of (1/(2 pi i)) d+ d- ln(-x+ x- + i eps) over [-R,R]^2 is 1" | V |
| SM 120 | `SMfullhessianclosure` | symmetry (G_s)_{IJ}(x,x′) = (G_s)_{JI}(x′,x) of the action-response matrix | — | S |
| SM 121 | `SMstresssourceexample` | Independent stress-source constraints | NRH06 "SMstresssourceexample: with only h^(0)_{om omb} the R stress constraint is d_- h^(2)_{op opb} = -(2L+ d+ + …", "SMstresssourceexample: on NR the third-derivative term is absent" | V |
| SM 122 | `SMgeneralWardoperators` | Boundary differential operators for stress and hair | NRH06 "SMstresssourceexample: with only h^(0)_{om omb} the R stress constraint is d_- h^(2)_{op opb} = -(2L+ d+ + …", "SMstresssourceexample: on NR the third-derivative term is absent" (+1 more) | V° |
| SM 123 | `SMstressparticular` | Particular stress kernels from an inverse derivative | NRH06 "SMstressparticular: -(1/l)/(64 pi G)/i = -kappa/(4 i) with kappa = 1/(16 pi G l)" | V |
| SM 124 | `SMindependentsourcegenerator` | Radial compensator keeping all other leading sources zero | NRH06 "SMindependentsourcegenerator: radial gauge h_{p ybar} = h_{y qbar} = 0 through u on both saddles, and delta…", "SMindependentsourcegenerator: the only leading source is h^(0)_{om omb} = -2 d_- alpha^+ (both saddles)" (+1 more) | V |
| SM 125 | `SMhairparticularsolution` | Full nonchiral hair response including its local term | NRH06 "SMhairparticularsolution: h^(2)_{op opb} = 2 V_+^s alpha^+ on both saddles", "SMhairparticularsolution: h^(2)_{op omb}\|R = 2 L- V_+^R alpha - (l^2/2) L+ d+ d-^2 alpha, \|NR = (1/2) W_+ a…" | V |
| SM 126 | `SMgeneralparticularkernels` | R and NR particular mixed kernels | NRH06 "SMgeneralparticularkernels: M_+^R = L_- A_+^R and M_+^NR = -(kappa/(16 i)) W_+ d_-^{-1} delta (normalizatio…" | V |
| SM 127 | `SMPBHkernel` | Local-plane Green-function convention | NRH06 "SMPBHkernel: d_-^{-1} delta^2 = (1/(2 pi i)) d_+ ln(-D+ D- + i eps) -> 1/(2 pi i D+) away from the light cone" | V |
| SM 128 | `SMgeneralpositionkernels` | Position-space particular A and M poles | NRH06 "SMgeneralpositionkernels, Rcorrelators: A_+^NR = [L+'/D+ - 2L+/D+^2]/(128 pi^2 G l)", "SMgeneralpositionkernels, Rcorrelators: A_+^R = A_+^NR + 3l/(256 pi^2 G D+^4)" (+1 more) | V |
| SM 129 | `SMRtwopt` | Riemannian plane-vacuum stress normalization | NRH06 "SMRtwopt: at L+ = 0, (8 pi)^2 A_+^R = (c/2)/D+^4 with c = 3l/(2G)" | V |
| SM 130 | `SMNRhairhessian` | Hair self-response as the derivative of an undetermined function | NRH06 "SMNRhairhessian: the hair-hair coefficient (1/l)/(64 pi G)/i = 1/(64 pi i G l)" | V |
| SM 131 | `SMRhairuniformization` | uniformization {f±; x±} = −4L±/l² (Schwarzian) and the thermal factor B±(x,x′) = f′(x)²f′(x′)²/[f(x) − f(x′)]⁴ | NRH06 "SMRhairuniformization: {e^{2 alpha x}; x} = -2 alpha^2 = -4 L/l^2 for alpha = Sqrt[2L]/l" | V |
| SM 132 | `SMRhairconditional` | connected R type-changing two-point function: (9l⁵/1024π²G)B₊B₋ + (3l/256π²G)[L₋L₋′B₊ + L₊L₊′B₋] | NRH06 "SMRhairconditional: the stress coefficient 3l/(256 pi^2 G) equals (8 pi)^{-2} c/2", "SMRhairconditional: 9 l^5/(64 pi^2 G r^8) with r^2 = -2 D+ D- is 9 l^5/(1024 pi^2 G) (-D+ D-)^{-4}" | V |
| SM 133 | `SMRhairconstantL` | constant-L thermal factor B± = [α±/sinh(α±Δ±)]⁴ with α± = √(2L±)/l | NRH06 "SMRhairconstantL: B(x,x') = [alpha/sinh(alpha D)]^4 for f = e^{2 alpha x}", "SMRhairconstantL: the L -> 0 limit of the thermal factor is 1/D^4" | V |
| SM 134 | `SMonepointequivalence` | equivalence (G_s)_{IJ}|action response = ½δÔ_I/δj_J = ¼δ²S_ren/δj_Iδj_J at j = 0 | NRH06 "SMonepointequivalence: (1/2) delta O^_I/delta j_J with O^_I = (1/2) delta S/delta j_I is (1/4) delta^2 S/de…" | V |
| SM 135 | `SMCPSform` | surface-charge one-form | NRH07 "lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary" | V |
| SM 136 | `SMRpotentialcomponents` | lim e^{−2d}K̂ components (the minus-sector check uses K̂^{y+} = −K̂^{+y}) | NRH07 "lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''" | V |
| SM 137 | `SMCPSresult` | Θ̂ → 0, k = (4/l)εδL | NRH07 "lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary", "k^{-y}[eps+] = (4/l) eps+ dL+" | V |
| SM 138 | `SMNRchargefalloffs` | state-dependent falloffs | NRH07 "state-dependent falloffs delta H^-_+ = 2 z dL+, delta H^+_- = -2 z dL-, delta H_{+-} = z dW1" | V |
| SM 139 | `SMNRchargecancellation` | componentwise W₁ cancellation | NRH07 "W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise" | V |
| SM 140 | `SMNRcharge` | δQ, Q | NRH07 "(i) k^{-y}[eps+] = delta[(4/l) eps+ L+] (the charge exists and is integrable)" | V |
| SM 141 | `SMCPSalgebra` | {Q,Q} = Q[[ε,η]], c_charge = 0 (+ footnote) | NRH07 "the same-chirality C-bracket closes up to a closed B-gauge parameter (slot x~+ only)", "(ii) {Q[e1+], Q[e2+]} - Q[[e1,e2]] is a total derivative: centerless plus sector" (+4 more) | V |
| SM 142 | `SMdyg` | doubled-yet-gauged action | — (defining action; its reductions are verified below) | D |
| SM 143 | `SMphysicalsectionA` | Â_αμ, D_αx^M on the section | — | D |
| SM 144 | `SMRfirstorder` | Riemannian first-order form | NRH08 "the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+", "eliminating the auxiliaries reproduces E_{mu nu} dx^mu dbar x^nu, E = g - B" | V |
| SM 145 | `SMceff` | c_eff² = 2F | NRH02 "c_eff^2 := 2F = 2(e^{2y/l} + L+ L- e^{-2y/l}) ~ 2 e^{2y/l}: the Gomis-Ooguri limit is reached without tuning"; NRH08 "c_eff^2 = 2F -> 4 Sqrt[L+L-] at the horizon u = Sqrt[L+L-]" | V |
| SM 146 | — | string momentum density P_μ and the long-string energy E = −Q_{∂t} = (1/2πα′)∫dσ(√−γ − wB_{tφ}) | NRH08 "before : gamma^{tau a} g_{t nu} d_a X^nu = gamma^{tau a} gamma_{a tau} = 1 on the static embedding", "before : B_{t phi} = l (e^{2y/l} + L+L- e^{-2y/l}) and B_{t nu} X'^nu = w B_{t phi}" (+1 more) | V |
| SM 147 | `SMlongstringE` | winding-string energy of the static probe (SM 5.1) | NRH08 "before : E = -Int_0^{2 pi} d sigma P_t with the displayed P_mu reproduces E(y) = -(2 w l/alpha') L+L- e^{-2…", "E(y) = (w l/alpha')[(e^{2y/l} - L+L- e^{-2y/l}) - (e^{2y/l} + L+L- e^{-2y/l})] = -(2 w l/alpha') L+L- e^{-2…" | V |
| SM 148 | `SNCtau` | SNC clock forms | NRH08 "SM: the dual vectors Y, Ybar exist at every radius (unit clock determinant)"; NRH09 "x rows = Sqrt[2] tau^pm, x~ rows = the dual vectors -Y/Sqrt[2], Ybar/Sqrt[2], and the W entries = -(W/(2 Sq…" | V |
| SM 149 | `SMdygconstraints` | τ⁺·∂̄x = 0 = τ⁻·∂x | NRH08 "SM: H^{mu nu} tau^pm_nu = 0 (two-dimensional longitudinal kernel)", "SM: the dual vectors Y, Ybar exist at every radius (unit clock determinant)" | V |
| SM 150 | `SMdygGO` | exact reduced Lagrangian | NRH08 "at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)" | V |
| SM 151 | `SMGO` | Gomis–Ooguri limit | NRH08 "at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)" | V |
| SM 152 | `SMvertex` | V_W = (1/4πα′)W∂x⁺∂̄x⁻ | NRH08 "with the 1/(2 pi alpha') prefactor this is V_W = (1/(4 pi alpha')) W dx+ dbar x-" | V |
| SM 153 | `SMWisB` | H(W) = Ω_bH(0)Ω_bᵀ, b_{+−} = −W/2 | NRH04 "SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched" | V |
| SM 154 | `SMdeltaL` | antisymmetric clock coupling | NRH08 "symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)" | V |
| SM 155 | `SMFTsector` | S_y, S_FT | — (definitions; used in SM 158) | D |
| SM 156 | `SMWZWlevel` | k = l²/α′ | NRH08 "k = \|Int_{S^3} H\| / (4 pi^2 alpha') = l^2/alpha'" | V |
| SM 157 | `SMFTweight` | T_y, h_y(a) | NRH08 "the two OPE contributions assemble to h_y(a)/(z-w)^2" | V |
| SM 158 | `SMBRSTradial` | (∂²_y + (2/l)∂_y)f = 0 | NRH08 "the marginal roots of h_y are a = 0 and a = -2/l (modes {1, e^{-2y/l}})", "(d_y^2 + (2/l) d_y) f = 0 for f = W0 + W1 e^{-2y/l}" | V |
| SM 159 | `SMBRSTmomentum` | h_y(p_y), roots p_y = 0, 2i/l | NRH08 "h_y(i p_y) = (alpha'/4) p_y (p_y - 2 i/l) and P_y = p_y - i/l gives (alpha'/4)(P_y^2 + 1/l^2)" | V |
| SM 160 | `SMBRSTvertex` | U_{W₁} is a weight-(1,1) primary | NRH08 "V_W x V_W is nonsingular: no contraction pair lies inside {x^+, x^-}^2", "<x^+ x^-> = 0 in the Gomis-Ooguri system (x^+ pairs only with beta)" | V° |
| SM 161 | `SMBRSTcentral` | c_y = 1 + 6α′/l², total 3(k+2)/k | NRH08 "c_{beta gamma} + c_y = 2 + (1 + 6/k) = 3(k+2)/k" | V |
| SM 162 | `SMWgaugeobstruction` | W₀ gauge, e^{−2y/l}W₁ not | NRH08 "Lhat_xi H^infty = D H^infty + H^infty D^T with D = ((-(dv)^T, 0), (b, dv)), b = d lambda~", "the conditions force d_y v^mu = 0, d_- v^+ = 0 = d_+ v^-, b_{+-} = -varpi/2, b_{+y} = d_+ v^y, b_{-y} = -d_…" | V |
| SM 163 | `SMBRSTfusion` | self-contraction, h_n, resonances | NRH08 "e^{a y(z)} e^{a y(0)} ~ \|z\|^{-alpha' a^2} = \|z\|^{-4 alpha'/l^2} for a = -2/l (from <y y> = -(alpha'/2) Log\|…", "the n-fold fused weight h_n = n + q n(1-n) equals n + h_y(-2n/l) with q = alpha'/l^2" (+1 more) | V |
| SM 164 | `SMupliftblocks` | H₁₀ = H₃ ⊕ H_{S³} ⊕ H_{R⁴}, d₁₀ = d₃ + d_{S³} | NRH09 "…: S_(0)^{(10)} = 0", "…: (P S Pbar)^{(10)} = 0" | V |
| SM 165 | `SMDFTKilling` | generalized Killing equations | NRH09 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm", "Lhat_xi d = 0 (the radial component compensates the divergence)" | V |
| SM 166 | `SMtypeIIKS` | type-II Killing-spinor systems | — (definition of the systems solved in NRH07) | D |
| SM 167 | `SMsusyclosure` | closure on ĥL_X, X = iε̄₂Γε₁ | NRH09 "C = i sigma_2 is the Majorana conjugation for the representation: C gamma^p C^-1 = -(gamma^p)^T", "Lhat_X H^infty = 0 and Lhat_X d = 0 for the Killing-spinor bilinear, arbitrary chiral f1, f2" (+1 more) | V |
| SM 168 | `SMsemicov` | semi-covariant connection | NRH05 "SMbackgroundconnection on R: nabla_C P_AB = 0 through z^2", "SMbackgroundconnection on NR: nabla_C P_AB = 0 through z^2" (+4 more); NRH09 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2 (pairwise cancellation with AdS3)", "S_(0)(S3 block) = +4/l^2 via the closed form (cancels -4/l^2 of either 3d saddle)" | V° |
| SM 169 | `SMuplift` | AdS₃×S³×R⁴ with flux | NRH09 "R(S3) = +6/l^2 and H^2(S3) = +24/l^2 (pairwise cancellation with AdS3)", "S_(0)(S3 block) = +4/l^2 via the closed form (cancels -4/l^2 of either 3d saddle)" (+1 more) | V |
| SM 170 | `SMRDFTKilling` | Ordinary-field form of the generalized Killing equations | Cited decomposition; the ordinary Killing vectors are not enumerated in the suite. | S |
| SM 171 | `SMRlocaliso` | k^{(ij)} = s_is_j stabilizers | NRH09 "k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' - (l^2/4) k''' = 0 (the stabilizer equation)" | V |
| SM 172 | — | vol₃, H₃ = −(2/l)vol₃, H_{y−+} = +(2/l)√|g₃| (flux orientation) | NRH09 "flux orientation: H_{y-+} = +(2/l) Sqrt[\|g_3\|] = (2/l)(e^{2y/l} - L+L- e^{-2y/l})" | V |
| SM 173 | — | AdS₃ Killing-spinor equation ∇_μ ε± = ±(1/2l)γ_μ ε± | — | S |
| SM 174 | `SMRcomponentHill` | first-order pair → Hill equation | NRH09 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s" | V |
| SM 175 | `SMRlocalKS` | Hill equations for s± | NRH09 "the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s", "massless BTZ (L = 0): the (u, v) pair shifts by exactly 2 pi v over one circuit (unipotent)" (+1 more) | V |
| SM 176 | — | N_local = 16 real = 8_{Spin(1,9)} + 8_{Spin(9,1)} | — | S |
| SM 177 | `SMvielbein` | exact non-Riemannian double vielbein | NRH05 "SMexactNRdata: the lowered SMvielbein frames reconstruct P, Pbar of the exact hairy metric and are orthogonal" | V |
| SM 178 | `SMvielbeincheck` | defining relations | NRH09 "V_M^p V_Np = P_MN and Vbar_M^pbar Vbar_Npbar = Pbar_MN for the exact frame , any hh, sigma, W", "V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)" (+1 more) | V |
| SM 179 | `SMNRlocalstabilizer` | local stabilizer system | NRH09 "eps = c/Sqrt[L] solves eps dL + 2 L d eps = 0 (the generic obstruction)", "exact one-sided identity Lhat_xi H = delta H, with weights (1,2) for (W0, W1)" | V |
| SM 180 | `SMexactiso` | vacuum isometries | NRH09 "Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm", "X^M = (omega_+, 0, -(l/2) v'; v^+, 0, -(l/2) v') with v^+ = 2 f1 f2, omega_+ = -2 l^2 f1' f2'" | V |
| SM 181 | `SMweighteddilaton` | weighted dilaton condition | NRH09 "Lhat_xi d = 0 (the radial component compensates the divergence)", "equivalently d_M(e^{-2d} xi^M) = 0" | V |
| SM 182 | `SMkillingspinor` | vacuum Killing spinor E = (√2f, l∂₊f) | NRH09 "SM: vacuum spin connection = displayed Phi_{~+ oplus y} = 1/(2l), Phi_{+ ominus y} = -1/l (raised slots)", "D_{pbar} E = 0 for E = (Sqrt[2] f(x+), l f'(x+)), arbitrary chiral f" (+2 more) | V |
| SM 183 | `SMcomplexblocks` | σ_i, τ_i, ρ_m | NRH09 "rep: {gamma^p, gamma^q} = 2 eta^{pq} (3d lightcone blocks)" | V |
| SM 184 | — | index convention ε_α ≡ ε_{α₁α₂α₃α₄}, α₁,α₂,α₃ = 1,2, α₄ = 1,…,4 | — | D |
| SM 185 | `SMgammaten` | ten-dimensional Γ^p̂, Γ₁₁ | NRH09 "{Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32", "Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0" | V |
| SM 186 | `SMgammaMajorana` | Majorana intertwiner | NRH09 "BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, same for Gamma_11, and BB10 BB10^* = 1" | V |
| SM 187 | `SMgammabarred` | barred Clifford algebra | NRH09 "{Gammabar, Gammabar} = -2 eta, Gammabar_{pq} = -Gamma_{pq}, same Majorana intertwiner" | V |
| SM 188 | `SMreducedDirac` | reduced two-component system | NRH09 "SM text after SMreducedDirac: D_{pbar} E has the components 0, -L+ u e0/l, em0/Sqrt2, em1/Sqrt2 - W1 u e0/(…", "SM text after SMreducedDirac: gamma^p D_p E + E/(Sqrt2 l) has the components ey0/Sqrt2 + Sqrt2 e0/l and ep0…" | V |
| SM 189 | `SMinternalprojectors` | complex product spinors with auxiliary/R⁴ signs (−,−) and (+,+) | only the historical complex rank count is checked; the product basis and the real Majorana combinations are not separately verified; NRH09 "Weyl + S^3-line + zeta_+ leave complex dimension 4 (the Xi_{+r} span)" | V° |
| SM 190 | — | real Majorana combinations Ψ + B₁₀^{−1}Ψ*, i(Ψ − B₁₀^{−1}Ψ*) | NRH09 "the Majorana condition leaves 32 real components" | D |
| SM 191 | `SMspinorcountchain` | Majorana, Weyl and background-projector count | NRH09 "the Majorana condition leaves 32 real components", "adding the Weyl condition leaves 16 real components" (+1 more) | V° |
| SM 192 | `SMhairyKS` | four arbitrary real chiral functions on each one-sided NR branch (hairy Killing spinors) | the reduced one-sided jet system in the printed frame is verified; the full ten-dimensional fermionic equations with the current internal basis and nonzero supercharges are not; NRH09 "the displayed system has rank six (seven relations, one dependent) and spans the derived equations", "the solution space is exactly {e_1, d_+ e_1}: e_0 = 0, e_1 = F_+(x^+) with d_+ e_1 unconstrained" (+1 more) | V° |
| SM 193 | — | Fourier expansion of the chiral functions F^r_±(x^±) with reality f_{−n} = f_n* | — | D |
| SM 194 | — | formal pairing δQ[F] = ∫dφ Σ_r F^r δS_r = 2π Σ f_n^r δS_{r,−n} | — | D |
| SM 195 | `SMcandidateD` | gauge derivative 𝔻 | — (definition; used in NRH08) | D |
| SM 196 | `SMcandidateaction` | candidate action and its reduction | NRH10 "the P^(0) term contains only D_+ and the Pbar^(0) term only D_-", "each term selects a single independent component of phi_A" | V |
| SM 197 | `SMcandidategauge` | internal gauge invariance | NRH10 "delta_Lambda L = 0 for the non-abelian bosonic sector" | V |
| SM 198 | `SMcandidateWitt` | Witt transformations | NRH10 "delta_v L = d_+(v^+ L) + d_-(v^- L) (unit-weight scalar density)", "Lhat_xi H^(0) = 0 <=> d_- v^+ = 0 = d_+ v^-" | V |
| SM 199 | `SMcandidateWittL` | δ_vL = ∂(vL) | NRH10 "delta_v L = d_+(v^+ L) + d_-(v^- L) (unit-weight scalar density)" | V |
| SM 200 | `SMcandidatefermionic` | chiral fermionic transformations | NRH10 "delta_{eps+} L = d_-( psi^+ delta_{eps+} psi^+ )", "the fermionic parameters carry arbitrary chiral profiles (one function per chirality)" | V |
| SM 201 | `SMcandidatefermionicL` | their total-derivative variations | NRH10 "delta_{eps+} L = d_-( psi^+ delta_{eps+} psi^+ )", "delta_{eps-} L = d_+( psi^- delta_{eps-} psi^- )" | V |
| SM 202 | `SMcandidateextra` | Grassmann-even chiral transformation | NRH10 "delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]" | V |
| SM 203 | `SMcandidateextraL` | its variation | NRH10 "delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]" | V |
| SM 204 | `SMcandidatetrivial` | equation-of-motion redundancies | NRH10 "the transformations are built from the fermion equations of motion d_- psi^+ = 0 = d_+ psi^-, hence vanish …" | V |
| SM 205 | `SMcandidatetrivialL` | δ_αL total derivative | NRH10 "the two-component fermions used here are genuinely Grassmann: psi^+ d_- psi^+ is nonzero off shell" | V |
