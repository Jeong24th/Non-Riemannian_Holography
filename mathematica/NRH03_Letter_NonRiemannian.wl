(* ::Package:: *)

(* ::Title:: *)
(*NRH03 — Letter, "Everywhere Non-Riemannian Saddle and Long-String Hair"*)


(* ::Text:: *)
(*This file verifies the displayed equations of the Letter's non-Riemannian section and*)
(*the associated exact radial statements of SM 4:*)
(**)
(*  Eq. (12) [NRvariables]   the variables (Pi, q, e^sigma, chi);*)
(*  Eq. (13) [NRHcompact]    the everywhere non-Riemannian generalized metric;*)
(*  Eq. (14) [NRdilaton]     the DFT dilaton  e^{-2d} = e^{2y/l}(1-q^2);*)
(*  Eq. (15) [NRWgeneral]    the exact hair profile W(x^+, x^-, chi), and the claim that*)
(*                           (13)-(14) with (15) solve  G_MN = 2 l^-2 J_MN  for arbitrary*)
(*                           chiral L_pm(x^pm) and arbitrary W_0(x), W_1(x);*)
(*  SM (61)-(67) [NRhill-NRGprofile] the exact radial reduction d^2 W/d chi^2 = F;*)
(*  SM (70)-(71) [SMBtransform, SMWshift] the finite B-shift that removes W_0;*)
(*  SM (77) [SNCtau]         the SNC clock forms and their Banados-frame asymptotics;*)
(*  SM (82) [SMWisB]         the pointwise W-as-B identity used on the worldsheet;*)
(*  Eq. (16) [NRasympt]      the asymptotic transformation laws: delta_eps L_pm without*)
(*                           anomaly and the inhomogeneous third derivative in delta_eps W_1*)
(*                           (verified in the near-boundary expansion, where W_1 is defined);*)
(*  SM (108) [SMNRlocalstabilizer], second line: the weight-one law of W_0 (two-sided, leading order);*)
(*  Eq. (12) and SM (61) [NRhill]: the sinh(chi/2) expansion and the Hill identity A = psi''/psi.*)
(**)
(*Strategy.  The radial dependence of the exact saddle enters only through chi, so we use*)
(*worksheet coordinates (x^+, x^-, chi) with the exact chain rules*)
(*   d chi/dy    = -(2 Sqrt[2]/l) Sinh[chi/Sqrt[2]] ,*)
(*   d chi/dx^pm = -Sqrt[2] Sinh[chi/Sqrt[2]] psi_pm'/psi_pm ,     psi_pm := L_pm^{-1/2},*)
(*(equivalent to SM (61)-(63); the rules are themselves re-derived below) and treat the*)
(*explicit -y/l term of the dilaton through a bookkeeping symbol Y with dY/dy = 1.*)
(*For the exact zero tests every hyperbolic function is rationalized by  chi -> 2 Sqrt[2] Log[T],*)
(*so that expressions become Laurent polynomials in the algebraically independent monomials*)
(*T^1 and T^(2 Sqrt[2]); Mathematica then decides each identity exactly.*)


ClearAll["Global`*"];
Get[FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH03_Letter_NonRiemannian.wl"];

JJ = ODDJ[3];

(* chain rules *)
chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
dpOp = Function[e, D[e, xp] + chp D[e, ch]];
dmOp = Function[e, D[e, xm] + chm D[e, ch]];
dyOp = Function[e, chy D[e, ch] + D[e, Ysym]];
xsNR = {dpOp, dmOp, dyOp};

(* exact zero test: rationalize chi -> 2 Sqrt[2] Log[T] *)
NRHZeroNR[label_, e_] := NRH`CheckZero[label,
   Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];

(* the exact non-Riemannian generalized metric, Eq. (13), with generic W *)
esig = psim[xm]/psip[xp];   (* e^sigma = Sqrt[L_+/L_-] = psi_-/psi_+ *)
HNRof[Wexpr_] := Module[{CC = Cosh[ch], SS = Sinh[ch]},
   {{0, 0, 0, CC, -SS/esig, 0},
    {0, 0, 0, esig SS, -CC, 0},
    {0, 0, 1, 0, 0, 0},
    {CC, esig SS, 0, -Wexpr esig SS, Wexpr CC, 0},
    {-SS/esig, -CC, 0, Wexpr CC, -Wexpr SS/esig, 0},
    {0, 0, 0, 0, 0, 1}}];
dNR = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];


(* ::Section:: *)
(*Eqs. (12)-(14): constraint, non-Riemannian type, dilaton, chain rules*)


HW = HNRof[W[xp, xm, ch]];
NRHZeroNR["Eq.(13): H J H = J for the exact non-Riemannian matrix (generic W)",
   HW . JJ . HW - JJ];
NRH`Check["Eq.(13): type (1,1) at every radius: upper-left block = diag(0,0,1)",
   HW[[1 ;; 3, 1 ;; 3]] === {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}];
NRHZeroNR["Eq.(14): e^{-2d} e^{-2y/l} = 1 - q^2 with q = tanh(chi/(2 Sqrt[2]))",
   Exp[-2 (dNR + Ysym/l)] - (1 - Tanh[ch/(2 Sqrt[2])]^2)];
(* the radial chain rule re-derived:  d chi/dy = -(4 Sqrt[Pi]/l) e^{2d}.  With
   e^{2d} = e^{-2y/l} cosh^2(chi/(2 Sqrt[2])) and e^{-2y/l} = Sqrt[2] q psi_+ psi_-,
   q = tanh(chi/(2 Sqrt[2])), Sqrt[Pi] = 1/(psi_+ psi_-), the psi's cancel and the claim
   becomes the pure hyperbolic identity checked here. *)
NRHZeroNR["SM(61)-(63): d chi/dy = -(4 Sqrt[Pi]/l) e^{2d}  (hyperbolic identity form)",
   chy + (4 Sqrt[2]/l) Tanh[ch/(2 Sqrt[2])] Cosh[ch/(2 Sqrt[2])]^2];
NRHZeroNR["SM: RG rapidity  mu d chi/d mu = -2 Sqrt[2] Sinh[chi/Sqrt[2]]",
   l chy + 2 Sqrt[2] Sinh[ch/Sqrt[2]]];


(* ::Section:: *)
(*SM (61)-(67) [NRhill-NRGprofile]: exact radial operator identities*)


chiq = 2 Sqrt[2] ArcTanh[q];
NRH`CheckZero["SM(62) [NRradialchange]: d chi/d q = 2 Sqrt[2]/(1-q^2)",
   Together[D[chiq, q] - 2 Sqrt[2]/(1 - q^2)]];
NRH`CheckZero["SM(63) [NRradialoperator]: transformed radial operator equals f''(chi(q))",
   Together[(1 - q^2)^2/8 (D[ff[chiq], {q, 2}] - 2 q/(1 - q^2) D[ff[chiq], q])
      - Derivative[2][ff][chiq]]];
rho = Together[4 Sqrt[2] D[Sinh[ch]/Sinh[ch/Sqrt[2]], ch]];
Gp[c_] := 4 Sqrt[2] (Sinh[c]/Sinh[c/Sqrt[2]] - Sqrt[2]);   (* G'(chi), with G(0) = 0 *)
NRHZeroNR["SM(66)-(67) [NRg, NRGprofile]: d^2 G/d chi^2 = rho(chi)", D[Gp[ch], ch] - rho];
NRH`CheckZero["SM(67) [NRGprofile]: G'(0) = 0  (the G-integrand vanishes at chi = 0)",
   Limit[Gp[ch], ch -> 0]];
NRH`CheckZero["SM(67) [NRGprofile]: d^2/dchi^2 (e^{s chi} - 1 - s chi) = e^{s chi},  s = +1, -1",
   {D[Exp[ch] - 1 - ch, {ch, 2}] - Exp[ch], D[Exp[-ch] - 1 + ch, {ch, 2}] - Exp[-ch]}];


(* ::Section:: *)
(*The full Einstein double field equation with GENERIC W(x^+, x^-, chi)*)


(* ::Text:: *)
(*This is the central bulk computation.  With W a completely generic function of*)
(*(x^+, x^-, chi) we evaluate the full doubled curvature of (13)-(14) and verify:*)
(*  (i)  the scalar equation holds identically:  S_(0) = -4/l^2  for ANY W;*)
(*  (ii) the whole tensor equation (P S Pbar)_MN = 0 collapses to the single radial ODE*)
(*       d^2 W/d chi^2 = F  of SM (64)-(65) - no other independent equation remains;*)
(*  (iii) with the ODE imposed,  G_MN = 2 l^-2 J_MN  holds exactly.*)
(*The source F carries the Hill data  A_pm = psi_pm''/psi_pm  of SM (61) [NRhill].*)


FF = l^2/(16 psip[xp] psim[xm]) (
     rho (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm])
     - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 Exp[ch]
     + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 Exp[-ch]);

curvNR = DFTCurvature[HW, dNR, xsNR];
NRHZeroNR["EDFE scalar: S_(0) = -4/l^2 for ARBITRARY W(x^+, x^-, chi)",
   curvNR["S0"] + 4/l^2];

odeRule = Derivative[0, 0, 2][W][xp, xm, ch] -> FF;
NRH`Check["the tensor equation is not empty: (P S Pbar) contains d^2W/dchi^2",
   ! FreeQ[curvNR["PSPbar"], Derivative[0, 0, 2][W]]];
NRHZeroNR["EDFE tensor: (P S Pbar)_MN = 0  <=>  d^2W/dchi^2 = F   [SM (64)-(65)]",
   curvNR["PSPbar"] /. odeRule];
NRHZeroNR["G_MN = 2 l^-2 J_MN on the ODE shell",
   (curvNR["G"] /. odeRule) - 2/l^2 JJ];


(* ::Section:: *)
(*Eq. (15) [NRWgeneral]: the exact profile solves the ODE, hence the full EDFE*)


Wexact = W0[xp, xm] + W1[xp, xm] ch psip[xp] psim[xm]/2 +
   l^2/(16 psip[xp] psim[xm]) (
      (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm]) GG[ch]
      - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 (Exp[ch] - 1 - ch)
      + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 (Exp[-ch] - 1 + ch));

NRHZeroNR["Eq.(15) [NRWgeneral] solves d^2W/dchi^2 = F  (via G'' = rho)",
   (D[Wexact, {ch, 2}] /. {Derivative[2][GG][ch] -> D[Gp[ch], ch]}) - FF];
NRH`Check["Eq.(15) [NRWgeneral]: W_0 and W_1 multiply the two homogeneous modes {1, chi/(2 Sqrt[Pi])}",
   {D[Wexact, W0[xp, xm]],
    Together[D[Wexact, W1[xp, xm]] - ch psip[xp] psim[xm]/2]} === {1, 0}];
NRH`CheckZero["near the boundary chi = 2 Sqrt[2] q + O(q^3): homogeneous modes ~ {1, e^{-2y/l}}",
   {SeriesCoefficient[chiq, {q, 0, 1}] - 2 Sqrt[2], SeriesCoefficient[chiq, {q, 0, 2}]}];


(* ::Section:: *)
(*SM (70)-(71) [SMBtransform, SMWshift] and SM (82) [SMWisB]: finite B-shifts*)


bshift = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
   {0, bpm, 0, 1, 0, 0}, {-bpm, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};
NRHZeroNR["SM(70)-(71) [SMBtransform, SMWshift]: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched",
   bshift . HNRof[Wf[xp, xm, ch]] . Transpose[bshift] - HNRof[Wf[xp, xm, ch] - 2 bpm]];
NRHZeroNR["SM(82) [SMWisB]: H(W) = Omega_b H(0) Omega_b^T with b_{+-} = -W/2 (pointwise identity)",
   (bshift /. bpm -> -Wf[xp, xm, ch]/2) . HNRof[0] . Transpose[bshift /. bpm -> -Wf[xp, xm, ch]/2]
      - HNRof[Wf[xp, xm, ch]]];
NRH`Check["Omega_b is O(3,3):  Omega J Omega^T = J",
   Together[bshift . JJ . Transpose[bshift] - JJ] === ConstantArray[0, {6, 6}]];


(* ::Section:: *)
(*SM (77) [SNCtau]: SNC clock forms*)


tauP = {Cosh[ch/2], -esig^-1 Sinh[ch/2], 0};   (* components (dx^+, dx^-, dy) *)
tauM = {-esig Sinh[ch/2], Cosh[ch/2], 0};
NRHZeroNR["W-part of the lower-right block = W (tau+ tau- + tau- tau+)  [half angles]",
   (HNRof[Wc][[4 ;; 5, 4 ;; 5]] - HNRof[0][[4 ;; 5, 4 ;; 5]])
      - Wc Table[tauP[[i]] tauM[[j]] + tauM[[i]] tauP[[j]], {i, 2}, {j, 2}]];
NRHZeroNR["unit clock determinant: tau+ ^ tau- = dx+ ^ dx-",
   tauP[[1]] tauM[[2]] - tauP[[2]] tauM[[1]] - 1];
(* Banados-frame asymptotics, expanded in q:  tau+ = dx+ - L_- e^{-2y/l} dx- + O(e^{-4y/l}),
   using e^{-2y/l} = Sqrt[2] q psi_+ psi_-  and  L_- = 1/psi_-^2. *)
tauPq = tauP /. ch -> chiq;
NRH`CheckZero["SM(77) [SNCtau]: tau+ = dx+ - L_- e^{-2y/l} dx- + O(e^{-4y/l})",
   {SeriesCoefficient[tauPq[[1]], {q, 0, 0}] - 1,
    SeriesCoefficient[tauPq[[1]], {q, 0, 1}],
    Together[SeriesCoefficient[tauPq[[2]], {q, 0, 1}] + Sqrt[2] psip[xp]/psim[xm]],
    SeriesCoefficient[tauPq[[2]], {q, 0, 0}],
    SeriesCoefficient[tauPq[[2]], {q, 0, 2}]}];


(* ::Section:: *)
(*Eq. (16) [NRasympt]: asymptotic transformation laws  (near-boundary expansion)*)


(* ::Text:: *)
(*For the asymptotic statement we work in the rational radial variable u = e^{2y/l} with*)
(*W_0 = 0 and the normalizable hair W = W_1(x^+, x^-)/u.  The homogeneous truncation is*)
(*legitimate here because the inhomogeneous parts of Eq. (15) start at O(u^-2), beyond*)
(*every order matched below.  The generator is Eq. (5) with the radial tails of xi^pm*)
(*removed and the dual tails retained, as stated in the Letter.  We verify in the*)
(*1/u expansion (through the normalizable order) that*)
(*   Lhat_xi H = family variation of H under  delta_eps L_pm = eps dL + 2 L d eps,*)
(*                                           delta_eps W_1 = eps^i d_i W_1 + 2 W_1 d_i eps^i*)
(*                                                          - l^2 (L_- d_+^3 eps^+ + L_+ d_-^3 eps^-),*)
(*   Lhat_xi d = O(u^-2):*)
(*so L_pm transforms WITHOUT anomaly and the third-derivative terms sit entirely in W_1.*)


(* To keep every series coefficient rational we parametrize L_pm = psi_pm^{-2}. *)
xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};
qOfu = 1/(Sqrt[2] psip[xp] psim[xm] u);
chOfu = 2 Sqrt[2] ArcTanh[qOfu];
esigU = psim[xm]/psip[xp];
HNRu = {{0, 0, 0, Cosh[chOfu], -Sinh[chOfu]/esigU, 0},
   {0, 0, 0, esigU Sinh[chOfu], -Cosh[chOfu], 0},
   {0, 0, 1, 0, 0, 0},
   {Cosh[chOfu], esigU Sinh[chOfu], 0, -WW esigU Sinh[chOfu], WW Cosh[chOfu], 0},
   {-Sinh[chOfu]/esigU, -Cosh[chOfu], 0, WW Cosh[chOfu], -WW Sinh[chOfu]/esigU, 0},
   {0, 0, 0, 0, 0, 1}} /. WW -> W1[xp, xm]/u;
dNRu = -1/2 Log[u] + Log[Cosh[chOfu/(2 Sqrt[2])]];
LpPsi = 1/psip[xp]^2;  LmPsi = 1/psim[xm]^2;

xiUpNR = {
   -l^2/(2 u) LpPsi D[em[xm], {xm, 2}],
   +l^2/(2 u) LmPsi D[ep[xp], {xp, 2}],
   -l/2 (D[ep[xp], xp] - D[em[xm], xm]),
   ep[xp],
   em[xm],
   -l/2 (D[ep[xp], xp] + D[em[xm], xm])};

lieHNR = GenLieH[xiUpNR, HNRu, xsU];
lieDNR = GenLieD[xiUpNR, dNRu, xsU];

(* Eq. (16) transformation laws, written in the psi parametrization:
   delta L = eps dL + 2 L d eps  <=>  delta psi = eps d psi - psi d eps. *)
dPsiP = ep[xp] D[psip[xp], xp] - psip[xp] D[ep[xp], xp];
dPsiM = em[xm] D[psim[xm], xm] - psim[xm] D[em[xm], xm];
dW1NR = (ep[xp] D[W1[xp, xm], xp] + em[xm] D[W1[xp, xm], xm] +
   2 W1[xp, xm] (D[ep[xp], xp] + D[em[xm], xm]) -
   l^2 (LmPsi D[ep[xp], {xp, 3}] + LpPsi D[em[xm], {xm, 3}]));

HNRuGen = HNRu /. {psip[xp] -> PSPv, psim[xm] -> PSMv, W1[xp, xm] -> W1v};
depsHNR = (D[HNRuGen, PSPv] dPsiP + D[HNRuGen, PSMv] dPsiM + D[HNRuGen, W1v] dW1NR) /.
   {PSPv -> psip[xp], PSMv -> psim[xm], W1v -> W1[xp, xm]};

seriesZero[m_, ord_] := Map[Function[e, Together[Normal[Series[e, {u, Infinity, ord}]]]], m, {2}];
NRH`CheckZero["Eq.(16) [NRasympt]: Lhat_xi H - delta_(L,W1) H = O(u^-2) componentwise",
   seriesZero[lieHNR - depsHNR, 1]];
NRH`CheckZero["Eq.(16) [NRasympt]: delta psi law is equivalent to delta L = eps dL + 2 L d eps (no anomaly)",
   Together[(D[1/PSPv^2, PSPv] dPsiP /. PSPv -> psip[xp])
      - (ep[xp] D[LpPsi, xp] + 2 LpPsi D[ep[xp], xp])]];
NRH`CheckZero["Eq.(16) [NRasympt]: Lhat_xi d = O(u^-2)",
   Together[Normal[Series[lieDNR, {u, Infinity, 1}]]]];
NRH`Check["Eq.(16) [NRasympt]: delta_eps L_pm carries NO third-derivative anomaly",
   FreeQ[{dPsiP, dPsiM}, Derivative[3][_][_]]];


(* ::Section:: *)
(*SM (108) [SMNRlocalstabilizer], second line: the weight-one law of W_0 (two-sided, leading order)*)


(* ::Text:: *)
(*NRH07 verifies the exact one-sided identity with weights (1, 2) for (W_0, W_1).  Here the*)
(*two-sided family with  W = W_0(x) + W_1(x)/u  is transformed by the generator of Eq. (16), and*)
(*the u^0 coefficient of Lhat_xi H is compared with the family variation under*)
(*   delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i      (SM (108), second line);*)
(*at this order neither the W_1 law nor the dual tails contribute.*)


HNRu0 = HNRu /. W1[xp, xm] -> u W0[xp, xm] + W1[xp, xm];   (* W = W_0 + W_1/u *)
lieHNR0 = GenLieH[xiUpNR, HNRu0, xsU];
dW0NR = ep[xp] D[W0[xp, xm], xp] + em[xm] D[W0[xp, xm], xm] + W0[xp, xm] (D[ep[xp], xp] + D[em[xm], xm]);
HNRu0Gen = HNRu0 /. {psip[xp] -> PSPv, psim[xm] -> PSMv, W0[xp, xm] -> W0v};
depsHNR0 = (D[HNRu0Gen, PSPv] dPsiP + D[HNRu0Gen, PSMv] dPsiM + D[HNRu0Gen, W0v] dW0NR) /.
   {PSPv -> psip[xp], PSMv -> psim[xm], W0v -> W0[xp, xm]};
NRH`CheckZero["SM(108) line 2: at O(u^0), Lhat_xi H = family variation with delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i (weight one)",
   seriesZero[lieHNR0 - depsHNR0, 0]];


(* ::Section:: *)
(*Eq. (12) and SM (61) [NRhill]: the sinh(chi/2) expansion and the Hill identity*)


(* ::Text:: *)
(*Below Eq. (12) the Letter states  e^{+-sigma} sinh(chi/2) = e^{-2y/l} L_+- + O(e^{-6y/l}),  and the*)
(*SM (text below SM (106)) writes the same combination as  u L_+- [1 + u^2 Pi/3 + O((u^2 Pi)^2)],*)
(*u = e^{-2y/l}.  With chi = 2 Sqrt[2] arctanh q and q = e^{-2y/l} Sqrt[Pi/2] both are the expansion*)
(*sinh(Sqrt[2] arctanh q) = Sqrt[2] q (1 + (2/3) q^2 + ...):  the q^2 term is absent and (2/3) q^2 = u^2 Pi/3.*)
(*SM (61) states  A_+- = (1/4)(d ln L)^2 - (1/2) d^2 ln L = psi''/psi  for psi = L^{-1/2}.*)
(*The last check records the orders behind the Letter's falloff statement for W:  G(chi) = O(chi^3)*)
(*and e^{s chi} - 1 - s chi = O(chi^2), so the derivative-dependent terms of Eq. (15) start at*)
(*O(e^{-4y/l}) in general and are absent for constant L_pm or in the one-sided limits.*)


sinhHalf = Sinh[chiq/2];
NRH`CheckZero["Eq.(12): e^{sigma} sinh(chi/2) = e^{-2y/l} L_+ + O(e^{-6y/l})  (q^0 and q^2 coefficients vanish; q^1 coefficient matches)",
   {SeriesCoefficient[sinhHalf, {q, 0, 0}], SeriesCoefficient[sinhHalf, {q, 0, 2}],
    Together[esig SeriesCoefficient[sinhHalf, {q, 0, 1}] - Sqrt[2] psip[xp] psim[xm]/psip[xp]^2]}];
NRH`CheckZero["SM text below (106): e^{sigma} sinh(chi/2) = u L_+ [1 + u^2 Pi/3 + ...]  (q^3 coefficient = (2/3) Sqrt[2])",
   Together[SeriesCoefficient[sinhHalf, {q, 0, 3}] - 2/3 Sqrt[2]]];
NRH`CheckZero["SM(61) [NRhill]: (1/4)(d ln L)^2 - (1/2) d^2 ln L = psi''/psi  with psi = L^{-1/2}",
   Together[PowerExpand[1/4 D[Log[LL[x]], x]^2 - 1/2 D[Log[LL[x]], {x, 2}]
      - D[LL[x]^(-1/2), {x, 2}]/LL[x]^(-1/2)]]];
NRH`CheckZero["falloff bookkeeping: G'(chi) = (2/3) chi^2 + O(chi^4)  and  e^{s chi} - 1 - s chi = chi^2/2 + O(chi^3)",
   {SeriesCoefficient[Gp[ch], {ch, 0, 0}], SeriesCoefficient[Gp[ch], {ch, 0, 1}],
    SeriesCoefficient[Gp[ch], {ch, 0, 2}] - 2/3,
    SeriesCoefficient[Exp[ch] - 1 - ch, {ch, 0, 2}] - 1/2}];

NRH`FileSummary[];
