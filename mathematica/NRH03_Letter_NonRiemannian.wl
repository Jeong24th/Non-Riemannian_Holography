(* ::Title:: *)
(*NRH03 Letter NonRiemannian*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH03_Letter_NonRiemannian.wl"];

JJ = ODDJ[3];

chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
dpOp = Function[e, D[e, xp] + chp D[e, ch]];
dmOp = Function[e, D[e, xm] + chm D[e, ch]];
dyOp = Function[e, chy D[e, ch] + D[e, Ysym]];
xsNR = {dpOp, dmOp, dyOp};

NRHZeroNR[label_, e_] := NRH`CheckZero[label,
   Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];

esig = psim[xm]/psip[xp];
HNRof[Wexpr_] := NonRiemannianH[ch, esig, Wexpr];
dNR = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];

HNR = HNRof[W[xp, xm, ch]];
NRHZeroNR["H J H = J for the exact non-Riemannian matrix (generic W)",
   HNR . JJ . HNR - JJ];
NRH`Check["type (1,1) at every radius: upper-left block = diag(0,0,1)",
   HNR[[1 ;; 3, 1 ;; 3]] === {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}];
NRHZeroNR["e^{-2d} e^{-2y/l} = 1 - q^2 with q = tanh(chi/(2 Sqrt[2]))",
   Exp[-2 (dNR + Ysym/l)] - (1 - Tanh[ch/(2 Sqrt[2])]^2)];

NRHZeroNR["d chi/dy = -(4 Sqrt[Pi]/l) e^{2d} (hyperbolic identity form)",
   chy + (4 Sqrt[2]/l) Tanh[ch/(2 Sqrt[2])] Cosh[ch/(2 Sqrt[2])]^2];
NRHZeroNR["SM: RG rapidity mu d chi/d mu = -2 Sqrt[2] Sinh[chi/Sqrt[2]]",
   l chy + 2 Sqrt[2] Sinh[ch/Sqrt[2]]];

chiq = 2 Sqrt[2] ArcTanh[q];
NRH`CheckZero["NRradialchange: d chi/d q = 2 Sqrt[2]/(1-q^2)",
   Together[D[chiq, q] - 2 Sqrt[2]/(1 - q^2)]];
NRH`CheckZero["NRradialoperator: transformed radial operator equals f''(chi(q))",
   Together[(1 - q^2)^2/8 (D[ff[chiq], {q, 2}] - 2 q/(1 - q^2) D[ff[chiq], q])
      - Derivative[2][ff][chiq]]];
rho = Together[4 Sqrt[2] D[Sinh[ch]/Sinh[ch/Sqrt[2]], ch]];
Gp[c_] := 4 Sqrt[2] (Sinh[c]/Sinh[c/Sqrt[2]] - Sqrt[2]);
NRHZeroNR["NRg, NRGprofile: d^2 G/d chi^2 = rho(chi)", D[Gp[ch], ch] - rho];
NRH`CheckZero["NRGprofile: G'(0) = 0 (the G-integrand vanishes at chi = 0)",
   Limit[Gp[ch], ch -> 0]];
NRH`CheckZero["NRGprofile: d^2/dchi^2 (e^{s chi} - 1 - s chi) = e^{s chi}, s = +1, -1",
   {D[Exp[ch] - 1 - ch, {ch, 2}] - Exp[ch], D[Exp[-ch] - 1 + ch, {ch, 2}] - Exp[-ch]}];

radialSource = l^2/(16 psip[xp] psim[xm]) (
     rho (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm])
     - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 Exp[ch]
     + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 Exp[-ch]);

curvNR = DFTCurvature[HNR, dNR, xsNR];
NRHZeroNR["EDFE scalar: S_(0) = -4/l^2 for ARBITRARY W(x^+, x^-, chi)",
   curvNR["S0"] + 4/l^2];

odeRule = Derivative[0, 0, 2][W][xp, xm, ch] -> radialSource;
NRH`Check["the tensor equation is not empty: (P S Pbar) contains d^2W/dchi^2",
   ! FreeQ[curvNR["PSPbar"], Derivative[0, 0, 2][W]]];
NRHZeroNR["EDFE tensor: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F []",
   curvNR["PSPbar"] /. odeRule];
NRHZeroNR["G_MN = 2 l^-2 J_MN on the ODE shell",
   (curvNR["G"] /. odeRule) - 2/l^2 JJ];

Wexact = W0[xp, xm] + W1[xp, xm] ch psip[xp] psim[xm]/2 +
   l^2/(16 psip[xp] psim[xm]) (
      (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm]) GG[ch]
      - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 (Exp[ch] - 1 - ch)
      + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 (Exp[-ch] - 1 + ch));

NRHZeroNR["NRWgeneral: solves d^2W/dchi^2 = F (via G'' = rho)",
   (D[Wexact, {ch, 2}] /. {Derivative[2][GG][ch] -> D[Gp[ch], ch]}) - radialSource];
NRH`Check["NRWgeneral: W_0 and W_1 multiply the two homogeneous modes {1, chi/(2 Sqrt[Pi])}",
   {D[Wexact, W0[xp, xm]],
    Together[D[Wexact, W1[xp, xm]] - ch psip[xp] psim[xm]/2]} === {1, 0}];
NRH`CheckZero["near the boundary chi = 2 Sqrt[2] q + O(q^3): homogeneous modes ~ {1, e^{-2y/l}}",
   {SeriesCoefficient[chiq, {q, 0, 1}] - 2 Sqrt[2], SeriesCoefficient[chiq, {q, 0, 2}]}];

bshift = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
   {0, bpm, 0, 1, 0, 0}, {-bpm, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};
NRHZeroNR["SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched",
   bshift . HNRof[Wf[xp, xm, ch]] . Transpose[bshift] - HNRof[Wf[xp, xm, ch] - 2 bpm]];
NRHZeroNR["SMWisB: H(W) = Omega_b H(0) Omega_b^T with b_{+-} = -W/2 (pointwise identity)",
   (bshift /. bpm -> -Wf[xp, xm, ch]/2) . HNRof[0] . Transpose[bshift /. bpm -> -Wf[xp, xm, ch]/2]
      - HNRof[Wf[xp, xm, ch]]];
NRH`Check["Omega_b is O(3,3): Omega J Omega^T = J",
   Together[bshift . JJ . Transpose[bshift] - JJ] === ConstantArray[0, {6, 6}]];

tauP = {Cosh[ch/2], -esig^-1 Sinh[ch/2], 0};
tauM = {-esig Sinh[ch/2], Cosh[ch/2], 0};
NRHZeroNR["W-part of the lower-right block = W (tau+ tau- + tau- tau+) [half angles]",
   (HNRof[Wc][[4 ;; 5, 4 ;; 5]] - HNRof[0][[4 ;; 5, 4 ;; 5]])
      - Wc Table[tauP[[i]] tauM[[j]] + tauM[[i]] tauP[[j]], {i, 2}, {j, 2}]];
NRHZeroNR["unit clock determinant: tau+ ^ tau- = dx+ ^ dx-",
   tauP[[1]] tauM[[2]] - tauP[[2]] tauM[[1]] - 1];

tauPq = tauP /. ch -> chiq;
NRH`CheckZero["SNCtau: tau+ = dx+ - L_- e^{-2y/l} dx- + O(e^{-4y/l})",
   {SeriesCoefficient[tauPq[[1]], {q, 0, 0}] - 1,
    SeriesCoefficient[tauPq[[1]], {q, 0, 1}],
    Together[SeriesCoefficient[tauPq[[2]], {q, 0, 1}] + Sqrt[2] psip[xp]/psim[xm]],
    SeriesCoefficient[tauPq[[2]], {q, 0, 0}],
    SeriesCoefficient[tauPq[[2]], {q, 0, 2}]}];

xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};
qu = 1/(Sqrt[2] psip[xp] psim[xm] u);
chU = 2 Sqrt[2] ArcTanh[qu];

HNRu = NonRiemannianH[chU, esig, W1[xp, xm]/u];
dNRu = -1/2 Log[u] + Log[Cosh[chU/(2 Sqrt[2])]];


xiUpNR = {
   -l^2/(2 u) (1/psip[xp]^2) D[em[xm], {xm, 2}],
   +l^2/(2 u) (1/psim[xm]^2) D[ep[xp], {xp, 2}],
   -l/2 (D[ep[xp], xp] - D[em[xm], xm]),
   ep[xp],
   em[xm],
   -l/2 (D[ep[xp], xp] + D[em[xm], xm])};

lieHNR = GenLieH[xiUpNR, HNRu, xsU];
lieDNR = GenLieD[xiUpNR, dNRu, xsU];

dPsiP = ep[xp] D[psip[xp], xp] - psip[xp] D[ep[xp], xp];
dPsiM = em[xm] D[psim[xm], xm] - psim[xm] D[em[xm], xm];
dW1NR = (ep[xp] D[W1[xp, xm], xp] + em[xm] D[W1[xp, xm], xm] +
   2 W1[xp, xm] (D[ep[xp], xp] + D[em[xm], xm]) -
   l^2 ((1/psim[xm]^2) D[ep[xp], {xp, 3}] + (1/psip[xp]^2) D[em[xm], {xm, 3}]));

HNRuGen = HNRu /. {psip[xp] -> PSPv, psim[xm] -> PSMv, W1[xp, xm] -> W1v};
depsHNR = (D[HNRuGen, PSPv] dPsiP + D[HNRuGen, PSMv] dPsiM + D[HNRuGen, W1v] dW1NR) /.
   {PSPv -> psip[xp], PSMv -> psim[xm], W1v -> W1[xp, xm]};

seriesZero[m_, ord_] := Map[Function[e, Together[Normal[Series[e, {u, Infinity, ord}]]]], m, {2}];
NRH`CheckZero["NRasympt: Lhat_xi H - delta_(L,W1) H = O(u^-2) componentwise",
   seriesZero[lieHNR - depsHNR, 1]];
NRH`CheckZero["NRasympt: delta psi law is equivalent to delta L = eps dL + 2 L d eps (no anomaly)",
   Together[(D[1/PSPv^2, PSPv] dPsiP /. PSPv -> psip[xp])
      - (ep[xp] D[(1/psip[xp]^2), xp] + 2 (1/psip[xp]^2) D[ep[xp], xp])]];
NRH`CheckZero["NRasympt: Lhat_xi d = O(u^-2)",
   Together[Normal[Series[lieDNR, {u, Infinity, 1}]]]];
NRH`Check["NRasympt: delta_eps L_pm carries NO third-derivative anomaly",
   FreeQ[{dPsiP, dPsiM}, Derivative[3][_][_]]];

HNRu0 = HNRu /. W1[xp, xm] -> u W0[xp, xm] + W1[xp, xm];
lieHNR0 = GenLieH[xiUpNR, HNRu0, xsU];
dW0NR = ep[xp] D[W0[xp, xm], xp] + em[xm] D[W0[xp, xm], xm] + W0[xp, xm] (D[ep[xp], xp] + D[em[xm], xm]);
HNRu0Gen = HNRu0 /. {psip[xp] -> PSPv, psim[xm] -> PSMv, W0[xp, xm] -> W0v};
depsHNR0 = (D[HNRu0Gen, PSPv] dPsiP + D[HNRu0Gen, PSMv] dPsiM + D[HNRu0Gen, W0v] dW0NR) /.
   {PSPv -> psip[xp], PSMv -> psim[xm], W0v -> W0[xp, xm]};
NRH`CheckZero["line 2: at O(u^0), Lhat_xi H = family variation with delta W_0 = eps^i d_i W_0 + W_0 d_i eps^i (weight one)",
   seriesZero[lieHNR0 - depsHNR0, 0]];

sinhHalf = Sinh[chiq/2];
NRH`CheckZero["e^{sigma} sinh(chi/2) = e^{-2y/l} L_+ + O(e^{-6y/l}) (q^0 and q^2 coefficients vanish; q^1 coefficient matches)",
   {SeriesCoefficient[sinhHalf, {q, 0, 0}], SeriesCoefficient[sinhHalf, {q, 0, 2}],
    Together[esig SeriesCoefficient[sinhHalf, {q, 0, 1}] - Sqrt[2] psip[xp] psim[xm]/psip[xp]^2]}];
NRH`CheckZero["SM text below (106): e^{sigma} sinh(chi/2) = u L_+ [1 + u^2 Pi/3 + ...] (q^3 coefficient = (2/3) Sqrt[2])",
   Together[SeriesCoefficient[sinhHalf, {q, 0, 3}] - 2/3 Sqrt[2]]];
NRH`CheckZero["NRhill: (1/4)(d ln L)^2 - (1/2) d^2 ln L = psi''/psi with psi = L^{-1/2}",
   Together[PowerExpand[1/4 D[Log[LL[x]], x]^2 - 1/2 D[Log[LL[x]], {x, 2}]
      - D[LL[x]^(-1/2), {x, 2}]/LL[x]^(-1/2)]]];
NRH`CheckZero["falloff bookkeeping: G'(chi) = (2/3) chi^2 + O(chi^4) and e^{s chi} - 1 - s chi = chi^2/2 + O(chi^3)",
   {SeriesCoefficient[Gp[ch], {ch, 0, 0}], SeriesCoefficient[Gp[ch], {ch, 0, 1}],
    SeriesCoefficient[Gp[ch], {ch, 0, 2}] - 2/3,
    SeriesCoefficient[Exp[ch] - 1 - ch, {ch, 0, 2}] - 1/2}];

NRH`FileSummary[];
