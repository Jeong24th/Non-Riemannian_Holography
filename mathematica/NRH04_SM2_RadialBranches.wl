(* ::Title:: *)
(*NRH04 SM2 Exact Radial Branches and Boundary Data*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH04_SM2_RadialBranches.wl"];

JJ = ODDJ[3];
xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};

(* ::Section:: *)
(*SM2.1 Riemannian saddle*)

gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];
NRH`Check["SM2.1: the Riemannian saddle is of type (0,0): the upper-left block H^{mu nu} = g^{-1} is invertible",
   Together[Det[HR[[1 ;; 3, 1 ;; 3]]]] =!= 0];
NRH`CheckZero["SM2.1: the fixed, L-independent e^{-2y/l} falloff in the type-changing channel H^{+-} (even at L_pm = 0)",
   {SeriesCoefficient[HR[[1, 2]], {u, Infinity, 1}] + 1, SeriesCoefficient[HR[[1, 2]] /. {Lp -> (0 &), Lm -> (0 &)}, {u, Infinity, 1}] + 1}];
NRH`CheckZero["SM2.1: exterior terminates at the Killing horizon e^{4y/l} = L+ L- where e^{-2d} = 0",
   Together[Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]]];

(* ::Section:: *)
(*SM2.2 Non-Riemannian saddle: exact radial equation*)

NRH`CheckZero["NRhill: A = (1/4)(d ln L)^2 - (1/2) d^2 ln L = psi''/psi with psi = L^{-1/2}",
   Together[PowerExpand[1/4 D[Log[LL[x]], x]^2 - 1/2 D[Log[LL[x]], {x, 2}] - D[LL[x]^(-1/2), {x, 2}]/LL[x]^(-1/2)]]];
chiq = 2 Sqrt[2] ArcTanh[q];
NRH`CheckZero["NRradialchange: d chi/d q = 2 Sqrt[2]/(1 - q^2) and d_y = -(2q/l) d_q for q = e^{-2y/l} Sqrt[Pi/2]",
   {Together[D[chiq, q] - 2 Sqrt[2]/(1 - q^2)], Together[D[qq0 Exp[-2 yv/l], yv] + (2/l) qq0 Exp[-2 yv/l]]}];
NRH`CheckZero["NRradialoperator: (1-q^2)^2/8 (d_q^2 - 2q/(1-q^2) d_q) = d_chi^2",
   Together[(1 - q^2)^2/8 (D[ff[chiq], {q, 2}] - 2 q/(1 - q^2) D[ff[chiq], q]) - Derivative[2][ff][chiq]]];
chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
NRHZeroNR[label_, e_] := NRH`CheckZero[label, Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];
NRHZeroNR["d_y chi = -4 Sqrt[Pi] e^{2d}/l (hyperbolic identity form)",
   chy + (4 Sqrt[2]/l) Tanh[ch/(2 Sqrt[2])] Cosh[ch/(2 Sqrt[2])]^2];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
xsNR = {Function[e, D[e, xp] + chp D[e, ch]], Function[e, D[e, xm] + chm D[e, ch]], Function[e, chy D[e, ch] + D[e, Ysym]]};
esig = psim[xm]/psip[xp];
HNR = NonRiemannianH[ch, esig, W[xp, xm, ch]];
dNR = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];
rho = Together[4 Sqrt[2] D[Sinh[ch]/Sinh[ch/Sqrt[2]], ch]];
radialSource = l^2/(16 psip[xp] psim[xm]) (rho (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm])
   - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 Exp[ch] + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 Exp[-ch]);
curvNR = DFTCurvature[HNR, dNR, xsNR];
NRH`Check["NRchiODE: the tensor EDFE contains d^2 W/d chi^2", ! FreeQ[curvNR["PSPbar"], Derivative[0, 0, 2][W]]];
NRHZeroNR["NRchiODE, NRsource: (P S Pbar)_MN = 0 <=> d^2W/dchi^2 = F with the displayed source", curvNR["PSPbar"] /. Derivative[0, 0, 2][W][xp, xm, ch] -> radialSource];
NRHZeroNR["EDFE scalar: S_(0) = -4/l^2 for arbitrary W", curvNR["S0"] + 4/l^2];
Gp[c_] := 4 Sqrt[2] (Sinh[c]/Sinh[c/Sqrt[2]] - Sqrt[2]);  (* Gp = I_radial', GG = I_radial: manuscript calligraphic I. *)
NRHZeroNR["NRg, NRGprofile: d^2 I_radial/d chi^2 = rho(chi)", D[Gp[ch], ch] - rho];
NRH`CheckZero["NRGprofile: I_radial'(0) = 0 and I_radial'(chi) = (2/3) chi^2 + O(chi^4); the integral normalization I_radial(0) = 0 gives I_radial = (2/9) chi^3 + O(chi^5)", {Limit[Gp[ch], ch -> 0], Normal[Series[Gp[ch], {ch, 0, 3}]] - 2/3 ch^2}];
NRH`CheckZero["NRGprofile: d^2/dchi^2 (e^{s chi} - 1 - s chi) = e^{s chi}", {D[Exp[ch] - 1 - ch, {ch, 2}] - Exp[ch], D[Exp[-ch] - 1 + ch, {ch, 2}] - Exp[-ch]}];
Wexact = W0[xp, xm] + W1[xp, xm] ch psip[xp] psim[xm]/2 + l^2/(16 psip[xp] psim[xm]) (
   (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm]) GG[ch]
   - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 (Exp[ch] - 1 - ch)
   + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 (Exp[-ch] - 1 + ch));
NRHZeroNR["NRWgeneral: twice integrating NRchiODE gives the exact solution (with I_radial'' = rho)",
   (D[Wexact, {ch, 2}] /. Derivative[2][GG][ch] -> D[Gp[ch], ch]) - radialSource];
NRH`Check["the two integration functions W_0, W_1 multiply {1, chi/(2 Sqrt[Pi])}",
   {D[Wexact, W0[xp, xm]], Together[D[Wexact, W1[xp, xm]] - ch psip[xp] psim[xm]/2]} === {1, 0}];
NRH`CheckZero["constant L: q = Sqrt[L+L-/2] mu_RG^{-2}; mu d chi/d mu = -2 Sqrt[2] sinh(chi/Sqrt[2])",
   Module[{qm = Sqrt[L0p L0m/2] muRG^-2, chim},
      chim = 2 Sqrt[2] ArcTanh[qm];
      Simplify[muRG D[chim, muRG] + 2 Sqrt[2] Sinh[chim/Sqrt[2]], L0p > 0 && L0m > 0 && muRG > 0 && qm < 1]]];
NRH`CheckZero["the endpoint q = 1 lies at e^{4y/l} = L+ L-/2", Together[(Sqrt[L0p L0m/2] muRG^-2 /. muRG -> (L0p L0m/2)^(1/4)) - 1]];

(* ::Section:: *)
(*SM2.2 The radial zero mode W_0 is locally pure generalized gauge*)

bshift = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, bpm, 0, 1, 0, 0}, {-bpm, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};
NRHZeroNR["SMBtransform, SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly, d untouched",
   bshift . NonRiemannianH[ch, esig, Wf[xp, xm, ch]] . Transpose[bshift] - NonRiemannianH[ch, esig, Wf[xp, xm, ch] - 2 bpm]];
NRH`Check["SMBtransform: Omega_b is O(3,3)", Together[bshift . JJ . Transpose[bshift] - JJ] === ConstantArray[0, {6, 6}]];
NRH`CheckZero["b = (1/2) W_0 dx+ ^ dx- = d lambda with lambda_- = (1/2) Int^{x+} W_0, lambda_+ = lambda_y = 0",
   Module[{lam = {0, 1/2 Integrate[W0[s, xm], {s, 0, xp}], 0}, b},
      b = Table[D[lam[[j]], {xp, xm, yv}[[i]]] - D[lam[[i]], {xp, xm, yv}[[j]]], {i, 3}, {j, 3}];
      {b[[1, 2]] - W0[xp, xm]/2, b[[1, 3]], b[[2, 3]]}]];
NRH`FileSummary[];
