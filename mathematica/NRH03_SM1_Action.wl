(* ::Title:: *)
(*NRH03 SM1 Bulk Action and On-Shell Saddle Values*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH03_SM1_Action.wl"];

JJ = ODDJ[3];
xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};

(* SMgamma2: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M) on both saddles *)
gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];
gammaR = GammaDFT[HR, dR, xsU];
NRH`CheckZero["SMgamma2 on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M) for arbitrary chiral L_pm",
   Together[Exp[-2 dR] ScalarS0[HR, dR, xsU] - Gamma2Density[HR, dR, gammaR, xsU]
      - Sum[DblD[Exp[-2 dR] GammaBVector[HR, dR, xsU][[m]], m, xsU], {m, 6}]]];
BvecR = GammaBVector[HR, dR, xsU];
NRH`CheckZero["SMgamma2flux on R: B^y = 4 d_y d and e^{-2d} B^y = -(4/l)(e^{2y/l} + L+L- e^{-2y/l})",
   {Together[BvecR[[6]] - 4 (2 u/l) D[dR, u]], Together[Exp[-2 dR] BvecR[[6]] + 4/l (u + Lp[xp] Lm[xm]/u)]}];
NRH`Check["SMgamma2 on R: the section fluxes B^{x pm} are built from the W-independent blocks and are total x-derivatives of periodic data",
   FreeQ[Together[BvecR[[4 ;; 5]]], W] ];

chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
xsNR = {Function[e, D[e, xp] + chp D[e, ch]], Function[e, D[e, xm] + chm D[e, ch]], Function[e, chy D[e, ch] + D[e, Ysym]]};
esig = psim[xm]/psip[xp];
HNRchi = NonRiemannianH[ch, esig, W[xp, xm, ch]];
dNRchi = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];
NRHZeroNR[label_, e_] := NRH`CheckZero[label, Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];
gammaNRc = GammaDFT[HNRchi, dNRchi, xsNR];
NRHZeroNR["SMgamma2 on NR (arbitrary W): e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)",
   Exp[-2 dNRchi] ScalarS0[HNRchi, dNRchi, xsNR] - Gamma2Density[HNRchi, dNRchi, gammaNRc, xsNR]
      - Sum[DblD[Exp[-2 dNRchi] GammaBVector[HNRchi, dNRchi, xsNR][[m]], m, xsNR], {m, 6}]];
BvecNR = GammaBVector[HNRchi, dNRchi, xsNR];
NRH`Check["SMgamma2 on NR: B^pm and B^y contain no W (the hair never enters the flux)", FreeQ[Together[BvecNR], W]];
NRHZeroNR["SMgamma2flux on NR: B^y = 4 d_y d", BvecNR[[6]] - 4 (chy D[dNRchi, ch] + D[dNRchi, Ysym])];
eNRu = u - Lp[xp] Lm[xm]/(2 u);
NRH`CheckZero["SMmudefinition on NR: e^{-2d} = e^{2y/l}(1 - q^2) = u - (L+L-/2)/u with q = Sqrt[L+L-/2] e^{-2y/l}, so -2 d_y e^{-2d} = -(4/l)(u + mu/u) with mu = L+L-/2",
   Module[{em2d = Exp[-2 dNRchi] /. ch -> 2 Sqrt[2] ArcTanh[qq] /. Ysym -> (l/2) Log[u] /. qq -> Sqrt[Lp[xp] Lm[xm]/2]/u},
      {Simplify[em2d - eNRu, Assumptions -> u > 0 && l > 0 && Lp[xp] > 0 && Lm[xm] > 0],
       Together[-2 (2 u/l) D[eNRu, u] + 4/l (u + Lp[xp] Lm[xm]/(2 u))]}]];
NRH`CheckZero["SMgamma2flux, SMmudefinition: the same formula e^{-2d} B^y = -(4/l)(e^{2y/l} + mu e^{-2y/l}) holds on both branches",
   {Together[Exp[-2 dR] BvecR[[6]] + 4/l (u + (Lp[xp] Lm[xm])/u)], Together[-2 (2 u/l) D[eNRu, u] + 4/l (u + (Lp[xp] Lm[xm]/2)/u)]}];
SrenY = 1/(16 Pi G) (4/l (Exp[2 Y/l] + mu Exp[-2 Y/l]) - 8/l Sqrt[mu] - 4/l (Exp[2 Y/l] - mu Exp[-2 Y/l]));
NRH`CheckZero["SMgamma2cutoff: the regulated flux plus the volume counterterm equals (16 pi G)^{-1}[(8 mu/l) e^{-2Y/l} - (8/l) Sqrt[mu]]",
   Together[SrenY - 1/(16 Pi G) (8 mu/l Exp[-2 Y/l] - 8/l Sqrt[mu])]];
NRH`CheckZero["SMgamma2cutoff: the interior endpoint e^{2y*/l} = Sqrt[mu] contributes (4/l)(Sqrt[mu] + Sqrt[mu]) = (8/l) Sqrt[mu]",
   Together[4/l (Exp[2 ys/l] + mu Exp[-2 ys/l]) - 8/l Sqrt[mu] /. ys -> l/2 Log[Sqrt[mu]]] /. Sqrt[mu]^2 -> mu];
NRH`CheckZero["SMgamma2value: Y -> Infinity gives S_ren = -(8 Sqrt[mu])/(16 pi G l) Int d^2x",
   Limit[SrenY, Y -> Infinity, Assumptions -> l > 0 && mu > 0] + 8 Sqrt[mu]/(16 Pi G l)];
NRH`CheckZero["endpoints: e^{-2d} = 0 at u^2 = L+L- (Riemannian horizon) and at u^2 = L+L-/2 (non-Riemannian q = 1)",
   {Together[Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]], Together[eNRu /. u -> Sqrt[Lp[xp] Lm[xm]/2]]}];
NRH`Check["the hair costs no on-shell action: neither the on-shell density, the flux nor the counterterm involves W",
   FreeQ[Together[Exp[-2 dNRchi] BvecNR[[6]]], W] && FreeQ[Together[Exp[-2 dNRchi] (-4/l^2)], W]];

NRH`FileSummary[];
