(* ::Title:: *)
(*NRH02 Letter Riemannian*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH02_Letter_Riemannian.wl"];

xs = {xp, xm, Function[e, (2 u/l) D[e, u]]};
JJ = ODDJ[3];
fR = u + Lp[xp] Lm[xm]/u;

gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];

NRH`CheckZero["H J H = J (O(3,3) constraint)",
   Map[Together, HR . JJ . HR - JJ, {2}]];
NRH`CheckZero["(e^{-2d})^2 = -det g (phi_0 = 0)",
   Together[Exp[-2 dR]^2 + Det[gR]]];
NRH`CheckZero["e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l})",
   Together[Exp[-2 dR] - u (1 - Lp[xp] Lm[xm]/u^2)]];

curvR = DFTCurvature[HR, dR, xs];
NRH`CheckZero["EDFE tensor part: (P S Pbar)_MN = 0, arbitrary chiral L_pm", curvR["PSPbar"]];
NRH`CheckZero["EDFE scalar part: S_(0) = -4/l^2", Together[curvR["S0"] + 4/l^2]];
NRH`CheckZero["S_(0) closed form = S_(0) from S_ABCD contraction",
   Together[curvR["S0"] - ScalarS0FromS4[curvR["Gamma"], curvR["R4"], HR, xs]]];
NRH`CheckZero["G_MN = 2 l^-2 J_MN", Map[Together, curvR["G"] - 2/l^2 JJ, {2}]];

Module[{P = (JJ + HR)/2, gamma = curvR["Gamma"], compat, tr},
   compat = Table[
      DblD[P, c, xs][[a, b]]
      + Sum[(gamma[[c]] . JJ)[[a, dd]] P[[dd, b]], {dd, 6}]
      + Sum[(gamma[[c]] . JJ)[[b, dd]] P[[a, dd]], {dd, 6}],
      {c, 6}, {a, 6}, {b, 6}];
   NRH`CheckZero["Gamma compatibility: nabla_C P_AB = 0", Map[Together, compat, {3}]];
   tr = Table[Sum[JJ[[b, e]] gamma[[e, b, a]], {b, 6}, {e, 6}], {a, 6}];
   NRH`CheckZero["Gamma dilaton trace: Gamma^B_{BA} = -2 partial_A d",
      Together[tr + 2 DblGrad[dR, xs]]];
   NRH`CheckZero["Gamma torsionless: cyclic sum Gamma_{[CAB]} = 0",
      Map[Together, Table[gamma[[c, a, b]] + gamma[[a, b, c]] + gamma[[b, c, a]], {c, 6}, {a, 6}, {b, 6}], {3}]];
];

NRH`CheckZero["boundaryH: H^infty matches the displayed 6x6 matrix",
   Map[Limit[#, u -> Infinity] &, HR, {2}] - Hinf];
NRH`CheckZero["H - H^infty = O(e^{-2y/l}) (every entry vanishes at the boundary)",
   Map[Limit[#, u -> Infinity] &, HR - Hinf, {2}]];
NRH`CheckZero["d + y/l = O(e^{-4y/l}): u * (d + y/l) still vanishes at the boundary",
   Limit[u (dR + 1/2 Log[u]), u -> Infinity]];

H0 = Hinf[[{1, 2, 4, 5}, {1, 2, 4, 5}]];
NRH`Check["boundaryH: induced boundary H^(0) has vanishing upper-left block (type (1,1))",
   H0[[1 ;; 2, 1 ;; 2]] === {{0, 0}, {0, 0}}];

J4 = ODDJ[2];
V0 = {{1/Sqrt[2], 0}, {0, 0}, {0, -Sqrt[2]}, {0, 0}};
Vb0 = {{0, 0}, {0, 1/Sqrt[2]}, {0, 0}, {Sqrt[2], 0}};
eta2 = {{0, -1}, {-1, 0}};   etab2 = {{0, 1}, {1, 0}};

NRH`CheckZero["Rboundaryframe: boundary vielbeins reproduce P^(0) and Pbar^(0)",
   {V0 . eta2 . Transpose[V0] - (J4 + H0)/2, Vb0 . etab2 . Transpose[Vb0] - (J4 - H0)/2}];

kmat = {{Kpp[xp, xm], Kpm[xp, xm]}, {Kmp[xp, xm], Kmm[xp, xm]}};
T0f = T0[xp, xm];

V0u = V0 . eta2;  Vb0u = Vb0 . etab2;
TAB = Table[
   2 Sum[(V0u[[a, i]] Vb0u[[b, j]] - V0u[[b, i]] Vb0u[[a, j]]) kmat[[i, j]], {i, 2}, {j, 2}]
   - 1/2 J4[[a, b]] T0f, {a, 4}, {b, 4}];
bD[e_, a_] := If[a <= 2, 0, D[e, {xp, xm}[[a - 2]]]];
divT = Table[Together[Sum[J4[[a, c]] bD[TAB[[a, b]], c], {a, 4}, {c, 4}]], {b, 4}];
ward1 = D[Kpp[xp, xm], xm] + 1/4 D[T0f, xp];
ward2 = D[Kmm[xp, xm], xp] + 1/4 D[T0f, xm];

NRH`CheckZero["div T = {d_- K_mp, -d_+ K_mp, -2 Ward_+, -2 Ward_-} exactly",
   Together[divT - {D[Kmp[xp, xm], xm], -D[Kmp[xp, xm], xp], -2 ward1, -2 ward2}]];
NRH`Check["no local condition on K_{op bom}", FreeQ[divT, Kpm]];

xiUp = {
   -l^2/(2 u) Lp[xp] D[em[xm], {xm, 2}],
   +l^2/(2 u) Lm[xm] D[ep[xp], {xp, 2}],
   -l/2 (D[ep[xp], xp] - D[em[xm], xm]),
   ep[xp] + l^2/(4 u) D[em[xm], {xm, 2}],
   em[xm] + l^2/(4 u) D[ep[xp], {xp, 2}],
   -l/2 (D[ep[xp], xp] + D[em[xm], xm])};

lieH = Map[Together, GenLieH[xiUp, HR, xs], {2}];
lieD = Together[GenLieD[xiUp, dR, xs]];

dLp = ep[xp] D[Lp[xp], xp] + 2 Lp[xp] D[ep[xp], xp] - l^2/4 D[ep[xp], {xp, 3}];
dLm = em[xm] D[Lm[xm], xm] + 2 Lm[xm] D[em[xm], xm] - l^2/4 D[em[xm], {xm, 3}];

HRgen = HR /. {Lp[xp] -> LPv, Lm[xm] -> LMv};
depsH = Map[Together,
   (D[HRgen, LPv] dLp + D[HRgen, LMv] dLm) /. {LPv -> Lp[xp], LMv -> Lm[xm]}, {2}];

orderCheck[m_, pow_] := Map[Function[e, Normal[Series[e, {u, Infinity, pow}]]], m, {2}];
NRH`CheckZero["Lhat_xi d = O(e^{-4y/l})",
   Normal[Series[lieD, {u, Infinity, 1}]]];
NRH`CheckZero["Lhat_xi H - delta_eps H = O(e^{-4y/l})",
   orderCheck[lieH - depsH, 1]];
NRH`CheckZero["delta_eps H = O(e^{-2y/l}) (leading falloff of the Banados variation)",
   orderCheck[depsH, 0]];

THol = 8 Pi (Lp[xp]/(16 Pi G l));
dT = (dLp/(2 G l)) /. {};
centralCharge = 3 l/(2 G);
NRH`CheckZero["c = 3l/2G reproduces delta_eps T = eps T' + 2 T eps' - (c/12) eps'''",
   Together[dT - (ep[xp] D[THol, xp] + 2 THol D[ep[xp], xp] - centralCharge/12 D[ep[xp], {xp, 3}])]];
NRH`CheckZero["(64 pi G l)^{-1} (3 l^2/(4 pi)) = (8 pi)^{-2} (c/2)",
   Together[1/(64 Pi G l) 3 l^2/(4 Pi) - 1/(8 Pi)^2 centralCharge/2]];

NRH`CheckZero["horizon: e^{-2d} = 0 at e^{4y/l} = L+ L- (u^2 = L+ L-)",
   Together[(Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]) ]];
EE2 = Map[Together, gR - bR, {2}];
NRH`Check["E_{+-} = 0, E_{-+} = -2F, E_{pm pm} = 2 L_pm (Gomis-Ooguri channel structure)",
   And[Together[EE2[[1, 2]]] === 0,
       Together[EE2[[2, 1]] + 2 fR] === 0,
       Together[EE2[[1, 1]] - 2 Lp[xp]] === 0,
       Together[EE2[[2, 2]] - 2 Lm[xm]] === 0]];
NRH`CheckZero["c_eff^2 := 2F = 2(e^{2y/l} + L+ L- e^{-2y/l}) and = 4 Sqrt[L+ L-] at the horizon",
   Together[(2 fR /. u -> Sqrt[Lp[xp] Lm[xm]]) - 4 Sqrt[Lp[xp] Lm[xm]]]];

NRH`CheckZero["RKdef: coefficient matching -2K = (16 pi G)^{-1} A^y and 2T_(0) = (16 pi G)^{-1} 2(B^y + 4/l)",
   {Together[-2 (-(1/(32 Pi G))) - 1/(16 Pi G)], Together[2 (1/(16 Pi G)) - 2/(16 Pi G)]}];
NRH`CheckZero["SMinvariantstress: T_{x+ x~-} = 2 K_{op bop} and T_{x- x~+} = 2 K_{om bom} in the frame of",
   {Together[TAB[[3, 2]] - 2 Kpp[xp, xm]], Together[TAB[[4, 1]] - 2 Kmm[xp, xm]]}];
NRH`Check["J_{x+ x~-} = 0 = J_{x- x~+}, so no T_(0) term enters these two components",
   J4[[3, 2]] === 0 && J4[[4, 1]] === 0];

AyR = MomentumAK[HR, curvR["Gamma"], xs][[6]];
AfixR = Map[Together, Transpose[JJ . Vinf] . AyR . (JJ . Vbinf), {2}];
KR = Map[Limit[-(1/(32 Pi G)) u #, u -> Infinity] &, AfixR[[1 ;; 2, 1 ;; 2]], {2}];
NRH`CheckZero["on the exact family: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) [, , mixed entries after ]",
   Map[Together, KR - {{Lp[xp], Lp[xp] Lm[xm]}, {1, Lm[xm]}}/(16 Pi G l), {2}]];
NRH`CheckZero["on the exact family: e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0",
   Limit[u (GammaBVector[HR, dR, xs][[6]] + 4/l), u -> Infinity]];

NRH`FileSummary[];
