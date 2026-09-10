(* ::Title:: *)
(*NRH06 SM Worldsheet*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH06_SM_Worldsheet.wl"];

fR = u + Lp Lm/u;
L1 = dy by + 2 Lp dxp bxp + 2 Lm dxm bxm + beta bxp + betab dxm + beta betab/(2 fR);
betaSol = First@Solve[{D[L1, beta] == 0, D[L1, betab] == 0}, {beta, betab}];
NRH`Check["the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+",
   Together[(beta /. betaSol) + 2 fR dxm] === 0 && Together[(betab /. betaSol) + 2 fR bxp] === 0];
L2 = Together[L1 /. betaSol];
EmatR = {{2 Lp, 0, 0}, {-2 fR, 2 Lm, 0}, {0, 0, 1}};
LE = Sum[EmatR[[m, n]] {dxp, dxm, dy}[[m]] {bxp, bxm, by}[[n]], {m, 3}, {n, 3}];
NRH`CheckZero["eliminating the auxiliaries reproduces E_{mu nu} dx^mu dbar x^nu, E = g - B",
   Together[L2 - LE]];
NRH`CheckZero["-det g_par = F^2 - 4 L+ L- = e^{-4 d_R}",
   Together[-Det[{{2 Lp, -fR}, {-fR, 2 Lm}}] - (fR^2 - 4 Lp Lm)]];
NRH`CheckZero["c_eff^2 = 2F -> 4 Sqrt[L+L-] at the horizon u = Sqrt[L+L-]",
   Together[(2 fR /. u -> Sqrt[Lp Lm]) - 4 Sqrt[Lp Lm]]];

tauP = {Cosh[ch/2], -esig^-1 Sinh[ch/2], 0};
tauM = {-esig Sinh[ch/2], Cosh[ch/2], 0};
Hupper = {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}};
NRH`CheckZero["SM: H^{mu nu} tau^pm_nu = 0 (two-dimensional longitudinal kernel)",
   {Hupper . tauP, Hupper . tauM}];

Ysol = First@Solve[{yv1 tauP[[1]] + yv2 tauP[[2]] == 1, yv1 tauM[[1]] + yv2 tauM[[2]] == 0}, {yv1, yv2}];
YbarSol = First@Solve[{w1 tauM[[1]] + w2 tauM[[2]] == 1, w1 tauP[[1]] + w2 tauP[[2]] == 0}, {w1, w2}];
NRH`CheckZero["SM: the dual vectors Y, Ybar exist at every radius (unit clock determinant)",
   {Together[(yv1 tauP[[1]] + yv2 tauP[[2]] /. Ysol) - 1],
    Together[yv1 tauM[[1]] + yv2 tauM[[2]] /. Ysol],
    Together[(w1 tauM[[1]] + w2 tauM[[2]] /. YbarSol) - 1],
    Together[w1 tauP[[1]] + w2 tauP[[2]] /. YbarSol],
    Together[tauP[[1]] tauM[[2]] - tauP[[2]] tauM[[1]] - 1]}];

tdotd = tauP[[1]] dxp + tauP[[2]] dxm;   tdotb = tauP[[1]] bxp + tauP[[2]] bxm;
mdotd = tauM[[1]] dxp + tauM[[2]] dxm;   mdotb = tauM[[1]] bxp + tauM[[2]] bxm;
symRoute = Wc/2 (tdotd mdotb + mdotd tdotb);
antisymRoute = Wc/2 (dxp bxm - bxp dxm);
NRH`CheckZero["symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)",
   Together[symRoute - antisymRoute - Wc mdotd tdotb]];

NRH`CheckZero["at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)",
   Together[(symRoute /. ch -> 0 /. esig -> 1) - Wc/2 (dxp + 0) (bxm + 0) - Wc/2 dxm bxp]];
NRH`Check["with the 1/(2 pi alpha') prefactor this is V_W = (1/(4 pi alpha')) W dx+ dbar x-",
   Together[1/(2 Pi alphaPrime) Wc/2 - Wc/(4 Pi alphaPrime)] === 0];

gR = RiemannianMetric[Lp, Lm, u];
bR = RiemannianB[Lp, Lm, u];
et = {1/Sqrt[2], 1/Sqrt[2], 0};
ephi = {l wN/Sqrt[2], -l wN/Sqrt[2], 0};
g2 = {{et . gR . et, et . gR . ephi}, {ephi . gR . et, ephi . gR . ephi}};

ephi1 = {l/Sqrt[2], -l/Sqrt[2], 0};
gtphi = {{et . gR . et, et . gR . ephi1}, {ephi1 . gR . et, ephi1 . gR . ephi1}};
Btphi = et . bR . ephi1;
NRH`Check["before : with constant L_pm the (t, phi) components of g and B are t-independent, so d_t is Killing and B is invariant",
   FreeQ[{gtphi, Btphi}, t] && FreeQ[{gR, bR}, xp] && FreeQ[{gR, bR}, xm]];
NRH`CheckZero["before : gamma^{tau a} g_{t nu} d_a X^nu = gamma^{tau a} gamma_{a tau} = 1 on the static embedding",
   Together[Sum[Inverse[g2][[1, a]] (et . gR . {et, ephi}[[a]]), {a, 2}] - 1]];
NRH`CheckZero["before : B_{t phi} = l (e^{2y/l} + L+L- e^{-2y/l}) and B_{t nu} X'^nu = w B_{t phi}",
   {Together[Btphi - l (u + Lp Lm/u)], Together[et . bR . ephi - wN Btphi]}];
sqrtMinusGamma = l wN (u - Lp Lm/u);
Pt = 1/(2 Pi alphaPrime) (-sqrtMinusGamma Sum[Inverse[g2][[1, a]] (et . gR . {et, ephi}[[a]]), {a, 2}] + et . bR . ephi);
NRH`CheckZero["before : E = -Int_0^{2 pi} d sigma P_t with the displayed P_mu reproduces E(y) = -(2 w l/alpha') L+L- e^{-2y/l}",
   Together[-2 Pi Pt + 2 wN l/alphaPrime Lp Lm/u]];
NRH`CheckZero["-det g_(t,phi) = l^2 (e^{2y/l} - L+L- e^{-2y/l})^2 per winding, i.e. the Nambu-Goto area density l(e^{2y/l} - L+L- e^{-2y/l}) (exact, arbitrary chiral L_pm)",
   Together[-Det[g2] - (l wN (u - Lp Lm/u))^2]];

BtphiPerW = l (u + Lp Lm/u);
NRH`CheckZero["E(y) = (w l/alpha')[(e^{2y/l} - L+L- e^{-2y/l}) - (e^{2y/l} + L+L- e^{-2y/l})] = -(2 w l/alpha') L+L- e^{-2y/l}",
   Together[wN l/alphaPrime ((u - Lp Lm/u) - (u + Lp Lm/u)) + 2 wN l/alphaPrime Lp Lm/u]];
NRH`CheckZero["remark: with phi_0 = 0 the area density equals l e^{-2d} (a coincidence of the gauge choice, not a property of the Nambu-Goto action)",
   Together[l (u - Lp Lm/u) - l u (1 - Lp Lm/u^2)]];

fluxH = Integrate[D[l^2 Cos[th]^2, th], {th, 0, Pi/2}] (2 Pi) (2 Pi);
NRH`CheckZero["k = |Int_{S^3} H| / (4 pi^2 alpha') = l^2/alpha'",
   Together[(-fluxH)/(4 Pi^2 alphaPrime) - l^2/alphaPrime]];

prop[zz_] := -alphaPrime/2 Log[zz];
doubleContraction = -(1/alphaPrime) aa^2 (D[prop[z - w], z])^2;
improvement = -(aa/l) D[prop[z - w], {z, 2}];
hy = -alphaPrime/4 aa (aa + 2/l);
NRH`CheckZero["the two OPE contributions assemble to h_y(a)/(z-w)^2",
   Together[doubleContraction + improvement - hy/(z - w)^2]];
NRH`CheckZero["the marginal roots of h_y are a = 0 and a = -2/l (modes {1, e^{-2y/l}})",
   {hy /. aa -> 0, hy /. aa -> -2/l}];
NRH`CheckZero["(d_y^2 + (2/l) d_y) f = 0 for f = W0 + W1 e^{-2y/l}",
   Module[{f = W0c + W1c Exp[-2 yy/l]}, Together[D[f, {yy, 2}] + 2/l D[f, yy]]]];
NRH`CheckZero["h_y(i p_y) = (alpha'/4) p_y (p_y - 2 i/l) and P_y = p_y - i/l gives (alpha'/4)(P_y^2 + 1/l^2)",
   {Together[(-alphaPrime/4 (I pY) (I pY + 2/l)) - alphaPrime/4 pY (pY - 2 I/l)],
    Together[alphaPrime/4 pY (pY - 2 I/l) - alphaPrime/4 ((pY - I/l)^2 + 1/l^2)]}];

fieldsGO = {betaF, xpF, betabF, xmF, yF};
contractionPairs = {{betaF, xpF}, {betabF, xmF}, {yF, yF}};
vertexContent = {xpF, xmF};
NRH`Check["V_W x V_W is nonsingular: no contraction pair lies inside {x^+, x^-}^2",
   ! AnyTrue[contractionPairs, SubsetQ[vertexContent, #] &]];
NRH`Check["<x^+ x^-> = 0 in the Gomis-Ooguri system (x^+ pairs only with beta)",
   ! MemberQ[contractionPairs, {xpF, xmF}] && ! MemberQ[contractionPairs, {xmF, xpF}]];

ddprop = D[prop[z - w], z, w];
d2d2prop = D[prop[z - w], {z, 2}, {w, 2}];
cOver2 = Together[(2 (1/alphaPrime)^2 ddprop^2 + (1/l)^2 d2d2prop) (z - w)^4];
NRH`CheckZero["c_y/2 = 1/2 + 3 alpha'/l^2 from the two TT contractions, i.e. c_y = 1 + 6 alpha'/l^2",
   Together[cOver2 - (1 + 6 alphaPrime/l^2)/2]];
NRH`CheckZero["c_{beta gamma} + c_y = 2 + (1 + 6/k) = 3(k+2)/k",
   Together[2 + 1 + 6/kk - 3 (kk + 2)/kk]];

JJ6 = ODDJ[3];

xsY = {xp, xm, yy};
lamT = {lt1[xp, xm, yy], lt2[xp, xm, yy], lt3[xp, xm, yy]};
vvT = {v1[xp, xm, yy], v2[xp, xm, yy], v3[xp, xm, yy]};
xiG = Join[lamT, vvT];
dvT = Table[D[vvT[[n]], xsY[[m]]], {m, 3}, {n, 3}];
bbT = Table[D[lamT[[n]], xsY[[m]]] - D[lamT[[m]], xsY[[n]]], {m, 3}, {n, 3}];
DmatT = ArrayFlatten[{{-Transpose[dvT], 0}, {bbT, dvT}}];
NRH`CheckZero["Lhat_xi H^infty = D H^infty + H^infty D^T with D = ((-(dv)^T, 0), (b, dv)), b = d lambda~",
   Map[Together, GenLieH[xiG, Hinf, xsY] - (DmatT . Hinf + Hinf . Transpose[DmatT]), {2}]];

bS = {{0, b12, b13}, {-b12, 0, b23}, {-b13, -b23, 0}};
dvS = Table[dvs[m, n], {m, 3}, {n, 3}];
DmatS = ArrayFlatten[{{-Transpose[dvS], 0}, {bS, dvS}}];
hW = ConstantArray[0, {6, 6}]; hW[[4, 5]] = ww; hW[[5, 4]] = ww;
eqsW = DeleteCases[Union[Flatten[DmatS . Hinf + Hinf . Transpose[DmatS] - hW]], 0];
dilW = vy0 (-1/l) - 1/2 (dvs[1, 1] + dvs[2, 2] + dvs[3, 3]);
unkW = {dvs[3, 1], dvs[3, 2], dvs[3, 3], dvs[2, 1], dvs[1, 2], b12, b13, b23, vy0};
solW = Solve[Join[eqsW, {dilW}] == 0, unkW];
NRH`Check["the conditions force d_y v^mu = 0, d_- v^+ = 0 = d_+ v^-, b_{+-} = -varpi/2, b_{+y} = d_+ v^y, b_{-y} = -d_- v^y, v^y = -(l/2)(d_+ v^+ + d_- v^-)",
   Length[solW] == 1 &&
   Together[(unkW /. First[solW]) - {0, 0, 0, 0, 0, -ww/2, dvs[1, 3], -dvs[2, 3], -l/2 (dvs[1, 1] + dvs[2, 2])}] === {0, 0, 0, 0, 0, 0, 0, 0, 0}];
NRH`CheckZero["(db)_{+-y} = d_+ b_{-y} + d_- b_{y+} + d_y b_{+-} = -2 d_+ d_- v^y - (1/2) d_y varpi",
   Module[{bp = {{0, -w2[xp, xm, yy]/2, D[vyf[xp, xm, yy], xp]}, {w2[xp, xm, yy]/2, 0, -D[vyf[xp, xm, yy], xm]},
       {-D[vyf[xp, xm, yy], xp], D[vyf[xp, xm, yy], xm], 0}}},
      Together[D[bp[[2, 3]], xp] + D[bp[[3, 1]], xm] + D[bp[[1, 2]], yy]
         + 2 D[vyf[xp, xm, yy], xp, xm] + 1/2 D[w2[xp, xm, yy], yy]]]];
NRH`CheckZero["v^y = -(l/2)(d_+ v^+(x^+) + d_- v^-(x^-)) has d_+ d_- v^y = 0, so closure forces d_y varpi = 0: W_0 is gauge, e^{-2y/l} W_1 is not",
   D[-l/2 (D[vpf[xp], xp] + D[vmf[xm], xm]), xp, xm]];

NRH`CheckZero["e^{a y(z)} e^{a y(0)} ~ |z|^{-alpha' a^2} = |z|^{-4 alpha'/l^2} for a = -2/l (from <y y> = -(alpha'/2) Log|z|^2)",
   Together[-alphaPrime (-2/l)^2 + 4 alphaPrime/l^2]];
hn = nn + qq nn (1 - nn);
NRH`CheckZero["the n-fold fused weight h_n = n + q n(1-n) equals n + h_y(-2n/l) with q = alpha'/l^2",
   Together[(hn - (nn + (-alphaPrime/4 (-2 nn/l) (-2 nn/l + 2/l)))) /. qq -> alphaPrime/l^2]];
NRH`Check["h_n = 1 exactly at n = 1 or at the resonant value q = 1/n",
   Solve[hn == 1, qq] === {{qq -> 1/nn}} && Together[(hn /. nn -> 1) - 1] === 0];

NRH`FileSummary[];
