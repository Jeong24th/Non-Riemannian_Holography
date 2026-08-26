(* ::Package:: *)

(* ::Title:: *)
(*NRH06 — SM 5, "Worldsheet Realizations"*)


(* ::Text:: *)
(*This file verifies the worldsheet statements of SM 5 in order:*)
(**)
(*  SM (75) [SMRfirstorder]  the exact first-order rewriting of the Riemannian sigma model:*)
(*        integrating out (beta, betabar) reproduces  E_{mu nu} d x^mu dbar x^nu  with*)
(*        E = g - B, i.e. E_{+-} = 0, E_{-+} = -2F, E_{pm pm} = 2 L_pm;*)
(*  SM (76) [SMceff]  c_eff^2 = 2F, its horizon value, and  -det g_par = F^2 - 4 L+ L- = e^{-4d};*)
(*  SM (78)-(79) [SMdygconstraints, SMdygGO]  the non-Riemannian reduction: the clock*)
(*        kernel of H^{mu nu}, the dual vectors Y, Ybar, and the exact reduced Lagrangian;*)
(*        the equality of the antisymmetric-clock route SM (94) [SMdeltaL] with the*)
(*        symmetric block descent, modulo the product of the two constraints;*)
(*  SM (80) [SMGO], SM (81) [SMvertex]  the Gomis-Ooguri limit and the hair vertex*)
(*        coefficient  (1/4 pi alpha') W d x^+ dbar x^-;*)
(*  SM (82) [SMlongstringE]  the winding-string potential, including the exact identity*)
(*        Sqrt[-det g_2] = l e^{-2d} for the static winding probe;*)
(*  SM (95) [SMWZWlevel]  the flux level  k = l^2/alpha';*)
(*  SM (96)-(98) [SMFTweight, SMBRSTradial, SMBRSTmomentum]  the Fradkin-Tseytlin-improved*)
(*        radial weight  h_y(a) = -(alpha'/4) a (a + 2/l)  and its two marginal roots*)
(*        {1, e^{-2y/l}} - the two homogeneous hair modes W_0, W_1;*)
(*  the Gomis-Ooguri contraction bookkeeping behind the exact marginality of V_W[W_0].*)


ClearAll["Global`*"];
Get[FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH06_SM_Worldsheet.wl"];


(* ::Section:: *)
(*SM (75)-(76): the Riemannian first-order form and c_eff*)


(* worldsheet bilinears are represented by commuting symbols dxp = d x^+, bxp = dbar x^+ etc. *)
FF = u + Lp Lm/u;
L1 = dy by + 2 Lp dxp bxp + 2 Lm dxm bxm + beta bxp + betab dxm + beta betab/(2 FF);
betaSol = First@Solve[{D[L1, beta] == 0, D[L1, betab] == 0}, {beta, betab}];
NRH`Check["SM(75): the auxiliary equations give beta = -2F d x^- and betabar = -2F dbar x^+",
   Together[(beta /. betaSol) + 2 FF dxm] === 0 && Together[(betab /. betaSol) + 2 FF bxp] === 0];
L2 = Together[L1 /. betaSol];
EmatR = {{2 Lp, 0, 0}, {-2 FF, 2 Lm, 0}, {0, 0, 1}};   (* E_{mu nu}, rows mu = (+,-,y): E_{-+} = -2F *)
LE = Sum[EmatR[[m, n]] {dxp, dxm, dy}[[m]] {bxp, bxm, by}[[n]], {m, 3}, {n, 3}];
NRH`CheckZero["SM(75): eliminating the auxiliaries reproduces E_{mu nu} dx^mu dbar x^nu, E = g - B",
   Together[L2 - LE]];
NRH`CheckZero["SM(76): -det g_par = F^2 - 4 L+ L- = e^{-4 d_R}",
   Together[-Det[{{2 Lp, -FF}, {-FF, 2 Lm}}] - (FF^2 - 4 Lp Lm)]];
NRH`CheckZero["SM(76): c_eff^2 = 2F -> 4 Sqrt[L+L-] at the horizon u = Sqrt[L+L-]",
   Together[(2 FF /. u -> Sqrt[Lp Lm]) - 4 Sqrt[Lp Lm]]];


(* ::Section:: *)
(*SM (78)-(79): the non-Riemannian clock kernel and the reduced Lagrangian*)


(* clock forms and dual vectors at radius chi *)
tauP = {Cosh[ch/2], -es^-1 Sinh[ch/2], 0};   (* es = e^sigma *)
tauM = {-es Sinh[ch/2], Cosh[ch/2], 0};
Hupper = {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}};   (* upper-left block of the NR generalized metric *)
NRH`CheckZero["SM: H^{mu nu} tau^pm_nu = 0  (two-dimensional longitudinal kernel)",
   {Hupper . tauP, Hupper . tauM}];
(* dual vectors: Y tau+ = 1, Y tau- = 0; Ybar tau- = 1, Ybar tau+ = 0 *)
Ysol = First@Solve[{yv1 tauP[[1]] + yv2 tauP[[2]] == 1, yv1 tauM[[1]] + yv2 tauM[[2]] == 0}, {yv1, yv2}];
YbarSol = First@Solve[{w1 tauM[[1]] + w2 tauM[[2]] == 1, w1 tauP[[1]] + w2 tauP[[2]] == 0}, {w1, w2}];
NRH`CheckZero["SM: the dual vectors Y, Ybar exist at every radius (unit clock determinant)",
   {Together[(yv1 tauP[[1]] + yv2 tauP[[2]] /. Ysol) - 1],
    Together[yv1 tauM[[1]] + yv2 tauM[[2]] /. Ysol],
    Together[(w1 tauM[[1]] + w2 tauM[[2]] /. YbarSol) - 1],
    Together[w1 tauP[[1]] + w2 tauP[[2]] /. YbarSol],
    Together[tauP[[1]] tauM[[2]] - tauP[[2]] tauM[[1]] - 1]}];

(* two descent routes for the W coupling:
   symmetric block:  (W/2)(tau+ . dx)(tau- . dbar x) + (W/2)(tau- . dx)(tau+ . dbar x)
   antisymmetric clock (SM (94)):  (W/2)(dx+ dbar x- - dbar x+ dx-)
   difference = W (tau- . dx)(tau+ . dbar x): the product of the two constraints. *)
tdotd = tauP[[1]] dxp + tauP[[2]] dxm;   tdotb = tauP[[1]] bxp + tauP[[2]] bxm;
mdotd = tauM[[1]] dxp + tauM[[2]] dxm;   mdotb = tauM[[1]] bxp + tauM[[2]] bxm;
symRoute = Wc/2 (tdotd mdotb + mdotd tdotb);
antisymRoute = Wc/2 (dxp bxm - bxp dxm);
NRH`CheckZero["SM(94): symmetric-block route - antisymmetric-clock route = W (tau- . dx)(tau+ . dbar x)",
   Together[symRoute - antisymRoute - Wc mdotd tdotb]];
(* at the boundary chi -> 0 the clocks become dx^pm and the coupling is the GO vertex *)
NRH`CheckZero["SM(80)-(81): at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^-  (+ constraint terms)",
   Together[(symRoute /. ch -> 0 /. es -> 1) - Wc/2 (dxp + 0) (bxm + 0) - Wc/2 dxm bxp]];
NRH`Check["SM(81): with the 1/(2 pi alpha') prefactor this is V_W = (1/(4 pi alpha')) W dx+ dbar x-",
   Together[1/(2 Pi ap) Wc/2 - Wc/(4 Pi ap)] === 0];


(* ::Section:: *)
(*SM (82): the static winding probe*)


(* embedding: t = tau, phi = w sigma; x^pm = (t pm l phi)/Sqrt[2]; fixed radius y *)
gBan = {{2 Lp, -(u + Lp Lm/u), 0}, {-(u + Lp Lm/u), 2 Lm, 0}, {0, 0, 1}};
et = {1/Sqrt[2], 1/Sqrt[2], 0};          (* d x^mu / dt *)
ephi = {l wN/Sqrt[2], -l wN/Sqrt[2], 0}; (* d x^mu / d sigma, phi = w sigma *)
g2 = {{et . gBan . et, et . gBan . ephi}, {ephi . gBan . et, ephi . gBan . ephi}};
NRH`CheckZero["SM(82): Sqrt[-det g_2] = l w e^{-2d}   (exact, arbitrary chiral L_pm)",
   Together[-Det[g2] - (l wN (u - Lp Lm/u))^2]];
(* B_{t phi} = l (e^{2y/l} + L+ L- e^{-2y/l}) per winding unit *)
BtphiPerW = l (u + Lp Lm/u);
NRH`CheckZero["SM(82): E(y) = (w l/alpha')[e^{-2d} - (e^{2y/l} + L+L- e^{-2y/l})] = -(2 w l/alpha') L+L- e^{-2y/l}",
   Together[wN l/ap ((u - Lp Lm/u) - (u + Lp Lm/u)) + 2 wN l/ap Lp Lm/u]];


(* ::Section:: *)
(*SM (95): the flux level*)


(* B_{phi1 phi2} = l^2 cos^2 theta on S^3 (theta in [0, Pi/2]):  |Int H| = 4 pi^2 l^2 *)
fluxH = Integrate[D[l^2 Cos[th]^2, th], {th, 0, Pi/2}] (2 Pi) (2 Pi);   (* = -4 pi^2 l^2 *)
NRH`CheckZero["SM(95): k = |Int_{S^3} H| / (4 pi^2 alpha') = l^2/alpha'",
   Together[(-fluxH)/(4 Pi^2 ap) - l^2/ap]];


(* ::Section:: *)
(*SM (96)-(98): the Fradkin-Tseytlin radial weight*)


(* ::Text:: *)
(*Free radial OPE:  y(z) y(w) ~ -(alpha'/2) Log(z - w).  The improved stress tensor is*)
(*T_y = -(1/alpha') :(dy)^2: - (1/l) d^2 y.  Acting on e^{a y(w)}:*)
(*  - the double contraction of :(dy)^2: gives  -(1/alpha') a^2 [d_z<y y>]^2 = -(alpha'/4) a^2/(z-w)^2;*)
(*  - the improvement gives  -(a/l) d_z^2 <y(z) y(w)> = -(alpha' a)/(2 l) /(z-w)^2.*)
(*Hence  h_y(a) = -(alpha'/4) a (a + 2/l),  whose marginal roots  h_y = 0  are  a = 0 and*)
(*a = -2/l:  exactly the two homogeneous hair modes  {1, e^{-2y/l}} = {W_0, W_1} dressing.*)


prop[zz_] := -ap/2 Log[zz];
doubleContraction = -(1/ap) aa^2 (D[prop[z - w], z])^2;
improvement = -(aa/l) D[prop[z - w], {z, 2}];
hy = -ap/4 aa (aa + 2/l);
NRH`CheckZero["SM(96): the two OPE contributions assemble to  h_y(a)/(z-w)^2",
   Together[doubleContraction + improvement - hy/(z - w)^2]];
NRH`CheckZero["SM(97): the marginal roots of h_y are a = 0 and a = -2/l  (modes {1, e^{-2y/l}})",
   {hy /. aa -> 0, hy /. aa -> -2/l}];
NRH`CheckZero["SM(97): (d_y^2 + (2/l) d_y) f = 0  for  f = W0 + W1 e^{-2y/l}",
   Module[{f = W0c + W1c Exp[-2 yy/l]}, Together[D[f, {yy, 2}] + 2/l D[f, yy]]]];
NRH`CheckZero["SM(98): h_y(i k) = (alpha'/4) k (k - 2 i/l) and P = k - i/l gives (alpha'/4)(P^2 + 1/l^2)",
   {Together[(-ap/4 (I k) (I k + 2/l)) - ap/4 k (k - 2 I/l)],
    Together[ap/4 k (k - 2 I/l) - ap/4 ((k - I/l)^2 + 1/l^2)]}];


(* ::Section:: *)
(*Gomis-Ooguri contraction bookkeeping and exact marginality of V_W[W_0]*)


(* ::Text:: *)
(*In the Gomis-Ooguri system SM (80) the only singular contractions are*)
(*beta(z) x^+(w) ~ 1/(z-w) and betabar(zbar) x^-(wbar) ~ 1/(zbar-wbar); in particular*)
(*<x^+ x^-> = 0.  The vertex V_W = (1/(4 pi alpha')) W(x^+, x^-) dx^+ dbar x^- contains*)
(*neither beta nor betabar, so V_W x V_W has no singular contraction: the y-independent*)
(*deformation W_0 is exactly marginal.  We encode the bookkeeping as a structural check*)
(*on the contraction table.*)


fieldsGO = {betaF, xpF, betabF, xmF, yF};
contractionPairs = {{betaF, xpF}, {betabF, xmF}, {yF, yF}};   (* the ONLY singular pairs *)
vertexContent = {xpF, xmF};   (* fields appearing in V_W (x^pm and their derivatives) *)
NRH`Check["V_W x V_W is nonsingular: no contraction pair lies inside {x^+, x^-}^2",
   ! AnyTrue[contractionPairs, SubsetQ[vertexContent, #] &]];
NRH`Check["<x^+ x^-> = 0 in the Gomis-Ooguri system  (x^+ pairs only with beta)",
   ! MemberQ[contractionPairs, {xpF, xmF}] && ! MemberQ[contractionPairs, {xmF, xpF}]];

NRH`FileSummary[];
