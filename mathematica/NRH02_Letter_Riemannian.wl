(* ::Package:: *)

(* ::Title:: *)
(*NRH02 — Letter, "Riemannian Reference: Banados Holography in DFT"*)


(* ::Text:: *)
(*This file verifies, in the order of the Letter, the displayed equations of the section*)
(*"Riemannian Reference: Banados Holography in DFT":*)
(**)
(*  Eq. (1)  [label Rfields]        NS-NS Banados family (metric, B-field, dilaton);*)
(*  Eq. (2)  [label RDFTfields]     its packaging into the DFT pair (H_MN, e^{-2d});*)
(*  Eq. (3)  [label Rboundary]      the asymptotic falloffs;*)
(*  Eq. (4a) [label boundaryH]      the constant boundary data H^infty, d^infty = -y/l;*)
(*  Eq. (4b) [label Rboundaryframe] the aligned D=2 boundary double-vielbein representative;*)
(*  Eq. (9b) [label RDFTconservation] the conserved tensor T_AB built from K_{a bbar}, T_(0);*)
(*  Eq. (10)  [label Rcontinuity]    the boundary Ward identities;*)
(*  Eq. (5)  [label Rkilling]       the asymptotic-symmetry generator;*)
(*  Eq. (6)  [label RVirasoro]      delta_epsilon L_pm with the anomalous -l^2/4 term;*)
(*  the Brown-Henneaux central charge c = 3l/2G implied by Eqs. (5), (6), and (9a);*)
(*  Eq. (11) [label Rcorrelators]   the normalization arithmetic of the two-point function;*)
(*  the BTZ Killing horizon location and the long-string effective speed c_eff;*)
(*  Eq. (8) -> Eq. (9a) [Rrenvariation, RKdef]  the coefficient matching behind the dictionary,*)
(*                                   and SM (37) [SMinvariantstress];*)
(*  Eq. (9a) applied nonlinearly to the exact family: the full one-point matrix of Eq. (11) and*)
(*                                   SM (36) from the DFT connection of the saddle, and <T_(0)> = 0.*)
(**)
(*Everything is exact and symbolic for ARBITRARY chiral functions L_+(x^+), L_-(x^-).*)
(*We use the rational radial variable u := Exp[2y/l]  (so e^{-2y/l} = 1/u), with the*)
(*radial derivative  d/dy = (2u/l) d/du.  This keeps all expressions rational and lets*)
(*Mathematica decide every identity exactly.  LaTeX labels are the stable identifiers;*)
(*the current number mapping and manuscript SHA-256 are recorded in MANUSCRIPT_MAP.md.*)


ClearAll["Global`*"];
Get[FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH02_Letter_Riemannian.wl"];

xs = {xp, xm, Function[e, (2 u/l) D[e, u]]};
JJ = ODDJ[3];


(* ::Section:: *)
(*Eq. (1) -> Eq. (2):  the Banados family and its DFT variables*)


(* ::Text:: *)
(*ds^2 = dy^2 - 2 theta^+ theta^-  with  theta^pm = e^{y/l} dx^pm - e^{-y/l} L_mp dx^mp,*)
(*B = (e^{2y/l} + e^{-2y/l} L_+ L_-) dx^- ^ dx^+,  phi = 0.   In components (x^+, x^-, y):*)
(*g_{++} = 2 L_+, g_{--} = 2 L_-, g_{+-} = -(u + L_+ L_-/u), g_{yy} = 1,  B_{-+} = u + L_+ L_-/u.*)


fB = u + Lp[xp] Lm[xm]/u;
gB = {{2 Lp[xp], -fB, 0}, {-fB, 2 Lm[xm], 0}, {0, 0, 1}};
BB = fB {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};      (* B_{+-} = -fB,  B_{-+} = +fB *)
HR = Map[Together, RiemannianH[gB, BB], {2}];
dR = -1/2 Log[u (1 - Lp[xp] Lm[xm]/u^2)];        (* e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l}) *)

NRH`CheckZero["Eq.(2): H J H = J  (O(3,3) constraint)",
   Map[Together, HR . JJ . HR - JJ, {2}]];
NRH`CheckZero["Eq.(2): (e^{-2d})^2 = -det g  (phi_0 = 0)",
   Together[Exp[-2 dR]^2 + Det[gB]]];
NRH`CheckZero["Eq.(2): e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l})",
   Together[Exp[-2 dR] - u (1 - Lp[xp] Lm[xm]/u^2)]];


(* ::Section:: *)
(*The Einstein double field equation  G_MN = 2 l^-2 J_MN   (quoted after Eq. (15))*)


(* ::Text:: *)
(*The Letter states that the Riemannian fields of Eq. (2) solve the three-dimensional*)
(*Einstein double field equation G_MN = (2/l^2) J_MN, i.e. the projected Ricci vanishes*)
(*and the scalar curvature equals S_(0) = -4/l^2 (the cosmological term produced by the*)
(*S^3 flux upon reduction).  We verify this exactly for arbitrary chiral L_pm, and also*)
(*cross-check the scalar curvature computed two independent ways: from the closed-form*)
(*expression in (H, d) and from the contraction of the semi-covariant Riemann tensor.*)


curvR = DFTCurvature[HR, dR, xs];
NRH`CheckZero["EDFE tensor part: (P S Pbar)_MN = 0, arbitrary chiral L_pm", curvR["PSPbar"]];
NRH`CheckZero["EDFE scalar part: S_(0) = -4/l^2", Together[curvR["S0"] + 4/l^2]];
NRH`CheckZero["S_(0) closed form = S_(0) from S_ABCD contraction",
   Together[curvR["S0"] - ScalarS0FromS4[curvR["Gamma"], curvR["R4"], HR, xs]]];
NRH`CheckZero["G_MN = 2 l^-2 J_MN", Map[Together, curvR["G"] - 2/l^2 JJ, {2}]];

(* Connection sanity on this background: defining properties of Gamma. *)
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


(* ::Section:: *)
(*Eq. (3) - Eq. (4b):  asymptotics, type-(1,1) boundary data, and aligned frame*)


Hinf = {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, -1, 0}, {0, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1}};

NRH`CheckZero["Eq.(4a) [boundaryH]: H^infty matches the displayed 6x6 matrix",
   Map[Limit[#, u -> Infinity] &, HR, {2}] - Hinf];
NRH`CheckZero["Eq.(3): H - H^infty = O(e^{-2y/l})  (every entry vanishes at the boundary)",
   Map[Limit[#, u -> Infinity] &, HR - Hinf, {2}]];
NRH`CheckZero["Eq.(3): d + y/l = O(e^{-4y/l}):  u * (d + y/l) still vanishes at the boundary",
   Limit[u (dR + 1/2 Log[u]), u -> Infinity]];

(* The boundary block (rows/cols x~+, x~-, x+, x-) is the type-(1,1) generalized metric:
   its upper-left 2x2 block (the would-be inverse metric) vanishes identically. *)
H0 = Hinf[[{1, 2, 4, 5}, {1, 2, 4, 5}]];
NRH`Check["Eq.(4a) [boundaryH]: induced boundary H^(0) has vanishing upper-left block (type (1,1))",
   H0[[1 ;; 2, 1 ;; 2]] === {{0, 0}, {0, 0}}];


(* ::Section:: *)
(*Eq. (9b) - Eq. (10):  the conserved tensor and the boundary Ward identities*)


(* ::Text:: *)
(*On the constant boundary representative Eq. (4a), the conservation of*)
(*T_AB = 4 V^(0)_[A^a Vbar^(0)_B]^bbar K_{a bbar} - (1/2) J_AB T_(0)  reduces to ordinary*)
(*divergences.  With the lower-index boundary vielbeins of Eq. (4b) [Rboundaryframe] we*)
(*evaluate  partial^A T_AB = J^{AC} partial_C T_AB  componentwise for arbitrary response*)
(*functions K_{a bbar}(x^+, x^-), T_(0)(x^+, x^-), and read off the displayed identities:*)
(*   partial_- K_{op bop} + (1/4) partial_+ T_(0) = 0,*)
(*   partial_+ K_{om bom} + (1/4) partial_- T_(0) = 0,*)
(*   K_{om bop} = const,   and no condition on K_{op bom}.*)


J4 = ODDJ[2];
V0 = {{1/Sqrt[2], 0}, {0, 0}, {0, -Sqrt[2]}, {0, 0}};        (* Eq. (4b), rows (x~+,x~-,x+,x-) *)
Vb0 = {{0, 0}, {0, 1/Sqrt[2]}, {0, 0}, {Sqrt[2], 0}};
eta2 = {{0, -1}, {-1, 0}};   etab2 = {{0, 1}, {1, 0}};       (* flat metrics stated below Eq. (4b) *)

NRH`CheckZero["Eq.(4b) [Rboundaryframe]: boundary vielbeins reproduce P^(0) and Pbar^(0)",
   {V0 . eta2 . Transpose[V0] - (J4 + H0)/2, Vb0 . etab2 . Transpose[Vb0] - (J4 - H0)/2}];

kmat = {{Kpp[xp, xm], Kpm[xp, xm]}, {Kmp[xp, xm], Kmm[xp, xm]}};  (* K_{a bbar}, rows (op,om), cols (bop,bom) *)
T0f = T0[xp, xm];
(* K_{a bbar} carries lower frame indices, so it pairs with the raised vielbeins
   V_A{}^a = V_{Ab} eta^{ba}: *)
V0u = V0 . eta2;  Vb0u = Vb0 . etab2;
TAB = Table[
   2 Sum[(V0u[[a, i]] Vb0u[[b, j]] - V0u[[b, i]] Vb0u[[a, j]]) kmat[[i, j]], {i, 2}, {j, 2}]
   - 1/2 J4[[a, b]] T0f, {a, 4}, {b, 4}];
bD[e_, a_] := If[a <= 2, 0, D[e, {xp, xm}[[a - 2]]]];
divT = Table[Together[Sum[J4[[a, c]] bD[TAB[[a, b]], c], {a, 4}, {c, 4}]], {b, 4}];
ward1 = D[Kpp[xp, xm], xm] + 1/4 D[T0f, xp];
ward2 = D[Kmm[xp, xm], xp] + 1/4 D[T0f, xm];
(* The four components of partial^A T_AB come out as
   { partial_- K_{om bop},  -partial_+ K_{om bop},  -2*ward1,  -2*ward2 } :
   setting them to zero is exactly Eq. (10) plus the constancy of K_{om bop}. *)
NRH`CheckZero["Eq.(10): div T = {d_- K_mp, -d_+ K_mp, -2 Ward_+, -2 Ward_-} exactly",
   Together[divT - {D[Kmp[xp, xm], xm], -D[Kmp[xp, xm], xp], -2 ward1, -2 ward2}]];
NRH`Check["Eq.(10): no local condition on K_{op bom}", FreeQ[divT, Kpm]];


(* ::Section:: *)
(*Eq. (5) - Eq. (6):  asymptotic symmetry and the Virasoro transformation law*)


(* ::Text:: *)
(*The doubled vector xi of Eq. (5) has upper components  xi^M = (xitilde_mu ; xi^mu).*)
(*We verify, exactly in series around the boundary:*)
(*  (i)   Lhat_xi d = O(e^{-4y/l});*)
(*  (ii)  Lhat_xi H - delta_eps H = O(e^{-4y/l}), where delta_eps H is the variation of the*)
(*        exact family under  L_pm -> L_pm + t delta_eps L_pm  with the displayed law (6),*)
(*        delta_eps L_pm = eps partial L + 2 L partial eps - (l^2/4) partial^3 eps;*)
(*  (iii) delta_eps H itself is O(e^{-2y/l}).*)


xiUp = {  (* (xitilde_+, xitilde_-, xitilde_y ; xi^+, xi^-, xi^y) *)
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
(* HR contains L_pm only algebraically (no derivatives), so the family variation is a
   plain chain rule in the two values Lp[xp], Lm[xm]. *)
HRgen = HR /. {Lp[xp] -> LPv, Lm[xm] -> LMv};
depsH = Map[Together,
   (D[HRgen, LPv] dLp + D[HRgen, LMv] dLm) /. {LPv -> Lp[xp], LMv -> Lm[xm]}, {2}];

orderCheck[m_, pow_] := Map[Function[e, Normal[Series[e, {u, Infinity, pow}]]], m, {2}];
NRH`CheckZero["Eq.(5): Lhat_xi d = O(e^{-4y/l})",
   Normal[Series[lieD, {u, Infinity, 1}]]];
NRH`CheckZero["Eq.(5)-(6): Lhat_xi H - delta_eps H = O(e^{-4y/l})",
   orderCheck[lieH - depsH, 1]];
NRH`CheckZero["delta_eps H = O(e^{-2y/l})   (leading falloff of the Banados variation)",
   orderCheck[depsH, 0]];


(* ::Section:: *)
(*Central charge  c = 3l/2G  and the two-point normalization of Eq. (11)*)


(* ::Text:: *)
(*The Letter identifies 8 pi K_{op bop} as the holomorphic stress component with*)
(*<K_{op bop}> = L_+/(16 pi G l)  [SM (36)], so  T_hol := 8 pi <K> = L_+/(2 G l).*)
(*The anomalous term of Eq. (6) then fixes the central charge through the CFT law*)
(*   delta_eps T = eps T' + 2 T eps' - (c/12) eps''' :*)
(*   (c/12) = (l^2/4) * (1/(2 G l))  =>  c = 3l/2G.*)
(*The connected two-point normalization quoted in Eq. (11) is the arithmetic identity*)
(*   (1/(64 pi G l)) * (3 l^2/(4 pi)) = (8 pi)^{-2} (c/2),  cf. SM (39)-(40).*)


THol = 8 Pi (Lp[xp]/(16 Pi G l));
dT = (dLp/(2 G l)) /. {};   (* transformation of T_hol induced by Eq. (6) *)
cval = 3 l/(2 G);
NRH`CheckZero["c = 3l/2G reproduces delta_eps T = eps T' + 2 T eps' - (c/12) eps'''",
   Together[dT - (ep[xp] D[THol, xp] + 2 THol D[ep[xp], xp] - cval/12 D[ep[xp], {xp, 3}])]];
NRH`CheckZero["Eq.(11): (64 pi G l)^{-1} (3 l^2/(4 pi)) = (8 pi)^{-2} (c/2)",
   Together[1/(64 Pi G l) 3 l^2/(4 Pi) - 1/(8 Pi)^2 cval/2]];


(* ::Section:: *)
(*BTZ horizon and the long-string effective speed*)


NRH`CheckZero["horizon: e^{-2d} = 0 at e^{4y/l} = L+ L-  (u^2 = L+ L-)",
   Together[(Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]) ]];
EE2 = Map[Together, gB - BB, {2}];   (* sigma-model tensor E = g - B *)
NRH`Check["E_{+-} = 0, E_{-+} = -2F, E_{pm pm} = 2 L_pm  (Gomis-Ooguri channel structure)",
   And[Together[EE2[[1, 2]]] === 0,
       Together[EE2[[2, 1]] + 2 fB] === 0,
       Together[EE2[[1, 1]] - 2 Lp[xp]] === 0,
       Together[EE2[[2, 2]] - 2 Lm[xm]] === 0]];
NRH`CheckZero["c_eff^2 := 2F = 2(e^{2y/l} + L+ L- e^{-2y/l}) and = 4 Sqrt[L+ L-] at the horizon",
   Together[(2 fB /. u -> Sqrt[Lp[xp] Lm[xm]]) - 4 Sqrt[Lp[xp] Lm[xm]]]];


(* ::Section:: *)
(*Eq. (8) -> Eq. (9a) [Rrenvariation, RKdef], and SM (37) [SMinvariantstress]*)


(* ::Text:: *)
(*Eq. (8) pairs the cutoff data as  (16 pi G)^{-1} e^{-2d} [h^{p qbar} A^y_{p qbar} + 2 delta d (B^y + 4/l)],*)
(*while the DFT response definition SM (7)-(8) pairs them as  2 e^{-2d} (-h_{p qbar} K^{p qbar} + delta d T_(0)).*)
(*With e^{-2d(Y)} ~ e^{2Y/l} e^{-2d^(0)}, matching the two pairings gives exactly the coefficients*)
(*displayed in Eq. (9a):  K = -(32 pi G)^{-1} e^{2Y/l} A^y  and  T_(0) = (16 pi G)^{-1} e^{2Y/l} (B^y + 4/l).*)
(*Eq. (9b) then represents the invariant coordinate components in the aligned frame of Eq. (4b);*)
(*SM (37) reads off  T_{x^+ x~_-} = 2 K_{op bop}  and  T_{x^- x~_+} = 2 K_{om bom}  from the same*)
(*tensor T_AB assembled above (rows/columns ordered x~_+, x~_-, x^+, x^-).*)


NRH`CheckZero["Eq.(9a) [RKdef]: coefficient matching  -2K = (16 pi G)^{-1} A^y  and  2T_(0) = (16 pi G)^{-1} 2(B^y + 4/l)",
   {Together[-2 (-(1/(32 Pi G))) - 1/(16 Pi G)], Together[2 (1/(16 Pi G)) - 2/(16 Pi G)]}];
NRH`CheckZero["SM(37) [SMinvariantstress]: T_{x+ x~-} = 2 K_{op bop} and T_{x- x~+} = 2 K_{om bom} in the frame of Eq. (4b)",
   {Together[TAB[[3, 2]] - 2 Kpp[xp, xm]], Together[TAB[[4, 1]] - 2 Kmm[xp, xm]]}];
NRH`Check["SM(37): J_{x+ x~-} = 0 = J_{x- x~+}, so no T_(0) term enters these two components",
   J4[[3, 2]] === 0 && J4[[4, 1]] === 0];


(* ::Section:: *)
(*Eq. (9a) applied nonlinearly to the exact family (1)-(2): the one-point functions of Eq. (11)*)


(* ::Text:: *)
(*Eq. (9a) is a nonlinear statement:  <K_{a bbar}> = -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar}(Y),  with A^K the*)
(*generalized-metric momentum of Eq. (8)/SM (12) built from the full DFT connection of the saddle, and*)
(*the frame indices referring to the aligned boundary frame of Eq. (4b)/SM (14).  Here A^y_{MN} is*)
(*evaluated on the exact Banados family with arbitrary chiral L_pm, projected on that fixed frame, and*)
(*the limit is taken: the result is the complete one-point matrix  {{L_+, L_+ L_-}, {1, L_-}}/(16 pi G l)*)
(*- Eq. (11), SM (36), and the two mixed entries quoted after SM (37) - obtained without any*)
(*linearization.  In the same way  e^{2Y/l}(B^y + 4/l) -> 0  gives  <T_(0)> = 0,  as stated below Eq. (10).*)


VinfL = {{1/Sqrt[2], 0, 0}, {0, 0, 0}, {0, 0, 1/Sqrt[2]},
   {0, -Sqrt[2], 0}, {0, 0, 0}, {0, 0, 1/Sqrt[2]}};       (* SM (14), lower local indices *)
VbinfL = {{0, 0, 0}, {0, 1/Sqrt[2], 0}, {0, 0, 1/Sqrt[2]},
   {0, 0, 0}, {Sqrt[2], 0, 0}, {0, 0, -1/Sqrt[2]}};
AyR = MomentumAK[HR, curvR["Gamma"], xs][[6]];
AfixR = Map[Together, Transpose[JJ . VinfL] . AyR . (JJ . VbinfL), {2}];
KR = Map[Limit[-(1/(32 Pi G)) u #, u -> Infinity] &, AfixR[[1 ;; 2, 1 ;; 2]], {2}];
NRH`CheckZero["Eq.(9a) on the exact family: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l)  [Eq.(11), SM(36), mixed entries after SM(37)]",
   Map[Together, KR - {{Lp[xp], Lp[xp] Lm[xm]}, {1, Lm[xm]}}/(16 Pi G l), {2}]];
NRH`CheckZero["Eq.(9a) on the exact family: e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0",
   Limit[u (GammaBVector[HR, dR, xs][[6]] + 4/l), u -> Infinity]];

NRH`FileSummary[];
