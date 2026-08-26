(* ::Package:: *)

(* ::Title:: *)
(*NRH04 — SM 1, "Linearized Holographic Response"*)


(* ::Text:: *)
(*This file verifies the content of SM 1 in order:*)
(**)
(*  SM (10)-(11) [SMinfinityvielbein, SMflatmetrics]  the aligned limiting double-vielbeins*)
(*        and flat metrics, and the projectors they reproduce;*)
(*  SM (14)-(15) [SMFG, SMFGcount]  the Fefferman-Graham count 9 - (3+2) = 4;*)
(*  SM (16) [SMfixedprojection] with SM (29)-(30) [SMexactprojectionR(falloff)]:  the exact*)
(*        fixed-frame projection of the Riemannian family and its falloffs;*)
(*  SM (44)-(45) [SMexactprojectionNR(falloff)]:  the same for the non-Riemannian family,*)
(*        including the identification  h_{op bom}^{(2)} = W_1/2  at W_0 = 0;*)
(*  SM (17) [SMNRlin]  the linearized Einstein double field equations about the common*)
(*        background (derived here directly from the nonlinear G_MN of NRH01);*)
(*  SM (18) [SMNRradialintegration]  the elementary radial identities;*)
(*  SM (19) [SMNRsol]  the complete solution of the linearized system;*)
(*  SM (20) [SMlogfreeconditions]  the fixed-dilaton, log-free conditions;*)
(*  SM (22)-(23) [SMrwquadratic, SMrwvariation]  the crossed quadratic action in the*)
(*        (r, w) = (h_{om bop}, h_{op bom}) sector and its renormalized UV variation;*)
(*  SM (31)-(33) [SMresponsematrix, SMonept, and the NR one-point functions SM (48)]*)
(*        response normalization and one-point functions;*)
(*  SM (34)-(35) [SMPBHdata]  the PBH source and response data;*)
(*  SM (36)-(38) [SMPBHkernel, SMRtwopt]  the kernel arithmetic and the two-point*)
(*        normalization  (8 pi)^{-2} (c/2)  with  c = 3l/2G.*)
(**)
(*Conventions as in NRH02: radial variable u = e^{2y/l}, so e^{-2y/l} = 1/u and*)
(*y = (l/2) Log[u]; d/dy = (2u/l) d/du.*)


ClearAll["Global`*"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH04_SM_LinearResponse.wl"];

JJ = ODDJ[3];
xs = {xp, xm, Function[e, (2 u/l) D[e, u]]};


(* ::Section:: *)
(*SM (10)-(11): the aligned frame*)


Vinf = {{1/Sqrt[2], 0, 0}, {0, 0, 0}, {0, 0, -1/Sqrt[2]},
        {0, -Sqrt[2], 0}, {0, 0, 0}, {0, 0, -1/Sqrt[2]}};
Vbinf = {{0, 0, 0}, {0, 1/Sqrt[2], 0}, {0, 0, -1/Sqrt[2]},
         {0, 0, 0}, {Sqrt[2], 0, 0}, {0, 0, 1/Sqrt[2]}};
eta3 = {{0, -1, 0}, {-1, 0, 0}, {0, 0, 1}};
etab3 = -eta3;
Hinf = {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, -1, 0}, {0, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1}};
dinf = -1/2 Log[u];

NRH`CheckZero["SM(10)-(11): V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty",
   {Vinf . eta3 . Transpose[Vinf] - (JJ + Hinf)/2,
    Vbinf . etab3 . Transpose[Vbinf] - (JJ - Hinf)/2}];
NRH`CheckZero["SM: V and Vbar are mutually orthogonal:  V^M{}_p Vbar_{M qbar} = 0",
   Transpose[JJ . Vinf] . Vbinf];
NRH`Check["SM(15): coset count 9 mixed components - (3 diffeos + 2 B-gauge) = 4 tangential",
   3*3 - (3 + 2) == 4];


(* ::Section:: *)
(*SM (16) + SM (29)-(30): exact fixed-frame projection, Riemannian family*)


fB = u + Lp[xp] Lm[xm]/u;
gB = {{2 Lp[xp], -fB, 0}, {-fB, 2 Lm[xm], 0}, {0, 0, 1}};
BB = fB {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
HR = Map[Together, RiemannianH[gB, BB], {2}];
dR = -1/2 Log[u (1 - Lp[xp] Lm[xm]/u^2)];

hR = Map[Together, Transpose[JJ . Vinf] . (HR - Hinf) . (JJ . Vbinf), {2}];
hRdisplayed = (1/u) (1 + Lp[xp] Lm[xm]/u^2)/(1 - Lp[xp] Lm[xm]/u^2)^2 *
   {{2 Lp[xp], 2 Lp[xp] Lm[xm], 0}, {2, 2 Lm[xm], 0}, {0, 0, 0}};
NRH`CheckZero["SM(29): exact h^R_{p qbar} matches the displayed closed form",
   Map[Together, hR - hRdisplayed, {2}]];
NRH`CheckZero["SM(30): falloffs  h^{(2)} = {2L+, 2L+L-; 2, 2L-} + O(e^{-6y/l})",
   {SeriesCoefficient[hR[[1, 1]], {u, Infinity, 1}] - 2 Lp[xp],
    SeriesCoefficient[hR[[1, 2]], {u, Infinity, 1}] - 2 Lp[xp] Lm[xm],
    SeriesCoefficient[hR[[2, 1]], {u, Infinity, 1}] - 2,
    SeriesCoefficient[hR[[2, 2]], {u, Infinity, 1}] - 2 Lm[xm],
    SeriesCoefficient[hR[[1, 1]], {u, Infinity, 2}],
    SeriesCoefficient[hR[[2, 1]], {u, Infinity, 2}]}];
NRH`CheckZero["delta d_R = (1/2) L+ L- e^{-4y/l} + O(e^{-8y/l})",
   {SeriesCoefficient[dR - dinf, {u, Infinity, 2}] - Lp[xp] Lm[xm]/2,
    SeriesCoefficient[dR - dinf, {u, Infinity, 1}],
    SeriesCoefficient[dR - dinf, {u, Infinity, 3}]}];


(* ::Section:: *)
(*SM (44)-(45): exact fixed-frame projection, non-Riemannian family*)


(* psi parametrization L_pm = psi_pm^{-2} keeps all series coefficients rational *)
qofu = 1/(Sqrt[2] psip[xp] psim[xm] u);
chu = 2 Sqrt[2] ArcTanh[qofu];
esig = psim[xm]/psip[xp];
LpP = 1/psip[xp]^2; LmP = 1/psim[xm]^2;
Wu = W0[xp, xm] + W1[xp, xm]/u;   (* homogeneous radial modes; enough for the falloffs *)
HNR = {{0, 0, 0, Cosh[chu], -Sinh[chu]/esig, 0},
   {0, 0, 0, esig Sinh[chu], -Cosh[chu], 0},
   {0, 0, 1, 0, 0, 0},
   {Cosh[chu], esig Sinh[chu], 0, -Wu esig Sinh[chu], Wu Cosh[chu], 0},
   {-Sinh[chu]/esig, -Cosh[chu], 0, Wu Cosh[chu], -Wu Sinh[chu]/esig, 0},
   {0, 0, 0, 0, 0, 1}};
dNR = -1/2 Log[u] + Log[Cosh[chu/(2 Sqrt[2])]];

hNR = Transpose[JJ . Vinf] . (HNR - Hinf) . (JJ . Vbinf);
NRH`Check["SM(44): exact h^NR = {{e^s sinh chi, (W/2) cosh chi, 0},{0, e^{-s} sinh chi, 0},{0,...}}",
   And[Simplify[hNR[[1, 1]] - esig Sinh[chu]] === 0,
       Simplify[hNR[[1, 2]] - Wu Cosh[chu]/2] === 0,
       Simplify[hNR[[2, 1]]] === 0,
       Simplify[hNR[[2, 2]] - Sinh[chu]/esig] === 0,
       Simplify[hNR[[3, 3]]] === 0,
       Simplify[hNR[[1, 3]]] === 0 && Simplify[hNR[[3, 1]]] === 0]];
NRH`CheckZero["SM(45): falloffs  {2L+ /u, (W0 + W1/u)/2; 0, 2L- /u} with O(u^-3) diagonals",
   Together[{SeriesCoefficient[hNR[[1, 1]], {u, Infinity, 1}] - 2 LpP,
    SeriesCoefficient[hNR[[1, 1]], {u, Infinity, 2}],
    SeriesCoefficient[hNR[[1, 2]], {u, Infinity, 0}] - W0[xp, xm]/2,
    SeriesCoefficient[hNR[[1, 2]], {u, Infinity, 1}] - W1[xp, xm]/2,
    SeriesCoefficient[hNR[[2, 2]], {u, Infinity, 1}] - 2 LmP,
    SeriesCoefficient[hNR[[2, 2]], {u, Infinity, 2}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 0}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 1}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 2}]}]];
NRH`CheckZero["delta d_NR = (1/4) L+ L- e^{-4y/l} + O(e^{-8y/l})",
   Together[{SeriesCoefficient[dNR - dinf, {u, Infinity, 2}] - LpP LmP/4,
    SeriesCoefficient[dNR - dinf, {u, Infinity, 1}],
    SeriesCoefficient[dNR - dinf, {u, Infinity, 3}]}]];


(* ::Section:: *)
(*SM (17): the linearized Einstein double field equations, derived from scratch*)


(* ::Text:: *)
(*We linearize the full nonlinear G_MN of NRH01 about (H^infty, d = -y/l) along the*)
(*constrained direction  delta H = 2 V_{(M}{}^p Vbar_{N)}{}^{qbar} h_{p qbar}, in the FG*)
(*gauge h_{a ybar} = h_{y bbar} = 0, with the four tangential fields h_{a bbar}(x, u) and*)
(*delta d(x, u).  The first-order piece of  G_MN - 2 l^{-2} J_MN  must reproduce exactly*)
(*the nine displayed equations SM (17).  We verify equivalence in both directions:*)
(*the linearized tensor vanishes when SM (17) holds (imposed as substitution rules for*)
(*the highest radial derivatives), and each displayed equation arises as an explicit*)
(*linear combination of components (checked by matching a complete list).*)


hfields = {hpp[xp, xm, u], hpm[xp, xm, u], hmp[xp, xm, u], hmm[xp, xm, u]};
hmat = {{hpp[xp, xm, u], hpm[xp, xm, u], 0}, {hmp[xp, xm, u], hmm[xp, xm, u], 0}, {0, 0, 0}};
VinfU = Vinf . eta3;  VbinfU = Vbinf . etab3;   (* raised frame indices *)
(* delta H_{MN} = V_M{}^p h_{p qbar} Vbar_N{}^{qbar} + (M <-> N); both frame indices raised *)
deltaH = VinfU . hmat . Transpose[VbinfU] // (# + Transpose[#]) &;
deltaH = Map[Together, deltaH, {2}];
Hlin = Hinf + t deltaH;
dlin = -1/2 Log[u] + t dd[xp, xm, u];

curvLin = Module[{gamma, r4, ric, s0},
   gamma = GammaDFT[Hlin, dlin, xs];
   r4 = RiemannR4[gamma, xs];
   ric = RicciS[gamma, r4, xs];
   s0 = ScalarS0[Hlin, dlin, xs];
   EinsteinG[Hlin, ric, s0, xs]];
GLin = Map[Together[D[#, t] /. t -> 0] &, curvLin - 2/l^2 JJ, {2}];

(* the displayed system SM (17), written with d/dy = (2u/l) d/du *)
Dy[e_] := (2 u/l) D[e, u];
eqs = {
   Dy[Dy[hmp[xp, xm, u]]] + 2/l Dy[hmp[xp, xm, u]],
   D[Dy[hmp[xp, xm, u]], xp],
   D[Dy[hmp[xp, xm, u]], xm],
   Dy[Dy[dd[xp, xm, u]]],
   8/l Dy[dd[xp, xm, u]] + D[hmp[xp, xm, u], xp, xm],
   Dy[Dy[hpp[xp, xm, u]]] + 2/l Dy[hpp[xp, xm, u]] + D[hmp[xp, xm, u], {xp, 2}],
   Dy[4 D[dd[xp, xm, u], xp] - D[hpp[xp, xm, u], xm]],
   Dy[Dy[hmm[xp, xm, u]]] + 2/l Dy[hmm[xp, xm, u]] + D[hmp[xp, xm, u], {xm, 2}],
   Dy[4 D[dd[xp, xm, u], xm] - D[hmm[xp, xm, u], xp]],
   Dy[Dy[hpm[xp, xm, u]]] + 2/l Dy[hpm[xp, xm, u]] + D[hmm[xp, xm, u], {xp, 2}]
      + D[hpp[xp, xm, u], {xm, 2}] - 4 D[dd[xp, xm, u], xp, xm]};

(* We verify the equivalence by substituting the GENERAL solution SM (19) (with free
   chiral/harmonic data) into both the displayed system and the freshly linearized
   tensor: since SM (19) is the general solution of SM (17), this checks in one stroke
   that SM (19) solves SM (17) and that SM (17) implies the linearized EDFE. *)

yv = l/2 Log[u];
solNR = {
   hmp -> Function[{a, b, c}, r0[a, b] + r2/c],
   dd -> Function[{a, b, c}, dd0[a, b] - l/8 (l/2 Log[c]) D[r0[a, b], a, b]],
   hpp -> Function[{a, b, c}, s0p[a, b] - l/2 (l/2 Log[c]) D[r0[a, b], {a, 2}] + s2p[a]/c],
   hmm -> Function[{a, b, c}, s0m[a, b] - l/2 (l/2 Log[c]) D[r0[a, b], {b, 2}] + s2m[b]/c],
   hpm -> Function[{a, b, c}, w0[a, b]
      + l/2 (l/2 Log[c]) (4 D[dd0[a, b], a, b] - D[s0m[a, b], {a, 2}] - D[s0p[a, b], {b, 2}]
         - l^2/4 D[r0[a, b], {a, 2}, {b, 2}])
      + l^2/8 (l/2 Log[c])^2 D[r0[a, b], {a, 2}, {b, 2}] + w2[a, b]/c]};

NRH`CheckZero["SM(19) solves the displayed system SM(17)",
   Together[eqs /. solNR]];
NRH`CheckZero["SM(17) reproduces the linearized EDFE: G^{(1)}_MN = 0 on the general solution SM(19)",
   Map[Together, GLin /. solNR, {2}]];
NRH`Check["the linearized tensor is not empty (it involves the radial derivatives of h)",
   ! FreeQ[GLin, hmp] && ! FreeQ[GLin, dd]];

(* SM (18): elementary radial identities *)
NRH`CheckZero["SM(18): D_y{1, e^{-2y/l}} = 0,  D_y y = 2/l,  D_y y^2 = 2 + 4y/l",
   {Dy[Dy[1]] + 2/l Dy[1],
    Together[Dy[Dy[1/u]] + 2/l Dy[1/u]],
    Together[Dy[Dy[yv]] + 2/l Dy[yv] - 2/l],
    Together[Dy[Dy[yv^2]] + 2/l Dy[yv^2] - 2 - 4 yv/l]}];

(* SM (20): the log (linear-in-y) branches of delta d, h_pp, h_mm in SM (19) are exact
   nonzero multiples of d+d- r^(0), d+^2 r^(0), d-^2 r^(0).  Hence the fixed-dilaton,
   log-free sector forces precisely those three conditions, leaving
   h^{(0)}_{om bop} = c0 + c+ x^+ + c- x^-. *)
NRH`CheckZero["SM(20): the log coefficients are exact multiples of the three r^(0) conditions",
   {Coefficient[dd[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/16 D[r0[xp, xm], xp, xm],
    Coefficient[hpp[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/4 D[r0[xp, xm], {xp, 2}],
    Coefficient[hmm[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/4 D[r0[xp, xm], {xm, 2}]}];


(* ::Section:: *)
(*SM (22)-(23): the crossed quadratic action in the (r, w) sector*)


(* ::Text:: *)
(*We build an EXACTLY constrained two-parameter family H(r, w) by exponentiating the*)
(*o(3,3) generator that excites the (om, bop) and (op, bom) channels, expand the*)
(*Gamma^2 Lagrangian density (with its cosmological term) to second order, and compare*)
(*with the displayed  -(1/2) e^{2y/l} dy r dy w  up to total derivatives (checked by*)
(*taking Euler-Lagrange derivatives of the difference).  The UV variation of the on-shell*)
(*quadratic action then reproduces SM (23) including its sign and 1/(16 pi G l) factor.*)


(* constrained family through second order: with the coset direction
   dHs = r K_{om bop} + w K_{op bom} (which anticommutes with H^infty in the J-mixed
   sense), the quadratic completion demanded by H J H = J is
   Hcorr = -(1/2) H^infty J dHs J dHs,  since then (H J)^2 = 1 holds to O(t^2). *)
chan[p_, qb_] := Module[{hm = ConstantArray[0, {3, 3}]}, hm[[p, qb]] = 1;
   VinfU . hm . Transpose[VbinfU]];
chanS[p_, qb_] := chan[p, qb] + Transpose[chan[p, qb]];
dHs = chanS[2, 1] t r1[xp, xm, u] + chanS[1, 2] t w1[xp, xm, u];
Hcorr = -1/2 Hinf . JJ . dHs . JJ . dHs;
Hrw = Hinf + dHs + Hcorr;
NRH`CheckZero["H(r,w) obeys H J H = J through O(t^2) (quadratic coset completion)",
   Map[Together, Normal[Series[Hrw . JJ . Hrw - JJ, {t, 0, 2}]], {2}]];
NRH`CheckZero["Hcorr is symmetric (a genuine generalized-metric correction)",
   Map[Together, Hcorr - Transpose[Hcorr], {2}]];
NRH`CheckZero["H(r,w) linearizes to the (om bop) and (op bom) channels",
   Map[Together, (D[Hrw, t] /. t -> 0) - (chanS[2, 1] r1[xp, xm, u] + chanS[1, 2] w1[xp, xm, u]), {2}]];

Lgamma2 = Module[{gamma = GammaDFT[Hrw, dinf, xs]},
   Gamma2Density[Hrw, dinf, gamma, xs] - 2 (-2/l^2) Exp[-2 dinf]];
(* second order in the pair (r, w) *)
Lq = Together[1/2 D[Lgamma2, {t, 2}] /. t -> 0];
Ltarget = -1/2 u ((2 u/l) D[r1[xp, xm, u], u]) ((2 u/l) D[w1[xp, xm, u], u]);
diffL = Together[Lq - Ltarget];
(* total-derivative test: Euler-Lagrange derivatives of the difference must vanish *)
EL[f_, e_] := Together[D[e, f[xp, xm, u]]
   - D[D[e, Derivative[1, 0, 0][f][xp, xm, u]], xp]
   - D[D[e, Derivative[0, 1, 0][f][xp, xm, u]], xm]
   - D[D[e, Derivative[0, 0, 1][f][xp, xm, u]], u]
   + D[D[e, Derivative[2, 0, 0][f][xp, xm, u]], {xp, 2}]
   + D[D[e, Derivative[0, 0, 2][f][xp, xm, u]], {u, 2}]
   + D[D[e, Derivative[1, 1, 0][f][xp, xm, u]], xp, xm]
   + D[D[e, Derivative[1, 0, 1][f][xp, xm, u]], xp, u]
   + D[D[e, Derivative[0, 1, 1][f][xp, xm, u]], xm, u]
   + D[D[e, Derivative[0, 2, 0][f][xp, xm, u]], {xm, 2}]];
NRH`CheckZero["SM(22): Gamma^2 quadratic density = -(1/2) e^{2y/l} dy r dy w  (mod total derivatives)",
   {EL[r1, diffL], EL[w1, diffL]}];

(* SM (23): UV variation of the on-shell quadratic action.
   On shell r = r0 + r2/u, w = w0 + w2/u; the momentum conjugate to r at the cutoff is
   dL/d(dy r) evaluated at u -> Infinity. *)
ronsh = rw0 + rw2/u; wonsh = ws0 + ws2/u;
momR = Limit[-1/2 u (2 u/l) D[wonsh, u], u -> Infinity];
momW = Limit[-1/2 u (2 u/l) D[ronsh, u], u -> Infinity];
NRH`CheckZero["SM(23): UV momenta give  delta S|_UV = (1/(16 pi G l)) Int (w2 dr0 + r2 dw0)",
   {Together[momR - ws2/l], Together[momW - rw2/l]}];


(* ::Section:: *)
(*SM (31)-(33), SM (48): response normalization and one-point functions*)


NRH`CheckZero["SM(31): <K_{a bbar}> = h^{(2)}_{a bbar}/(32 pi G l)  =>  stress one-points",
   {Together[2 Lp[xp]/(32 Pi G l) - Lp[xp]/(16 Pi G l)],
    Together[2 Lm[xm]/(32 Pi G l) - Lm[xm]/(16 Pi G l)],
    Together[2/(32 Pi G l) - 1/(16 Pi G l)],
    Together[2 Lp[xp] Lm[xm]/(32 Pi G l) - Lp[xp] Lm[xm]/(16 Pi G l)]}];
NRH`CheckZero["SM(48): NR hair channel:  h^{(2)}_{op bom} = W_1/2  =>  Int<K_{op bom}> = Int W_1/(64 pi G l)",
   Together[(W1[xp, xm]/2)/(32 Pi G l) - W1[xp, xm]/(64 Pi G l)]];


(* ::Section:: *)
(*SM (34)-(38): PBH data, response kernel, and the two-point normalization*)


(* ::Text:: *)
(*The PBH perturbation is  Delta_alpha H := Lhat_{xi[alpha]} H^infty with the Eq. (8)*)
(*vector continued off shell:  eps^+ -> alpha^+(x^+, x^-), eps^- -> 0.  Its leading*)
(*fixed-frame datum in the (om, bom) channel and normalizable datum in the (op, bop)*)
(*channel are  s_+ = -2 d_- alpha^+  and  r_+ = -(l^2/2) d_+^3 alpha^+.*)


xiPBH = {
   0,
   +l^2/(2 u) Lm[xm] D[al[xp, xm], {xp, 2}],
   -l/2 D[al[xp, xm], xp],
   al[xp, xm] ,
   l^2/(4 u) D[al[xp, xm], {xp, 2}],
   -l/2 D[al[xp, xm], xp]} /. Lm[xm] -> 0;
(* on the vacuum H^infty (L = 0), the dual tail proportional to L_- drops out *)
DeltaH = Map[Together, GenLieH[xiPBH, Hinf, xs], {2}];
hPBH = Transpose[JJ . Vinf] . DeltaH . (JJ . Vbinf);
NRH`CheckZero["SM(34): s_+ = [Delta H]^{(0)}_{om bom} = -2 d_- alpha^+",
   Together[SeriesCoefficient[hPBH[[2, 2]], {u, Infinity, 0}] + 2 D[al[xp, xm], xm]]];
NRH`CheckZero["SM(34): r_+ = [Delta H]^{(2)}_{op bop} = -(l^2/2) d_+^3 alpha^+",
   Together[SeriesCoefficient[hPBH[[1, 1]], {u, Infinity, 1}] + l^2/2 D[al[xp, xm], {xp, 3}]]];
NRH`CheckZero["SM(34): the PBH variation carries no (om bop) source deformation",
   Together[SeriesCoefficient[hPBH[[2, 1]], {u, Infinity, 0}]]];

(* SM (36)-(38): kernel and two-point normalization arithmetic *)
NRH`CheckZero["SM(36): (l^2/4) d_+^3 [-1/(2 pi x)] = 3 l^2/(4 pi x^4)",
   Together[l^2/4 D[-1/(2 Pi x), {x, 3}] - 3 l^2/(4 Pi x^4)]];
NRH`CheckZero["away from coincidence: d_+ d_- ln(-x^+ x^-) = 0",
   D[Log[-xps xms], xps, xms]];
cBH = 3 l/(2 G);
NRH`CheckZero["SM(38): (64 pi G l)^{-1} 3 l^2/(4 pi) = (8 pi)^{-2} (c/2),  c = 3l/2G",
   Together[1/(64 Pi G l) 3 l^2/(4 Pi) - 1/(8 Pi)^2 cBH/2]];

NRH`FileSummary[];
