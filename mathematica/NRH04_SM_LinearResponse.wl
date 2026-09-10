(* ::Title:: *)
(*NRH04 SM LinearResponse*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH04_SM_LinearResponse.wl"];

JJ = ODDJ[3];
xs = {xp, xm, Function[e, (2 u/l) D[e, u]]};

dinf = -1/2 Log[u];

NRH`CheckZero["V eta V^T = P^infty and Vbar etabar Vbar^T = Pbar^infty",
   {Vinf . eta3 . Transpose[Vinf] - (JJ + Hinf)/2,
    Vbinf . etab3 . Transpose[Vbinf] - (JJ - Hinf)/2}];
NRH`CheckZero["SM: V and Vbar are mutually orthogonal: V^M{}_p Vbar_{M qbar} = 0",
   Transpose[JJ . Vinf] . Vbinf];
NRH`Check["coset count 9 mixed components - (3 diffeos + 2 B-gauge) = 4 tangential",
   3*3 - (3 + 2) == 4];

hgen = Table[hh[p, q], {p, 3}, {q, 3}];
dVlow = 1/2 (Vbinf . etab3) . Transpose[hgen];
dVblow = -1/2 (Vinf . eta3) . hgen;
dHgen = (Vinf . eta3) . hgen . Transpose[Vbinf . etab3] // (# + Transpose[#]) &;
NRH`CheckZero["delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2 for generic h_{p qbar}",
   {dVlow . eta3 . Transpose[Vinf] + Vinf . eta3 . Transpose[dVlow] - dHgen/2,
    dVblow . etab3 . Transpose[Vbinf] + Vbinf . etab3 . Transpose[dVblow] + dHgen/2}];

w2u = Sqrt[u];
E3btz = {{w2u, -Lm0/w2u, 0}, {-Lp0/w2u, w2u, 0}, {0, 0, 1}};
g3btz = Map[Together, Transpose[E3btz] . eta3 . E3btz, {2}];
B3btz = (u + Lp0 Lm0/u) {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
Eibtz = Map[Together, Inverse[E3btz], {2}];
Vbtz = Map[Together, 1/Sqrt[2] ArrayFlatten[{{Eibtz}, {(g3btz + B3btz) . Eibtz}}], {2}];
Vbbtz = Map[Together, 1/Sqrt[2] ArrayFlatten[{{Eibtz}, {(B3btz - g3btz) . Eibtz}}], {2}];
Hbtz = Map[Together, RiemannianH[g3btz, B3btz], {2}];
dbtz = -1/2 Log[u (1 - Lp0 Lm0/u^2)];
NRH`CheckZero["setup: the BTZ Riemannian double vielbein obeys V eta V^T = P, Vbar etabar Vbar^T = Pbar, V^T J Vbar = 0",
   {Map[Together, Vbtz . eta3 . Transpose[Vbtz] - (JJ + Hbtz)/2, {2}],
    Map[Together, Vbbtz . etab3 . Transpose[Vbbtz] - (JJ - Hbtz)/2, {2}],
    Map[Together, Transpose[JJ . Vbtz] . Vbbtz, {2}]}];
hmatV = {{hpp[xp, xm, u], hpm[xp, xm, u], 0}, {hmp[xp, xm, u], hmm[xp, xm, u], 0}, {0, 0, 0}};
deltaHV = Map[Together, (Vbtz . eta3) . hmatV . Transpose[Vbbtz . etab3] // (# + Transpose[#]) &, {2}];
ddV = dd[xp, xm, u];
gammaTV = GammaDFT[Hbtz + t deltaHV, dbtz + t ddV, xs];
gamma0V = Map[Together, gammaTV /. t -> 0, {3}];
gamma1V = Map[Together, D[gammaTV, t] /. t -> 0, {3}];
LHS11 = Together[D[Gamma2Density[Hbtz + t deltaHV, dbtz + t ddV, gamma0V + t gamma1V, xs], t] /. t -> 0];
r4V = RiemannR4[gamma0V, xs]; ricV = RicciS[gamma0V, r4V, xs]; s0V = Together[ScalarS0[Hbtz, dbtz, xs]];
VupV = JJ . Vbtz; VbupV = JJ . Vbbtz;
SpqV = Map[Together, Transpose[VupV] . ricV . VbupV, {2}];
hupV = eta3 . hmatV . etab3;
bulk11 = 2 Exp[-2 dbtz] (Sum[hupV[[p, q]] SpqV[[p, q]], {p, 3}, {q, 3}] - ddV s0V);
trvecV = Table[Sum[JJ[[L, A]] gamma0V[[A, L, N]], {L, 6}, {A, 6}], {N, 6}];
NRH`CheckZero["setup: Gamma^L_{LN} = -2 d_N d on the BTZ background (the trace entering A^K)",
   Together[trvecV + 2 DblGrad[dbtz, xs]]];
trpV = Table[Sum[trvecV[[M]] VupV[[M, p]], {M, 6}], {p, 3}];
trqV = Table[Sum[trvecV[[N]] VbupV[[N, q]], {N, 6}], {q, 3}];
GpqV = Table[Sum[VupV[[M, p]] JJ[[K, A]] gamma0V[[M, A, N]] VbupV[[N, q]], {M, 6}, {A, 6}, {N, 6}], {K, 6}, {p, 3}, {q, 3}];
GqpV = Table[Sum[VbupV[[N, q]] JJ[[K, A]] gamma0V[[N, A, M]] VupV[[M, p]], {M, 6}, {A, 6}, {N, 6}], {K, 6}, {p, 3}, {q, 3}];
AKV = Table[VupV[[K, p]] trqV[[q]] + VbupV[[K, q]] trpV[[p]] - GpqV[[K, p, q]] - GqpV[[K, p, q]], {K, 6}, {p, 3}, {q, 3}];
BKV = GammaBVector[Hbtz, dbtz, xs];
flux11 = Sum[DblD[Exp[-2 dbtz] Sum[hupV[[p, q]] AKV[[K, p, q]], {p, 3}, {q, 3}] + 2 ddV Exp[-2 dbtz] BKV[[K]], K, xs], {K, 6}];
NRH`CheckZero["delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) + d_K(h e^{-2d} A^K + 2 delta d e^{-2d} B^K) exactly, generic tangential h and delta d on BTZ",
   Together[LHS11 - bulk11 - flux11]];
NRH`Check["the identity is not vacuous - the flux term carries the fluctuations and S_(0) = -4/l^2 on this on-shell background",
   ! FreeQ[flux11, hpp] && ! FreeQ[flux11, dd] && Together[s0V + 4/l^2] === 0 && NRH`ZeroQ[SpqV]];

NRH`CheckZero["on shell (S_{p qbar} = 0, S_(0) = 2 Lambda = -4/l^2) the bulk term 2e^{-2d}(-delta d S_(0)) cancels delta(-2 Lambda e^{-2d})",
   Together[2 (-dd0) (-4/l^2) + 4 (-2/l^2) dd0]];
NRH`CheckZero["delta S_ct = -(16 pi G)^{-1}(4/l) delta e^{-2d} = (16 pi G)^{-1} e^{-2d} 2 delta d (4/l): the +4/l shift of B^y",
   Together[-(4/l) (-2 dd0) - 2 dd0 (4/l)]];

gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];

hR = Map[Together, Transpose[JJ . Vinf] . (HR - Hinf) . (JJ . Vbinf), {2}];
hRdisplayed = (1/u) (1 + Lp[xp] Lm[xm]/u^2)/(1 - Lp[xp] Lm[xm]/u^2)^2 *
   {{2 Lp[xp], 2 Lp[xp] Lm[xm], 0}, {2, 2 Lm[xm], 0}, {0, 0, 0}};
NRH`CheckZero["exact h^R_{p qbar} matches the displayed closed form",
   Map[Together, hR - hRdisplayed, {2}]];
NRH`CheckZero["falloffs h^{(2)} = {2L+, 2L+L-; 2, 2L-} + O(e^{-6y/l})",
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

qu = 1/(Sqrt[2] psip[xp] psim[xm] u);
chU = 2 Sqrt[2] ArcTanh[qu];
esig = psim[xm]/psip[xp];

Wu = W0[xp, xm] + W1[xp, xm]/u;
HNR = NonRiemannianH[chU, esig, Wu];
dNR = -1/2 Log[u] + Log[Cosh[chU/(2 Sqrt[2])]];

hNR = Transpose[JJ . Vinf] . (HNR - Hinf) . (JJ . Vbinf);
NRH`Check["exact h^NR = {{e^s sinh chi, (W/2) cosh chi, 0},{0, e^{-s} sinh chi, 0},{0,...}}",
   And[Simplify[hNR[[1, 1]] - esig Sinh[chU]] === 0,
       Simplify[hNR[[1, 2]] - Wu Cosh[chU]/2] === 0,
       Simplify[hNR[[2, 1]]] === 0,
       Simplify[hNR[[2, 2]] - Sinh[chU]/esig] === 0,
       Simplify[hNR[[3, 3]]] === 0,
       Simplify[hNR[[1, 3]]] === 0 && Simplify[hNR[[3, 1]]] === 0]];
NRH`CheckZero["falloffs {2L+ /u, (W0 + W1/u)/2; 0, 2L- /u} with O(u^-3) diagonals",
   Together[{SeriesCoefficient[hNR[[1, 1]], {u, Infinity, 1}] - 2 (1/psip[xp]^2),
    SeriesCoefficient[hNR[[1, 1]], {u, Infinity, 2}],
    SeriesCoefficient[hNR[[1, 2]], {u, Infinity, 0}] - W0[xp, xm]/2,
    SeriesCoefficient[hNR[[1, 2]], {u, Infinity, 1}] - W1[xp, xm]/2,
    SeriesCoefficient[hNR[[2, 2]], {u, Infinity, 1}] - 2 (1/psim[xm]^2),
    SeriesCoefficient[hNR[[2, 2]], {u, Infinity, 2}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 0}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 1}],
    SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 2}]}]];
NRH`CheckZero["delta d_NR = (1/4) L+ L- e^{-4y/l} + O(e^{-8y/l})",
   Together[{SeriesCoefficient[dNR - dinf, {u, Infinity, 2}] - (1/psip[xp]^2) (1/psim[xm]^2)/4,
    SeriesCoefficient[dNR - dinf, {u, Infinity, 1}],
    SeriesCoefficient[dNR - dinf, {u, Infinity, 3}]}]];

hfields = {hpp[xp, xm, u], hpm[xp, xm, u], hmp[xp, xm, u], hmm[xp, xm, u]};
hmat = {{hpp[xp, xm, u], hpm[xp, xm, u], 0}, {hmp[xp, xm, u], hmm[xp, xm, u], 0}, {0, 0, 0}};
VinfU = Vinf . eta3;  VbinfU = Vbinf . etab3;

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

NRH`CheckZero["solves the displayed system",
   Together[eqs /. solNR]];
NRH`CheckZero["reproduces the linearized EDFE: G^{(1)}_MN = 0 on the general solution",
   Map[Together, GLin /. solNR, {2}]];
NRH`Check["the linearized tensor is not empty (it involves the radial derivatives of h)",
   ! FreeQ[GLin, hmp] && ! FreeQ[GLin, dd]];

NRH`CheckZero["D_y{1, e^{-2y/l}} = 0, D_y y = 2/l, D_y y^2 = 2 + 4y/l",
   {Dy[Dy[1]] + 2/l Dy[1],
    Together[Dy[Dy[1/u]] + 2/l Dy[1/u]],
    Together[Dy[Dy[yv]] + 2/l Dy[yv] - 2/l],
    Together[Dy[Dy[yv^2]] + 2/l Dy[yv^2] - 2 - 4 yv/l]}];

NRH`CheckZero["the log coefficients are exact multiples of the three r^(0) conditions",
   {Coefficient[dd[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/16 D[r0[xp, xm], xp, xm],
    Coefficient[hpp[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/4 D[r0[xp, xm], {xp, 2}],
    Coefficient[hmm[xp, xm, u] /. solNR /. Log[u] -> LG, LG] + l^2/4 D[r0[xp, xm], {xm, 2}]}];

NRH`CheckZero["SMcontinuitycheck: d_- h^(2)_{op bop} = 0, d_+ h^(2)_{om bom} = 0, d_pm h^(2)_{om bop} = 0",
   {D[Coefficient[(hpp[xp, xm, u] /. solNR) /. Log[u] -> LG, u, -1], xm],
    D[Coefficient[(hmm[xp, xm, u] /. solNR) /. Log[u] -> LG, u, -1], xp],
    D[Coefficient[(hmp[xp, xm, u] /. solNR) /. Log[u] -> LG, u, -1], xp],
    D[Coefficient[(hmp[xp, xm, u] /. solNR) /. Log[u] -> LG, u, -1], xm]}];

NRH`CheckZero["SMscalarcutoffresponse: 2(16 pi G)^{-1} e^{-2d}(B^y + 4/l) with B^y = 4 d_y d = -4/l + 4 d_y delta d equals (1/(2 pi G)) e^{-2d} d_y delta d",
   Together[2/(16 Pi G) (4 (-1/l + dyd) + 4/l) - 1/(2 Pi G) dyd]];
NRH`CheckZero["d_y delta d = -(l/8) d_+ d_- h^(0)_{om bop} in the general solution, so the coefficient is -(l/(16 pi G)) e^{2y/l} d_+ d_- h^(0)_{om bop}",
   {Together[Dy[dd[xp, xm, u] /. solNR] + l/8 D[r0[xp, xm], xp, xm]],
    Together[1/(2 Pi G) (-(l/8) r0pm) + l/(16 Pi G) r0pm]}];

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

Lq = Together[1/2 D[Lgamma2, {t, 2}] /. t -> 0];
Ltarget = -1/2 u ((2 u/l) D[r1[xp, xm, u], u]) ((2 u/l) D[w1[xp, xm, u], u]);
diffL = Together[Lq - Ltarget];

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
NRH`CheckZero["Gamma^2 quadratic density = -(1/2) e^{2y/l} dy r dy w (mod total derivatives)",
   {EL[r1, diffL], EL[w1, diffL]}];

ronsh = rw0 + rw2/u; wonsh = ws0 + ws2/u;
momR = Limit[-1/2 u (2 u/l) D[wonsh, u], u -> Infinity];
momW = Limit[-1/2 u (2 u/l) D[ronsh, u], u -> Infinity];
NRH`CheckZero["UV momenta give delta S|_UV = (1/(16 pi G l)) Int (w2 dr0 + r2 dw0)",
   {Together[momR - ws2/l], Together[momW - rw2/l]}];

Needs["VariationalMethods`"];
dHs4 = t (chanS[1, 1] hpp[xp, xm, u] + chanS[1, 2] hpm[xp, xm, u] + chanS[2, 1] hmp[xp, xm, u] + chanS[2, 2] hmm[xp, xm, u]);
H4 = Hinf + dHs4 - 1/2 Hinf . JJ . dHs4 . JJ . dHs4;
NRH`CheckZero["the four-channel coset family obeys H J H = J through O(t^2)",
   Map[Together, Normal[Series[H4 . JJ . H4 - JJ, {t, 0, 2}]], {2}]];
L4 = Module[{gamma = GammaDFT[H4, dinf, xs]},
   Gamma2Density[H4, dinf, gamma, xs] - 2 (-2/l^2) Exp[-2 dinf]];
L1of4 = Together[D[L4, t] /. t -> 0];
Lq4 = Together[1/2 D[L4, {t, 2}] /. t -> 0];
fields4 = {hpp, hpm, hmp, hmm};
mom[Lden_, f_] := (l/(2 u)) D[Lden, Derivative[0, 0, 1][f][xp, xm, u]];
NRH`CheckZero["the vacuum carries no first-order radial momentum (no vacuum one-point function)",
   Table[Together[mom[L1of4, f]], {f, fields4}]];

onsh4 = {hmp -> Function[{a, b, c}, r0c + r2c/c],
   hpp -> Function[{a, b, c}, s0p[a, b] + s2p[a]/c],
   hmm -> Function[{a, b, c}, s0m[a, b] + s2m[b]/c],
   hpm -> Function[{a, b, c}, w0[a, b] - l/2 (l/2 Log[c]) (D[s0m[a, b], {a, 2}] + D[s0p[a, b], {b, 2}]) + w2[a, b]/c]};
NRH`CheckZero["the restricted on-shell sector solves the displayed system (delta d = 0)",
   Together[(eqs /. dd -> Function[{a, b, c}, 0]) /. onsh4]];

srcVar = {r0c -> r0c + ee dr0c, s0p -> Function[{a, b}, s0p[a, b] + ee ds0p[a, b]],
   s0m -> Function[{a, b}, s0m[a, b] + ee ds0m[a, b]], w0 -> Function[{a, b}, w0[a, b] + ee dw0[a, b]]};
hOn[f_] := f[xp, xm, u] /. onsh4;
dhOn[f_] := D[hOn[f] /. srcVar, ee] /. ee -> 0;
ThetaY = Together[Sum[(mom[Lq4, f] /. onsh4) dhOn[f], {f, fields4}]] /. Log[u] -> LG;
ThetaSer = Normal[Series[ThetaY, {u, Infinity, 0}]];
NRH`Check["Theta^y at the cutoff is a Laurent polynomial in u with at most linear-in-y (log) pieces",
   PolynomialQ[Together[ThetaSer u^2], {u, LG}] && Exponent[ThetaSer, LG] <= 1 && Exponent[ThetaSer, u] <= 1];
varFns = {ds0p[xp, xm], ds0m[xp, xm], dw0[xp, xm]};
totalDerivQ[dens_] := NRH`ZeroQ[Together[VariationalD[dens, varFns, {xp, xm}]]] &&
   NRH`ZeroQ[Together[VariationalD[Coefficient[dens, dr0c, 1], {s0p[xp, xm], s0m[xp, xm], w0[xp, xm], w2[xp, xm]}, {xp, xm}]]];
divPieces = {Coefficient[Coefficient[ThetaSer, u, 1], LG, 0], Coefficient[Coefficient[ThetaSer, u, 1], LG, 1],
   Coefficient[Coefficient[ThetaSer, u, 0], LG, 1]};
NRH`Check["the divergent and linear-in-y pieces of Theta^y are total x-derivatives in the restricted sector",
   AllTrue[divPieces, totalDerivQ]];
finitePiece = Coefficient[Coefficient[ThetaSer, u, 0], LG, 0];
target26 = (1/l) (s2m[xm] ds0p[xp, xm] + s2p[xp] ds0m[xp, xm] + w2[xp, xm] dr0c + r2c dw0[xp, xm]);
NRH`Check["the finite piece equals (1/l)[h^(2)_{om bom} dh^(0)_{op bop} + h^(2)_{op bop} dh^(0)_{om bom} + h^(2)_{op bom} dh^(0)_{om bop} + h^(2)_{om bop} dh^(0)_{op bom}] (times (16 pi G)^{-1}), mod total x-derivatives",
   totalDerivQ[Together[finitePiece - target26]]];
NRH`Check["the finite piece is NOT a total derivative by itself (the pairing is genuine)",
   ! totalDerivQ[finitePiece]];

NRH`CheckZero["<K_{a bbar}> = h^{(2)}_{a bbar}/(32 pi G l) => stress one-points",
   {Together[2 Lp[xp]/(32 Pi G l) - Lp[xp]/(16 Pi G l)],
    Together[2 Lm[xm]/(32 Pi G l) - Lm[xm]/(16 Pi G l)],
    Together[2/(32 Pi G l) - 1/(16 Pi G l)],
    Together[2 Lp[xp] Lm[xm]/(32 Pi G l) - Lp[xp] Lm[xm]/(16 Pi G l)]}];
NRH`CheckZero["NR hair channel: h^{(2)}_{op bom} = W_1/2 => <K_{op bom}> = W_1/(64 pi G l)",
   Together[(W1[xp, xm]/2)/(32 Pi G l) - W1[xp, xm]/(64 Pi G l)]];
NRH`CheckZero["h^(2)_{om bop} = 0 on the NR family gives <K_{om bop}> = 0; (1/(64 pi G l)) delta(W_1/2) = (1/(128 pi G l)) delta W_1",
   {Together[SeriesCoefficient[hNR[[2, 1]], {u, Infinity, 1}]/(32 Pi G l)], Together[1/(64 Pi G l)/2 - 1/(128 Pi G l)]}];

xsZ = {xp, xm, Function[e, -(2 z/l) D[e, z]]};
HNRz = NRBoundaryH[Lp[xp], Lm[xm], W1[xp, xm], z];
dNRz = NRBoundaryD[Lp[xp], Lm[xm], z];
NRH`Check["NR setup: the z^2 truncation obeys H J H = J through O(z^2)",
   Module[{c = Expand[HNRz . JJ . HNRz - JJ]},
      AllTrue[Flatten[c], PossibleZeroQ[Coefficient[#, z, 0]] && PossibleZeroQ[Coefficient[#, z, 1]] && PossibleZeroQ[Coefficient[#, z, 2]] &]]];
AyNR = MomentumAK[HNRz, GammaDFT[HNRz, dNRz, xsZ], xsZ][[6]];
AfixNR = Map[Together, Transpose[JJ . Vinf] . AyNR . (JJ . Vbinf), {2}];
KNR = Map[Limit[-(1/(32 Pi G)) #/z, z -> 0] &, AfixNR[[1 ;; 2, 1 ;; 2]], {2}];
NRH`CheckZero["on the exact NR family: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, W1/4},{0, L-}}/(16 pi G l) [ without linearization]",
   Map[Together, KNR - {{Lp[xp], W1[xp, xm]/4}, {0, Lm[xm]}}/(16 Pi G l), {2}]];
NRH`CheckZero["on the exact NR family: e^{2Y/l}(B^y + 4/l) -> 0, hence <T_(0)> = 0",
   Limit[(GammaBVector[HNRz, dNRz, xsZ][[6]] + 4/l)/z, z -> 0]];

eta2 = eta3[[1 ;; 2, 1 ;; 2]]; etab2 = etab3[[1 ;; 2, 1 ;; 2]];
hlow2 = {{h0pp, h0pm}, {h0mp, h0mm}};
hup2 = eta2 . hlow2 . etab2;
NRH`CheckZero["J^{op bop} = -2 h^(0)op bop = 2 h^(0)_{om bom} and J^{om bom} = 2 h^(0)_{op bop}; hence the Hessian factor (32 pi G l)^{-1}/2 = (64 pi G l)^{-1}",
   {Together[-2 hup2[[1, 1]] - 2 h0mm], Together[-2 hup2[[2, 2]] - 2 h0pp], Together[1/(32 Pi G l)/2 - 1/(64 Pi G l)]}];

xiPBH = {
   0,
   +l^2/(2 u) Lm[xm] D[al[xp, xm], {xp, 2}],
   -l/2 D[al[xp, xm], xp],
   al[xp, xm] ,
   l^2/(4 u) D[al[xp, xm], {xp, 2}],
   -l/2 D[al[xp, xm], xp]} /. Lm[xm] -> 0;

DeltaH = Map[Together, GenLieH[xiPBH, Hinf, xs], {2}];
hPBH = Transpose[JJ . Vinf] . DeltaH . (JJ . Vbinf);
NRH`CheckZero["s_+ = [Delta H]^{(0)}_{om bom} = -2 d_- alpha^+",
   Together[SeriesCoefficient[hPBH[[2, 2]], {u, Infinity, 0}] + 2 D[al[xp, xm], xm]]];
NRH`CheckZero["r_+ = [Delta H]^{(2)}_{op bop} = -(l^2/2) d_+^3 alpha^+",
   Together[SeriesCoefficient[hPBH[[1, 1]], {u, Infinity, 1}] + l^2/2 D[al[xp, xm], {xp, 3}]]];
NRH`CheckZero["the PBH variation carries no (om bop) source deformation",
   Together[SeriesCoefficient[hPBH[[2, 1]], {u, Infinity, 0}]]];

NRH`CheckZero["(l^2/4) d_+^3 [-1/(2 pi x)] = 3 l^2/(4 pi x^4)",
   Together[l^2/4 D[-1/(2 Pi x), {x, 3}] - 3 l^2/(4 Pi x^4)]];
NRH`CheckZero["away from coincidence: d_+ d_- ln(-x^+ x^-) = 0",
   D[Log[-xps xms], xps, xms]];
centralCharge = 3 l/(2 G);
NRH`CheckZero["(64 pi G l)^{-1} 3 l^2/(4 pi) = (8 pi)^{-2} (c/2), c = 3l/2G",
   Together[1/(64 Pi G l) 3 l^2/(4 Pi) - 1/(8 Pi)^2 centralCharge/2]];

(* ::Section:: *)
(*Current source variations and particular two-point kernels*)

(* SMmomentumvariation: differentiate the unprojected connection tensor and both frames. *)
Module[{core, gamma, v, vb, dv, dvb, a0, da, moving, frozen, frames, b},
   core[g_] := MomentumCore[g, xs];
   gamma = GammaDFT[Hlin, dlin, xs];
   v = JJ . Vinf; vb = JJ . Vbinf;
   dv = vb . etab3 . Transpose[hmat]/2; dvb = -v . eta3 . hmat/2;
   a0 = core[gamma /. t -> 0]; da = D[core[gamma], t] /. t -> 0;
   frozen = Transpose[v] . da[[6]] . vb;
   frames = Transpose[dv] . a0[[6]] . vb + Transpose[v] . a0[[6]] . dvb;
   moving = D[Transpose[v + t dv] . core[gamma][[6]] . (vb + t dvb), t] /. t -> 0;
   NRH`CheckZero["SMmomentumvariation: moving-frame product rule", moving - frozen - frames];
   NRH`CheckZero["SMlinearizedmomenta: vacuum frame term = -h/l", frames + hmat/l];
   NRH`CheckZero["SMlinearizedmomenta: vacuum tangential connection term = dy h/2 + h/l", (frozen - Dy[hmat]/2 - hmat/l)[[1 ;; 2, 1 ;; 2]]];
   NRH`CheckZero["SMlinearizedmomenta: full vacuum tangential delta A = dy h/2", (moving - Dy[hmat]/2)[[1 ;; 2, 1 ;; 2]]];
   NRH`Check["negative control: freezing the bulk frames changes the response", ! NRH`ZeroQ[frozen - moving]];
   b = D[GammaBVector[Hlin, dlin, xs], t] /. t -> 0;
   NRH`CheckZero["SMmomentumvariation: delta B product rule",
      b - 4 (JJ . deltaH . JJ) . DblGrad[dinf, xs] - 4 (JJ . Hinf . JJ) . DblGrad[dd[xp, xm, u], xs]
        + Table[Sum[DblD[(JJ . deltaH . JJ)[[k, j]], j, xs], {j, 6}], {k, 6}]];
   NRH`CheckZero["SMlinearizedmomenta: vacuum delta B^y = 4 dy delta d", b[[6]] - 4 Dy[dd[xp, xm, u]]];
];

(* SMsecondvariation: independent flat h arguments in local coset coordinates. *)
Module[{j, s, h1, h2, generator, t1, t2, h12, dh1, dh2, cross, v, vb},
   j = ArrayFlatten[{{eta3, 0}, {0, etab3}}];
   s = ArrayFlatten[{{eta3, 0}, {0, -etab3}}];
   h1 = Array[a, {3, 3}]; h2 = Array[b, {3, 3}];
   generator[h_] := ArrayFlatten[{{ConstantArray[0, {3, 3}], -eta3 . h/2},
      {etab3 . Transpose[h]/2, ConstantArray[0, {3, 3}]}}];
   t1 = generator[h1]; t2 = generator[h2]; cross = (t1 . t2 + t2 . t1)/2;
   dh1 = t1 . s + s . Transpose[t1]; dh2 = t2 . s + s . Transpose[t2];
   h12 = cross . s + t1 . s . Transpose[t2] + t2 . s . Transpose[t1] + s . Transpose[cross];
   v = IdentityMatrix[6][[All, 1 ;; 3]]; vb = IdentityMatrix[6][[All, 4 ;; 6]];
   NRH`CheckZero["SMsecondvariation: mixed coset constraint at second order",
      h12 . j . s + s . j . h12 + dh1 . j . dh2 + dh2 . j . dh1];
   NRH`CheckZero["SMsecondvariation: h1 projection and delta2 h1 = 0",
      {Transpose[v] . j . dh1 . j . vb - h1,
       Transpose[t2 . v] . j . dh1 . j . vb + Transpose[v] . j . h12 . j . vb
        + Transpose[v] . j . dh1 . j . (t2 . vb)}];
   NRH`Check["negative control: independent h does not mean delta2 delta1 H = 0", ! NRH`ZeroQ[h12]];
];

Module[{h, d, metricFinite, scalarFinite},
   h = Exp[-2 y/l] (h2 + y h2L + y^2 h2LL);
   d = Exp[-2 y/l] (d2 + y d2L);
   metricFinite = Coefficient[Expand[-Exp[2 y/l] D[h, y]/(64 Pi G)], y, 0];
   scalarFinite = Coefficient[Expand[Exp[2 y/l] D[d, y]/(4 Pi G)], y, 0];
   NRH`CheckZero["SMgeneralresponsevariation: finite metric and scalar logarithmic coefficients",
      {metricFinite - h2/(32 Pi G l) + h2L/(64 Pi G),
       scalarFinite + d2/(2 Pi G l) - d2L/(4 Pi G)}];
];

(* Rcorrelators, SMgeneralWardoperators: apply verified Ward laws to the chosen inverse derivative.
   Delta is the separation; Lp[xp], W1[xp,xm] remain evaluated at the first point.
   The result checks particular kernels, not interior/zero-mode or local-contact completions. *)
Module[{inverse, wardNR, wardR, wardHair, aNR, aR, mNR, mR, targetA, targetM},
   inverse = -1/(2 Pi delta);
   wardNR[f_] := D[Lp[xp], xp] f + 2 Lp[xp] D[f, delta];
   wardR[f_] := wardNR[f] - l^2 D[f, {delta, 3}]/4;
   wardHair[f_] := D[W1[xp, xm], xp] f + 2 W1[xp, xm] D[f, delta] - l^2 Lm[xm] D[f, {delta, 3}];
   aNR = -wardNR[inverse]/(64 Pi G l); aR = -wardR[inverse]/(64 Pi G l);
   mNR = -wardHair[inverse]/(256 Pi G l);
   mR = -Lm[xm] wardR[inverse]/(64 Pi G l);
   targetA = (D[Lp[xp], xp]/delta - 2 Lp[xp]/delta^2)/(128 Pi^2 G l);
   targetM = (D[W1[xp, xm], xp]/delta - 2 W1[xp, xm]/delta^2 + 6 l^2 Lm[xm]/delta^4)/(512 Pi^2 G l);
   NRH`CheckZero["Rcorrelators: NR particular diagonal kernel", aNR - targetA];
   NRH`CheckZero["Rcorrelators: R particular diagonal kernel", aR - targetA - 3 l/(256 Pi^2 G delta^4)];
   NRH`CheckZero["Rcorrelators: NR particular non-diagonal kernel M", mNR - targetM];
   NRH`CheckZero["Rcorrelators: R particular non-diagonal kernel M", mR - Lm[xm] aR];
   NRH`Check["negative control: reversing the source sign fails the NR kernel", ! NRH`ZeroQ[aNR + targetA]];
   NRH`CheckZero["Rcorrelators: NR fourth-order pole is state dependent", D[Coefficient[Expand[mNR], delta, -4], Lm[xm]] - 3 l/(256 Pi^2 G)];
];

(* SMlinearizedmomenta: keep the finite-state backgrounds through u^-1. *)
Module[{firstOrder, backgroundH, v, vb, dv, dvb, dh, gamma, a, moving, frame, scalar, hBackground},
   firstOrder[e_] := Map[Together, Normal[Series[e, {s, 0, 1}]], {ArrayDepth[e]}];
   Do[
      hBackground = data[[2]]/u;
      v = JJ . (Vinf + s Vbinf . etab3 . Transpose[hBackground]/2);
      vb = JJ . (Vbinf - s Vinf . eta3 . hBackground/2);
      backgroundH = Hinf + s ((Vinf . eta3) . hBackground . Transpose[Vbinf . etab3]
         + (Vbinf . etab3) . Transpose[hBackground] . Transpose[Vinf . eta3]);
      dh = firstOrder[(JJ . v . eta3) . hmat . Transpose[JJ . vb . etab3]];
      dh = dh + Transpose[dh];
      gamma = GammaDFT[backgroundH + t dh, dinf + t dd[xp, xm, u], xs];
      a = MomentumCore[gamma, xs][[6]];
      dv = vb . etab3 . Transpose[hmat]/2; dvb = -v . eta3 . hmat/2;
      moving = firstOrder[D[Transpose[v + t dv] . a . (vb + t dvb), t] /. t -> 0];
      frame = firstOrder[Transpose[dv] . (a /. t -> 0) . vb + Transpose[v] . (a /. t -> 0) . dvb];
      scalar = firstOrder[D[GammaBVector[backgroundH + t dh, dinf + t dd[xp, xm, u], xs][[6]], t] /. t -> 0];
      NRH`CheckZero["SMlinearizedmomenta: " <> data[[1]] <> " tangential delta A through u^-1",
         (moving - Dy[hmat]/2)[[1 ;; 2, 1 ;; 2]]];
      NRH`CheckZero["SMlinearizedmomenta: " <> data[[1]] <> " frame term through u^-1", (frame + hmat/l)[[1 ;; 2, 1 ;; 2]]];
      NRH`CheckZero["SMlinearizedmomenta: " <> data[[1]] <> " scalar delta B through u^-1", scalar - 4 Dy[dd[xp, xm, u]]],
      {data, {{"R", {{2 Lp[xp], 2 Lp[xp] Lm[xm], 0}, {2, 2 Lm[xm], 0}, {0, 0, 0}}},
              {"NR", {{2 Lp[xp], W1[xp, xm]/2, 0}, {0, 2 Lm[xm], 0}, {0, 0, 0}}}}}];
];

NRH`FileSummary[];
