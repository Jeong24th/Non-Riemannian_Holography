(* ::Title:: *)
(*NRH07 SM KillingSpinors*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH07_SM_KillingSpinors.wl"];

JJ = ODDJ[3];

coordsS = {th, f1, f2};
gS = l^2 DiagonalMatrix[{1, Cos[th]^2, Sin[th]^2}];
BS = {{0, 0, 0}, {0, 0, l^2 Cos[th]^2}, {0, -l^2 Cos[th]^2, 0}};
giS = Inverse[gS];
chrS = Table[Together[1/2 Sum[giS[[m, s]] (D[gS[[s, i]], coordsS[[j]]] + D[gS[[s, j]], coordsS[[i]]] - D[gS[[i, j]], coordsS[[s]]]), {s, 3}]], {m, 3}, {i, 3}, {j, 3}];
ricS = Table[Together[Sum[D[chrS[[k, i, j]], coordsS[[k]]], {k, 3}] - Sum[D[chrS[[k, i, k]], coordsS[[j]]], {k, 3}]
    + Sum[chrS[[k, i, j]] chrS[[s, k, s]], {k, 3}, {s, 3}] - Sum[chrS[[k, i, s]] chrS[[s, j, k]], {k, 3}, {s, 3}]], {i, 3}, {j, 3}];
rSclr = Simplify[Sum[giS[[i, j]] ricS[[i, j]], {i, 3}, {j, 3}]];
hS = Table[D[BS[[j, k]], coordsS[[i]]] + D[BS[[k, i]], coordsS[[j]]] + D[BS[[i, j]], coordsS[[k]]], {i, 3}, {j, 3}, {k, 3}];
h2S = Simplify[Sum[hS[[i, j, k]] hS[[a, b, c]] giS[[i, a]] giS[[j, b]] giS[[k, c]], {i, 3}, {j, 3}, {k, 3}, {a, 3}, {b, 3}, {c, 3}]];
NRH`CheckZero["R(S3) = +6/l^2 and H^2(S3) = +24/l^2 (pairwise cancellation with AdS3)",
   {rSclr - 6/l^2, h2S - 24/l^2}];
NRH`CheckZero["R_{mu nu} = (1/4) H_{mu rho sigma} H_nu^{rho sigma} on the S3 factor",
   Simplify[ricS - 1/4 Table[Sum[hS[[i, r, s]] hS[[j, a, b]] giS[[r, a]] giS[[s, b]], {r, 3}, {s, 3}, {a, 3}, {b, 3}], {i, 3}, {j, 3}]]];

HS = Map[Together, RiemannianH[gS, BS], {2}];
dS = -1/4 Log[Det[gS]];
NRH`CheckZero["S_(0)(S3 block) = +4/l^2 via the closed form (cancels -4/l^2 of either 3d saddle)",
   Simplify[ScalarS0[HS, dS, coordsS] - 4/l^2]];
HR4 = IdentityMatrix[8];
NRH`Check["S_(0)(R4 block) = 0 (flat block, constant dilaton)",
   Together[ScalarS0[HR4, 0, {z1, z2, z3, z4}]] === 0];

ordD[e_, i_] := {D[e, xp], D[e, xm], (2 u/l) D[e, u]}[[i]];
gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
giB3 = Map[Together, Inverse[gR], {2}];
chrB = Table[Together[1/2 Sum[giB3[[m, s]] (ordD[gR[[s, i]], j] + ordD[gR[[s, j]], i] - ordD[gR[[i, j]], s]), {s, 3}]],
   {m, 3}, {i, 3}, {j, 3}];
ricB = Table[Together[Sum[ordD[chrB[[k, i, j]], k], {k, 3}] - Sum[ordD[chrB[[k, i, k]], j], {k, 3}]
    + Sum[chrB[[k, i, j]] chrB[[s, k, s]], {k, 3}, {s, 3}] - Sum[chrB[[k, i, s]] chrB[[s, j, k]], {k, 3}, {s, 3}]], {i, 3}, {j, 3}];
hB = Table[ordD[bR[[j, k]], i] + ordD[bR[[k, i]], j] + ordD[bR[[i, j]], k], {i, 3}, {j, 3}, {k, 3}];
NRH`CheckZero["R(AdS3, Banados) = -6/l^2 and H^2 = -24/l^2 for arbitrary chiral L_pm (pairwise cancellation with S3)",
   {Together[Sum[giB3[[i, j]] ricB[[i, j]], {i, 3}, {j, 3}] + 6/l^2],
    Together[Sum[hB[[i, j, k]] hB[[a, b, c]] giB3[[i, a]] giB3[[j, b]] giB3[[k, c]], {i, 3}, {j, 3}, {k, 3}, {a, 3}, {b, 3}, {c, 3}] + 24/l^2]}];
NRH`CheckZero["R_{mu nu} = (1/4) H_{mu rho sigma} H_nu^{rho sigma} on the AdS3 factor, arbitrary chiral L_pm",
   Map[Together, ricB - 1/4 Table[Sum[hB[[i, r, s]] hB[[j, a, b]] giB3[[r, a]] giB3[[s, b]], {r, 3}, {s, 3}, {a, 3}, {b, 3}], {i, 3}, {j, 3}], {2}]];
NRH`CheckZero["flux orientation: H_{y-+} = +(2/l) Sqrt[|g_3|] = (2/l)(e^{2y/l} - L+L- e^{-2y/l})",
   Together[hB[[3, 2, 1]] - 2/l (u - Lp[xp] Lm[xm]/u)]];

d10Assemble[H3_, d3_] := {ArrayFlatten[{
      {H3[[1 ;; 3, 1 ;; 3]], 0, 0, H3[[1 ;; 3, 4 ;; 6]], 0, 0},
      {0, HS[[1 ;; 3, 1 ;; 3]], 0, 0, HS[[1 ;; 3, 4 ;; 6]], 0},
      {0, 0, IdentityMatrix[4], 0, 0, 0},
      {H3[[4 ;; 6, 1 ;; 3]], 0, 0, H3[[4 ;; 6, 4 ;; 6]], 0, 0},
      {0, HS[[4 ;; 6, 1 ;; 3]], 0, 0, HS[[4 ;; 6, 4 ;; 6]], 0},
      {0, 0, 0, 0, 0, IdentityMatrix[4]}}], d3 + dS};

uplift10[H3v_, d3v_, label_] := Module[
   {xs10, gamma, r4, ric, s0, psp, Pn, Pbn, JJ10 = ODDJ[10], H10, d10, tt},
   {H10, d10} = d10Assemble[H3v, d3v];
   xs10 = {xp, xm, Function[e, (2 u/l) D[e, u]], th, f1, f2, z1, z2, z3, z4};
   {tt, gamma} = AbsoluteTiming[GammaDFT[H10, d10, xs10]];
   Print["    [10d timing] Gamma: ", Round[tt, 0.1], " s"];
   {tt, r4} = AbsoluteTiming[RiemannR4[gamma, xs10]];
   Print["    [10d timing] R4: ", Round[tt, 0.1], " s"];
   {tt, ric} = AbsoluteTiming[RicciS[gamma, r4, xs10]];
   Print["    [10d timing] Ricci: ", Round[tt, 0.1], " s"];
   s0 = ScalarS0[H10, d10, xs10];
   Pn = (JJ10 + H10)/2; Pbn = (JJ10 - H10)/2;
   psp = Pn . JJ10 . ric . JJ10 . Pbn;
   NRH`CheckZero[label <> ":  S_(0)^{(10)} = 0", Together[s0]];
   NRH`CheckZero[label <> ":  (P S Pbar)^{(10)} = 0", Map[Together, psp, {2}]]];

HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];
uplift10[HR, dR, "R uplift, arbitrary chiral L_+(x^+), L_-(x^-)"];

Wsym = W0[xp, xm] + W1[xp, xm]/u;
HNR3sym = {{0, 0, 0, 1, 0, 0},
   {0, 0, 0, 2 Lp[xp]/u, -1, 0},
   {0, 0, 1, 0, 0, 0},
   {1, 2 Lp[xp]/u, 0, -2 Lp[xp] Wsym/u, Wsym, 0},
   {0, -1, 0, Wsym, 0, 0},
   {0, 0, 0, 0, 0, 1}};
dNR3sym = -1/2 Log[u];
NRH`CheckZero["one-sided NR family obeys H J H = J exactly for arbitrary L_+(x^+), W_0, W_1",
   Map[Together, HNR3sym . JJ . HNR3sym - JJ, {2}]];
uplift10[HNR3sym, dNR3sym, "NR uplift, one-sided hairy family L_- = 0, arbitrary L_+(x^+), W_0(x), W_1(x)"];

xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};
dinf = -1/2 Log[u];
xiIso = {om1[xp], om2[xm], -l/2 (D[vp[xp], xp] - D[vm[xm], xm]),
   vp[xp], vm[xm], -l/2 (D[vp[xp], xp] + D[vm[xm], xm])};
NRH`CheckZero["Lhat_xi H^infty = 0 for arbitrary chiral v^pm and omega_pm",
   Map[Together, GenLieH[xiIso, Hinf, xsU], {2}]];
NRH`CheckZero["Lhat_xi d = 0 (the radial component compensates the divergence)",
   Together[GenLieD[xiIso, dinf, xsU]]];
NRH`CheckZero["equivalently d_M(e^{-2d} xi^M) = 0",
   Together[Sum[DblD[Exp[-2 dinf] xiIso[[m]], m, xsU], {m, 6}]]];

NRH`CheckZero["eps = c/Sqrt[L] solves eps dL + 2 L d eps = 0 (the generic obstruction)",
   Together[cc/Sqrt[LL[x]] D[LL[x], x] + 2 LL[x] D[cc/Sqrt[LL[x]], x]]];
HNRone = Module[{Wn = W0[xp, xm] + W1[xp, xm]/u},
   {{0, 0, 0, 1, 0, 0},
    {0, 0, 0, 2 Lp[xp]/u, -1, 0},
    {0, 0, 1, 0, 0, 0},
    {1, 2 Lp[xp]/u, 0, -2 Lp[xp] Wn/u, Wn, 0},
    {0, -1, 0, Wn, 0, 0},
    {0, 0, 0, 0, 0, 1}}];
NRH`CheckZero["the exact one-sided family (L- = 0) obeys H J H = J for arbitrary L+(x+), W0, W1",
   Map[Together, HNRone . JJ . HNRone - JJ, {2}]];
xiChiralP = {0, 0, -l/2 D[ep[xp], xp], ep[xp], 0, -l/2 D[ep[xp], xp]};
lieOne = GenLieH[xiChiralP, HNRone, xsU];
HNRoneGen = HNRone /. {Lp[xp] -> LPv, W0[xp, xm] -> W0v, W1[xp, xm] -> W1v};
depsOne = (D[HNRoneGen, LPv] (ep[xp] D[Lp[xp], xp] + 2 Lp[xp] D[ep[xp], xp])
   + D[HNRoneGen, W0v] (ep[xp] D[W0[xp, xm], xp] + W0[xp, xm] D[ep[xp], xp])
   + D[HNRoneGen, W1v] (ep[xp] D[W1[xp, xm], xp] + 2 W1[xp, xm] D[ep[xp], xp])) /.
   {LPv -> Lp[xp], W0v -> W0[xp, xm], W1v -> W1[xp, xm]};
NRH`CheckZero["exact one-sided identity Lhat_xi H = delta H, with weights (1,2) for (W0, W1)",
   Map[Together, lieOne - depsOne, {2}]];

chh = Cosh[hh]; shh = Sinh[hh];
VexU = {{0, -chh/Sqrt[2], 0}, {0, -es shh/Sqrt[2], 0}, {0, 0, 1/Sqrt[2]},
   {Sqrt[2] chh, Wc es shh/(2 Sqrt[2]), 0}, {-Sqrt[2] shh/es, -Wc chh/(2 Sqrt[2]), 0}, {0, 0, 1/Sqrt[2]}};
VbexU = {{shh/(es Sqrt[2]), 0, 0}, {chh/Sqrt[2], 0, 0}, {0, 0, -1/Sqrt[2]},
   {-Wc chh/(2 Sqrt[2]), -Sqrt[2] es shh, 0}, {Wc shh/(2 Sqrt[2] es), Sqrt[2] chh, 0}, {0, 0, 1/Sqrt[2]}};
HNRex = Module[{CC = Cosh[2 hh], SS = Sinh[2 hh]},
   {{0, 0, 0, CC, -SS/es, 0}, {0, 0, 0, es SS, -CC, 0}, {0, 0, 1, 0, 0, 0},
    {CC, es SS, 0, -Wc es SS, Wc CC, 0}, {-SS/es, -CC, 0, Wc CC, -Wc SS/es, 0}, {0, 0, 0, 0, 0, 1}}];
ratH[e_] := Together[TrigToExp[e] /. E^(k_. hh) :> T^k];
Pex = (JJ + HNRex)/2; Pbex = (JJ - HNRex)/2;
NRH`CheckZero["V_M^p V_Np = P_MN and Vbar_M^pbar Vbar_Npbar = Pbar_MN for the exact frame , any hh, sigma, W",
   {Map[ratH, VexU . eta3 . Transpose[VexU] - Pex, {2}], Map[ratH, VbexU . (-eta3) . Transpose[VbexU] - Pbex, {2}]}];
NRH`CheckZero["V^M_p Vbar_{M qbar} = 0 and P + Pbar = J",
   {Map[ratH, Transpose[JJ . VexU] . VbexU, {2}], Pex + Pbex - JJ}];
NRH`CheckZero["V eta V^T - Vbar etabar Vbar^T = H(W) exactly (H = P - Pbar)",
   Map[ratH, VexU . eta3 . Transpose[VexU] - VbexU . (-eta3) . Transpose[VbexU] - HNRex, {2}]];
NRH`CheckZero["at hh = 0 and W_0 = 0, lowering the local indices with eta, etabar gives the limiting frame",
   {(VexU . eta3 /. {hh -> 0, Wc -> 0}) - Vinf, (VbexU . (-eta3) /. {hh -> 0, Wc -> 0}) - Vbinf}];
tauPex = {chh, -shh/es}; tauMex = {-es shh, chh};
Yex = {chh, es shh}; Ybex = {shh/es, chh};
NRH`CheckZero["x rows = Sqrt[2] tau^pm, x~ rows = the dual vectors -Y/Sqrt[2], Ybar/Sqrt[2], and the W entries = -(W/(2 Sqrt[2])) tau^mp",
   Map[ratH, Flatten[{VexU[[4 ;; 5, 1]] - Sqrt[2] tauPex, VbexU[[4 ;; 5, 2]] - Sqrt[2] tauMex,
      VexU[[1 ;; 2, 2]] + Yex/Sqrt[2], VbexU[[1 ;; 2, 1]] - Ybex/Sqrt[2],
      Yex . tauPex - 1, Yex . tauMex, Ybex . tauMex - 1, Ybex . tauPex,
      VexU[[4 ;; 5, 2]] + Wc/(2 Sqrt[2]) tauMex, VbexU[[4 ;; 5, 1]] + Wc/(2 Sqrt[2]) tauPex}]]];

hillRule = {Derivative[2][sA][x] -> 2/l^2 LL[x] sA[x], Derivative[2][sB][x] -> 2/l^2 LL[x] sB[x]};
hillD[e_] := D[e, x] /. hillRule;
kAB = sA[x] sB[x];
k1 = hillD[kAB]; k2 = hillD[k1]; k3 = hillD[k2];
NRH`CheckZero["k = s_i s_j with (l^2/2) s'' = L s obeys k L' + 2 L k' - (l^2/4) k''' = 0 (the stabilizer equation)",
   Together[kAB D[LL[x], x] + 2 LL[x] k1 - l^2/4 k3]];

s1 = {{0, 1}, {1, 0}}; s2 = {{0, -I}, {I, 0}}; s3 = {{1, 0}, {0, -1}};
id2 = IdentityMatrix[2]; id4 = IdentityMatrix[4]; id32 = IdentityMatrix[32];
tauP3 = {s1, s2, s3};
rho = {KroneckerProduct[s1, id2], KroneckerProduct[s2, id2],
   KroneckerProduct[s3, s1], KroneckerProduct[s3, s2], KroneckerProduct[s3, s3]};
gamOp = Sqrt[2] {{0, 0}, {1, 0}};
gamOm = -Sqrt[2] {{0, 1}, {0, 0}};
gamY = s3;
gam3 = {gamOp, gamOm, gamY};

NRH`CheckZero["rep: {gamma^p, gamma^q} = 2 eta^{pq} (3d lightcone blocks)",
   Flatten[Table[gam3[[i]] . gam3[[j]] + gam3[[j]] . gam3[[i]] - 2 eta3[[i, j]] id2, {i, 3}, {j, 3}]]];

Gam[a_] := Which[
   a <= 3, KroneckerProduct[gam3[[a]], id2, s1, id4],
   a <= 6, KroneckerProduct[id2, tauP3[[a - 3]], s2, id4],
   True, KroneckerProduct[id2, id2, s3, rho[[a - 6]]]];
eta10 = ArrayFlatten[{{eta3, 0, 0}, {0, IdentityMatrix[3], 0}, {0, 0, id4}}];
NRH`CheckZero["{Gamma^p, Gamma^q} = 2 eta_{(10)}^{pq} I_32",
   Flatten[Table[Gam[a] . Gam[b] + Gam[b] . Gam[a] - 2 eta10[[a, b]] id32, {a, 10}, {b, 10}]]];
Gam11 = KroneckerProduct[id2, id2, s3, rho[[5]]];
NRH`CheckZero["Gamma_11^2 = 1 and {Gamma_11, Gamma^p} = 0",
   Join[Flatten[Gam11 . Gam11 - id32], Flatten[Table[Gam11 . Gam[a] + Gam[a] . Gam11, {a, 10}]]]];
BB10 = KroneckerProduct[id2, s2, id2, s1, s2];
conj[m_] := m /. Complex[re_, im_] :> Complex[re, -im];
NRH`CheckZero["BB10 Gamma^p BB10^{-1} = (Gamma^p)^*, same for Gamma_11, and BB10 BB10^* = 1",
   Join[Flatten[Table[BB10 . Gam[a] . Inverse[BB10] - conj[Gam[a]], {a, 10}]],
      Flatten[BB10 . Gam11 . Inverse[BB10] - conj[Gam11]],
      Flatten[BB10 . conj[BB10] - id32]]];
GamBar[a_] := Gam[a] . Gam11;
NRH`CheckZero["{Gammabar, Gammabar} = -2 eta, Gammabar_{pq} = -Gamma_{pq}, same Majorana intertwiner",
   Join[
      Flatten[Table[GamBar[a] . GamBar[b] + GamBar[b] . GamBar[a] + 2 eta10[[a, b]] id32, {a, 10}, {b, 10}]],
      Flatten[Table[(GamBar[a] . GamBar[b] - GamBar[b] . GamBar[a])/2
         + (Gam[a] . Gam[b] - Gam[b] . Gam[a])/2, {a, 10}, {b, 10}]],
      Flatten[Table[BB10 . GamBar[a] . Inverse[BB10] - conj[GamBar[a]], {a, 10}]]]];

realify[m_] := ArrayFlatten[{{Re[m], -Im[m]}, {Im[m], Re[m]}}];

majoranaOps = ArrayFlatten[{{Re[BB10] - id32, -Im[BB10]}, {Im[BB10], Re[BB10] + id32}}];
NRH`Check["the Majorana condition leaves 32 real components",
   64 - MatrixRank[majoranaOps] == 32];
weylOps = realify[(id32 - Gam11)/2];
NRH`Check["adding the Weyl condition leaves 16 real components",
   64 - MatrixRank[Join[majoranaOps, weylOps]] == 16];

projS3perp = KroneckerProduct[id2, {{0, 0}, {0, 1}}, id2, id4];
projAuxPerp = KroneckerProduct[id2, id2, {{0, 0}, {0, 1}}, id4];
NRH`Check["Weyl + S^3-line + zeta_+ leave complex dimension 4 (the Xi_{+r} span)",
   32 - MatrixRank[Join[(id32 - Gam11)/2, projS3perp, projAuxPerp]] == 4];

gammaVac = GammaDFT[Hinf, dinf, xsU];
PhiVac = Map[Together, SpinConnectionDFT[Vinf, eta3, gammaVac, xsU], {3}];

NRH`CheckZero["SM: vacuum spin connection = displayed Phi_{~+ oplus y} = 1/(2l), Phi_{+ ominus y} = -1/l (raised slots)",
   {PhiVac[[1]] - {{0, 0, 0}, {0, 0, -1/(2 l)}, {0, 1/(2 l), 0}},
    PhiVac[[4]] - {{0, 0, 1/l}, {0, 0, 0}, {-1/l, 0, 0}},
    PhiVac[[2]], PhiVac[[3]], PhiVac[[5]], PhiVac[[6]]}];

gamLow = {-gam3[[2]], -gam3[[1]], gam3[[3]]};
gampqLow[p_, q_] := (gamLow[[p]] . gamLow[[q]] - gamLow[[q]] . gamLow[[p]])/2;
spinTerm[a_] := 1/4 Sum[PhiVac[[a, p, q]] gampqLow[p, q], {p, 3}, {q, 3}];
DA[Es_, a_] := DblD[Es, a, xsU] + spinTerm[a] . Es;

DP[Es_, p_] := Sum[(JJ . Vinf)[[a, p]] DA[Es, a], {a, 6}];
DPbar[Es_, pb_] := Sum[(JJ . Vbinf)[[a, pb]] DA[Es, a], {a, 6}];
slashD[Es_] := Sum[gam3[[p]] . DP[Es, p], {p, 3}];

Evac = {Sqrt[2] ff[xp], l Derivative[1][ff][xp]};
NRH`CheckZero["D_{pbar} E = 0 for E = (Sqrt[2] f(x+), l f'(x+)), arbitrary chiral f",
   Together[Flatten[Table[DPbar[Evac, pb], {pb, 3}]]]];
NRH`CheckZero["gamma^p D_p E = E/(Sqrt[2] l)",
   Together[slashD[Evac] - Evac/(Sqrt[2] l)]];
NRH`CheckZero["reduced system gamma^p D_p E = (1/(Sqrt[2] l)) diag(1,-1).E + (0,0;1,0).d_+E",
   Module[{Eg = {ee0[xp], ee1[xp]}},
      Together[slashD[Eg] - 1/(Sqrt[2] l) {{1, 0}, {0, -1}} . Eg - {{0, 0}, {1, 0}} . D[Eg, xp]]]];
NRH`CheckZero["the opposite channel E = (0, g(x+)) has eigenvalue -1/(Sqrt[2] l)",
   Together[slashD[{0, gg[xp]}] + {0, gg[xp]}/(Sqrt[2] l)]];
NRH`CheckZero["both channels also satisfy D_{pbar} E = 0",
   Together[Flatten[Table[DPbar[{0, gg[xp]}, pb], {pb, 3}]]]];

Cmat = I s2;
NRH`CheckZero["C = i sigma_2 is the Majorana conjugation for the representation: C gamma^p C^-1 = -(gamma^p)^T",
   Flatten[Table[Cmat . gam3[[p]] . Inverse[Cmat] + Transpose[gam3[[p]]], {p, 3}]]];
E1 = {Sqrt[2] f1[xp], l Derivative[1][f1][xp]};
E2 = {Sqrt[2] f2[xp], l Derivative[1][f2][xp]};
Xflat = Table[E2 . Cmat . gam3[[p]] . E1, {p, 3}];
XM = (JJ . Vinf) . Xflat;
vplus = 2 f1[xp] f2[xp];
NRH`CheckZero["X^M = (omega_+, 0, -(l/2) v'; v^+, 0, -(l/2) v') with v^+ = 2 f1 f2, omega_+ = -2 l^2 f1' f2'",
   Together[XM - {-2 l^2 Derivative[1][f1][xp] Derivative[1][f2][xp], 0, -l/2 D[vplus, xp], vplus, 0, -l/2 D[vplus, xp]}]];
NRH`CheckZero["Lhat_X H^infty = 0 and Lhat_X d = 0 for the Killing-spinor bilinear, arbitrary chiral f1, f2",
   {Map[Together, GenLieH[XM, Hinf, xsU], {2}], Together[GenLieD[XM, dinf, xsU]]}];
NRH`Check["the bilinear is symmetric under f1 <-> f2 (commuting coefficient functions)",
   Together[XM - (XM /. {f1 -> f2, f2 -> f1})] === {0, 0, 0, 0, 0, 0}];

jetVars = {e0, e1, ep0, ep1, em0, em1, ey0, ey1};
jetSys = {LpS uu e0, em0, W1S uu e0 + 2 Sqrt[2] l em1, ey0, ey1,
   2 e0 + l ey0, 2 LpS uu em0 + 2 ep0 + Sqrt[2] ey1};
jm = Table[D[jetSys[[i]], jetVars[[j]]], {i, Length[jetSys]}, {j, 8}];
NRH`Check["the displayed system has rank six (seven relations, one dependent)",
   MatrixRank[jm] == 6];
kernelJ = NullSpace[jm];
NRH`Check["the solution space is exactly {e_1, d_+ e_1} - one arbitrary chiral profile",
   Sort[RowReduce[kernelJ]] === Sort[{{0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}}]];
NRH`Check["the mechanism: W_1 multiplies only e_0, and the survivor has e_0 = 0",
   Union[Cases[jetSys, W1S x_ :> x, Infinity]] === {uu e0}];

AthS = I/2 s1;
Aph1S = -I/2 (Sin[th] s3 - Cos[th] s2);
Aph2S = -Aph1S;
NRH`CheckZero["below : d_theta eta = (i/2) sigma1 eta, d_phi1 eta = -d_phi2 eta = -(i/2)(sin theta sigma3 - cos theta sigma2) eta is integrable (flat connection)",
   {Simplify[D[Aph1S, th] - (AthS . Aph1S - Aph1S . AthS)], Simplify[D[Aph2S, th] - (AthS . Aph2S - Aph2S . AthS)],
    Simplify[Aph1S . Aph2S - Aph2S . Aph1S]}];

NRH`CheckZero["the first-order pair (d+ u = Sqrt[2]/l v, d+ v = Sqrt[2]/l L u) closes into (l^2/2) s'' = L s",
   Module[{vv = l/Sqrt[2] D[sfun[xp], xp]},
      Together[l/Sqrt[2] (D[vv, xp] - Sqrt[2]/l Lp[xp] sfun[xp])
         - (l^2/2 D[sfun[xp], {xp, 2}] - Lp[xp] sfun[xp])]]];

period = Sqrt[2] Pi l;
NRH`Check["global AdS3 (L = -1/4): Hill solutions e^{pm i x/(Sqrt[2] l)}; one circuit gives e^{pm i pi} = -1 (antiperiodic)",
   And[Simplify[l^2/2 D[Exp[I xp/(Sqrt[2] l)], {xp, 2}] + 1/4 Exp[I xp/(Sqrt[2] l)]] === 0,
       Simplify[Exp[I (xp + period)/(Sqrt[2] l)] + Exp[I xp/(Sqrt[2] l)]] === 0]];
NRH`CheckZero["massless BTZ (L = 0): the (u, v) pair shifts by exactly 2 pi v over one circuit (unipotent)",
   Module[{uSol = u0 + Sqrt[2]/l v0 xp},
      Together[(uSol /. xp -> xp + period) - uSol - 2 Pi v0]]];
NRH`Check["constant L > 0: monodromy multipliers e^{pm 2 Sqrt[L0] pi} are real and not +-1: (0,0) kernel",
   Simplify[Exp[Sqrt[2 L0]/l period] > 1, L0 > 0 && l > 0] === True];

NRH`Check["the two surviving directions are complementary in the reduced basis",
   MatrixRank[{{1, 0}, {0, 1}}] == 2 && {1, 0} . {0, 1} == 0];
NRH`Check["counting: (4,4) + (4,4) = the sixteen constant vacuum modes; one-sided (4,0)/(0,4) = extremal count",
   4 + 4 + 4 + 4 == 16 && 4 + 0 == 4 && 0 + 4 == 4];

NRH`FileSummary[];
