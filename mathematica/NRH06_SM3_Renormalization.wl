(* ::Title:: *)
(*NRH06 SM3 Radial Momenta, Counterterms and Two-Point Functions*)

(* Current SM3.2: action variations and counterterms; SM3.3: the vacuum normalization;
   SM3.4-3.5: general R/NR responses and kernels. Expanded ordered action calculations and
   contact rows are retained as ancillary checks, with stable pre-rewrite assertion identifiers.
   Both variation slots are extended to near-boundary bulk solutions. These checks provide
   indirect support for the compact variation formulas, not a direct equality test of the
   new curved-index first variation or symmetrically polarized Hessian.
   The fixed-flux K3 master equation and interior matching of H_R are not tested here.
   Historical conditional uniformization checks establish algebraic factors only; they do
   not derive the vacuum hair response or prescribe the general NR interior state.
   z = e^{-2y/l} is the manuscript's u; yy is the explicit y; the cutoff is y = Y = yy, e^{2Y/l} = 1/z. *)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH06_SM3_Renormalization.wl"];
Needs["VariationalMethods`"];
JJ = ODDJ[3];

(* ::Section:: *)
(* Current SM3.2: radial momenta and ancillary ordered second-variation checks *)

(* GammaDFT, variation, defB: the Gamma^2 variation identity with generic tangential h and delta d on constant BTZ *)
xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};
w2u = Sqrt[u];
E3btz = {{w2u, -Lm0/w2u, 0}, {-Lp0/w2u, w2u, 0}, {0, 0, 1}};
g3btz = Map[Together, Transpose[E3btz] . eta3 . E3btz, {2}];
B3btz = (u + Lp0 Lm0/u) {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
Eibtz = Map[Together, Inverse[E3btz], {2}];
Vbtz = Map[Together, 1/Sqrt[2] ArrayFlatten[{{Eibtz}, {(g3btz + B3btz) . Eibtz}}], {2}];
Vbbtz = Map[Together, 1/Sqrt[2] ArrayFlatten[{{Eibtz}, {(B3btz - g3btz) . Eibtz}}], {2}];
Hbtz = Map[Together, RiemannianH[g3btz, B3btz], {2}];
dbtz = -1/2 Log[u (1 - Lp0 Lm0/u^2)];
hmatV = {{hpp[xp, xm, u], hpm[xp, xm, u], 0}, {hmp[xp, xm, u], hmm[xp, xm, u], 0}, {0, 0, 0}};
deltaHV = Map[Together, MixedFluctuationH[Vbtz, Vbbtz, hmatV], {2}];
ddV = dd[xp, xm, u];
gammaTV = GammaDFT[Hbtz + t deltaHV, dbtz + t ddV, xsU];
gamma0V = Map[Together, gammaTV /. t -> 0, {3}];
gamma1V = Map[Together, D[gammaTV, t] /. t -> 0, {3}];
LHS = Together[D[Gamma2Density[Hbtz + t deltaHV, dbtz + t ddV, gamma0V + t gamma1V, xsU], t] /. t -> 0];
r4V = RiemannR4[gamma0V, xsU]; ricV = RicciS[gamma0V, r4V, xsU]; s0V = Together[ScalarS0[Hbtz, dbtz, xsU]];
VupV = JJ . Vbtz; VbupV = JJ . Vbbtz;
SpqV = Map[Together, Transpose[VupV] . ricV . VbupV, {2}];
hupV = eta3 . hmatV . etab3;
bulk = 2 Exp[-2 dbtz] (Sum[hupV[[p, q]] SpqV[[p, q]], {p, 3}, {q, 3}] - ddV s0V);
trvecV = Table[Sum[JJ[[L, A]] gamma0V[[A, L, N]], {L, 6}, {A, 6}], {N, 6}];
trpV = Table[Sum[trvecV[[M]] VupV[[M, p]], {M, 6}], {p, 3}];
trqV = Table[Sum[trvecV[[N]] VbupV[[N, q]], {N, 6}], {q, 3}];
GpqV = Table[Sum[VupV[[M, p]] JJ[[K, A]] gamma0V[[M, A, N]] VbupV[[N, q]], {M, 6}, {A, 6}, {N, 6}], {K, 6}, {p, 3}, {q, 3}];
GqpV = Table[Sum[VbupV[[N, q]] JJ[[K, A]] gamma0V[[N, A, M]] VupV[[M, p]], {M, 6}, {A, 6}, {N, 6}], {K, 6}, {p, 3}, {q, 3}];
AKV = Table[VupV[[K, p]] trqV[[q]] + VbupV[[K, q]] trpV[[p]] - GpqV[[K, p, q]] - GqpV[[K, p, q]], {K, 6}, {p, 3}, {q, 3}];
BKV = GammaBVector[Hbtz, dbtz, xsU];
flux = Sum[DblD[Exp[-2 dbtz] Sum[hupV[[p, q]] AKV[[K, p, q]], {p, 3}, {q, 3}] + 2 ddV Exp[-2 dbtz] BKV[[K]], K, xsU], {K, 6}];
NRH`CheckZero["variation, SMBdefinition, SMAdefinition: delta L_Gamma2 = 2e^{-2d}(h^{p qbar} S_{p qbar} - delta d S_(0)) + d_K(h e^{-2d} A^K + 2 delta d e^{-2d} B^K), generic tangential h and delta d on BTZ",
   Together[LHS - bulk - flux]];
NRH`Check["the identity is not vacuous: the flux carries the fluctuations, S_{p qbar} = 0 and S_(0) = -4/l^2 on shell",
   ! FreeQ[flux, hpp] && ! FreeQ[flux, dd] && Together[s0V + 4/l^2] === 0 && NRH`ZeroQ[SpqV]];
NRH`CheckZero["SMrenvariation, SMvolumecounterterm: on shell the bulk term cancels delta(-2 Lambda e^{-2d}) and the counterterm shifts B^y by 4/l",
   {Together[2 (-dd0) (-4/l^2) + 4 (-2/l^2) dd0], Together[-(4/l) (-2 dd0) - 2 dd0 (4/l)]}];
NRH`CheckZero["SMmetricradialdictionary, SMscalarradialdictionary: coefficient matching -2K = (16 pi G)^{-1} A^y and 2 T_(0) = (16 pi G)^{-1} 2 (B^y + 4/l)",
   {Together[-2 (-(1/(32 Pi G))) - 1/(16 Pi G)], Together[2 (1/(16 Pi G)) - 2/(16 Pi G)]}];

(* SMsecondvariation: independent flat arguments; delta_2 h_1^{p qbar} = 0 in local coset coordinates *)
Module[{j, s, h1, h2, generator, t1, t2, h12, dh1, dh2, cross, v, vb},
   j = ArrayFlatten[{{eta3, 0}, {0, etab3}}];
   s = ArrayFlatten[{{eta3, 0}, {0, -etab3}}];
   h1 = Array[a, {3, 3}]; h2 = Array[b, {3, 3}];
   generator[h_] := ArrayFlatten[{{ConstantArray[0, {3, 3}], -eta3 . h/2}, {etab3 . Transpose[h]/2, ConstantArray[0, {3, 3}]}}];
   t1 = generator[h1]; t2 = generator[h2]; cross = (t1 . t2 + t2 . t1)/2;
   dh1 = t1 . s + s . Transpose[t1]; dh2 = t2 . s + s . Transpose[t2];
   h12 = cross . s + t1 . s . Transpose[t2] + t2 . s . Transpose[t1] + s . Transpose[cross];
   v = IdentityMatrix[6][[All, 1 ;; 3]]; vb = IdentityMatrix[6][[All, 4 ;; 6]];
   NRH`CheckZero["SMsecondvariation: mixed coset constraint at second order", h12 . j . s + s . j . h12 + dh1 . j . dh2 + dh2 . j . dh1];
   NRH`CheckZero["SMsecondvariation: h_1 projection and delta_2 h_1^{p qbar} = 0",
      {Transpose[v] . j . dh1 . j . vb - h1,
       Transpose[t2 . v] . j . dh1 . j . vb + Transpose[v] . j . h12 . j . vb + Transpose[v] . j . dh1 . j . (t2 . vb)}];
   NRH`Check["negative control: independent h does not mean delta_2 delta_1 H = 0", ! NRH`ZeroQ[h12]]];

(* general z-series saddles, fluctuations of (xp, xm, yy) *)
sdR = SaddleSeries[RiemannianSaddleExact[], 1];
sdNR = SaddleSeries[NonRiemannianSaddleExact[W1[xp, xm] z + NRW2 z^2], 1];
hmat = {{hpp[xp, xm, yy], hpm[xp, xm, yy], 0}, {hmp[xp, xm, yy], hmm[xp, xm, yy], 0}, {0, 0, 0}};
ddf = dd[xp, xm, yy];
momentumData[bg_] := Module[{dH, Hlin, dlin, gamma, gamma0, gamma1, core0, core1, Vup, Vbup, dV, dVb, Aproj, dAmoving, dAfrozen, dAframes, Bv, dB, Aun, Abar},
   dH = SeriesZM[MixedFluctuationH[bg["V"], bg["Vb"], hmat], 1];
   Hlin = bg["H"] + t dH; dlin = bg["d"] + t ddf;
   gamma = GammaDFT[Hlin, dlin, xsZY];
   gamma0 = Map[SeriesZ[# /. t -> 0, 1] &, gamma, {3}];
   gamma1 = Map[SeriesZ[LinearT[#], 1] &, gamma, {3}];
   core0 = MomentumCore[gamma0, xsZY][[6]]; core1 = MomentumCore[gamma1, xsZY][[6]];
   Vup = JJ . bg["V"]; Vbup = JJ . bg["Vb"];
   {dV, dVb} = FrameVariation[bg["V"], bg["Vb"], hmat];
   Aproj = SeriesZM[Transpose[Vup] . core0 . Vbup, 1];
   dAfrozen = SeriesZM[Transpose[Vup] . core1 . Vbup, 1];
   dAframes = SeriesZM[Transpose[JJ . dV] . core0 . Vbup + Transpose[Vup] . core0 . (JJ . dVb), 1];
   dAmoving = dAfrozen + dAframes;
   Aun = SeriesZM[Transpose[Vup] . core0 . Vup, 1]; Abar = SeriesZM[Transpose[Vbup] . core0 . Vbup, 1];
   Bv = SeriesZ[GammaBVector[bg["H"], bg["d"], xsZY][[6]], 1];
   dB = SeriesZ[LinearT[GammaBVector[Hlin, dlin, xsZY][[6]]], 1];
   <|"A" -> Aproj, "dA" -> dAmoving, "dAfrozen" -> dAfrozen, "dAframes" -> dAframes, "Aun" -> Aun, "Abar" -> Abar, "B" -> Bv, "dB" -> dB|>];
Print["  computing background momenta and their variations on both saddles ..."];
{tmR, momR} = AbsoluteTiming[momentumData[sdR]]; Print["  [timing] R: ", Round[tmR, 0.1], " s"];
{tmNR, momNR} = AbsoluteTiming[momentumData[sdNR]]; Print["  [timing] NR: ", Round[tmNR, 0.1], " s"];
Do[
   Module[{m = data[[2]], s = data[[1]], target, Jsym},
      target = If[s == "R", {{Lp[xp], Lp[xp] Lm[xm]}, {1, Lm[xm]}}, {{1/psip[xp]^2, W1[xp, xm]/4}, {0, 1/psim[xm]^2}}];
      NRH`CheckZero["SMbackgroundmomenta on " <> s <> ": A^y_{a bbar} = -(2u/l) {{L+, J_s},{eps_s, L-}} + O(u^2), B^y = -4/l + O(u^2)",
         {Together[m["A"][[1 ;; 2, 1 ;; 2]] + (2 z/l) target], Together[m["B"] + 4/l]}];
      NRH`CheckZero["SMsamechiralitymomenta on " <> s <> ": tangential A^y_{ab} = eta/l, A^y_{abar bbar} = -etabar/l + O(u^2)",
         {Together[m["Aun"][[1 ;; 2, 1 ;; 2]] - eta3[[1 ;; 2, 1 ;; 2]]/l], Together[m["Abar"][[1 ;; 2, 1 ;; 2]] + etab3[[1 ;; 2, 1 ;; 2]]/l]}];
      NRH`CheckZero["SMframemomentumcheck on " <> s <> ": the two frame terms give -(1/l) h_{a bbar} + O(u^2)",
         Together[m["dAframes"][[1 ;; 2, 1 ;; 2]] + hmat[[1 ;; 2, 1 ;; 2]]/l]];
      NRH`CheckZero["SMconnectionmomentumcheck on " <> s <> ": V^M_a Vbar^N_bbar delta A^y_MN = (1/2) d_y h_{a bbar} + h_{a bbar}/l + O(u^2)",
         Together[m["dAfrozen"][[1 ;; 2, 1 ;; 2]] - Map[DyZY, hmat[[1 ;; 2, 1 ;; 2]]]/2 - hmat[[1 ;; 2, 1 ;; 2]]/l]];
      NRH`CheckZero["SMlinearizedmomenta on " <> s <> ": delta A^y_{a bbar} = (1/2) d_y h_{a bbar} + O(u^2), delta B^y = 4 d_y delta d + O(u^2)",
         {Together[m["dA"][[1 ;; 2, 1 ;; 2]] - Map[DyZY, hmat[[1 ;; 2, 1 ;; 2]]]/2], Together[m["dB"] - 4 DyZY[ddf]]}];
      NRH`Check["negative control on " <> s <> ": freezing the bulk frames changes the tangential momentum variation",
         ! NRH`ZeroQ[Together[m["dAfrozen"][[1 ;; 2, 1 ;; 2]] - m["dA"][[1 ;; 2, 1 ;; 2]]]]];
      (* SMonept / SMNRonept: the one-point matrices from the radial dictionary *)
      NRH`CheckZero["SMonept/SMNRonept on " <> s <> ": <K_{a bbar}> = -(32 pi G)^{-1} FP e^{2Y/l} A^y = target/(16 pi G l), <T_(0)> = 0",
         {Together[-(1/(32 Pi G)) Coefficient[m["A"][[1 ;; 2, 1 ;; 2]], z, 1] - target/(16 Pi G l)], Together[Coefficient[m["B"] + 4/l, z, 1]]}]],
   {data, {{"R", momR}, {"NR", momNR}}}];

(* SMstatefluctuations: source-free tangent variations of the two exact families projected on the saddle frames *)
Module[{rEx = RiemannianSaddleExact[], dHR, hR, nrF, dHNR, hNR, dW1},
   dHR = D[rEx["H"] /. {Lp -> Function[x, Lp[x] + s dLp[x]], Lm -> Function[x, Lm[x] + s dLm[x]]}, s] /. s -> 0;
   hR = SeriesZM[Transpose[JJ . rEx["V"]] . dHR . (JJ . rEx["Vb"]), 1];
   NRH`CheckZero["SMstatefluctuations on R: h_{a bbar} = u {{2 dL+, 2 d(L+L-)},{0, 2 dL-}} + O(u^2), delta d = O(u^2)",
      {Together[hR[[1 ;; 2, 1 ;; 2]] - z {{2 dLp[xp], 2 (dLp[xp] Lm[xm] + Lp[xp] dLm[xm])}, {0, 2 dLm[xm]}}],
       SeriesZ[D[rEx["d"] /. {Lp -> Function[x, Lp[x] + s dLp[x]]}, s] /. s -> 0, 1]}];
   nrF = NonRiemannianSaddleExact[(W1[xp, xm] + s dW1[xp, xm]) z + NRW2 z^2];
   dHNR = D[nrF["H"] /. {psip -> Function[x, psip[x] + s dpsip[x]], psim -> Function[x, psim[x] + s dpsim[x]]}, s] /. s -> 0;
   hNR = SeriesZM[Transpose[JJ . nrF["V"]] . dHNR . (JJ . nrF["Vb"]) /. s -> 0, 1];
   NRH`CheckZero["SMstatefluctuations on NR: h_{a bbar} = u {{2 dL+, dW1/2},{0, 2 dL-}} + O(u^2), delta d = O(u^2)",
      {Together[hNR[[1 ;; 2, 1 ;; 2]] - z {{2 D[1/psip[xp]^2, psip[xp]] dpsip[xp], dW1[xp, xm]/2}, {0, 2 D[1/psim[xm]^2, psim[xm]] dpsim[xm]}}],
       SeriesZ[D[nrF["d"] /. {psip -> Function[x, psip[x] + s dpsip[x]]}, s] /. s -> 0, 1]}]];

(* ::Section:: *)
(* Current SM3.2-3.5: ordered bulk-solution Hessian, counterterms, ancillary finite rows, and kernels *)

(* Expanded near-boundary solutions checked in NRH05 and underlying current SM3.4-3.5.
   Slot k = 1, 2 carries sources a_k, b_k, r_k, c_k, v_k
   = (h^(0)_{om omb}, h^(0)_{op opb}, h^(0)_{om opb}, h^(0)_{op omb}, delta d^(0)).
   Rp_k,Rm_k are the full h_pp^(2,0),h_mm^(2,0), including compact b/omega/local terms;
   H_k is the full hair response H_s, not a partial Omega coefficient. On R, cs_k = -4 f_0. *)
(* Historical machine-transcribed expanded targets; the current compact a/b/omega/U formulas require a separate direct bridge. *)
(* Fields: hpp,hpm,hmp,hmm,dd of [xp,xm,yy]; sources a0,b0,r0,c0,v0 of [xp,xm]; responses Rp,Rm of [xp,xm]. *)
E0target["pp"] := ReleaseHold[-1/4*Hold[D[hmp[xp, xm, yy], {xp, 2}]] - 1/4*Hold[D[hpp[xp, xm, yy], {yy, 2}]] - 1/2*Hold[D[hpp[xp, xm, yy], yy]]/l];
E0target["pm"] := ReleaseHold[-1/4*Hold[D[hmm[xp, xm, yy], {xp, 2}]] - 1/4*Hold[D[hpm[xp, xm, yy], {yy, 2}]] - 1/4*Hold[D[hpp[xp, xm, yy], {xm, 2}]] + Hold[D[dd[xp, xm, yy], xm, xp]] - 1/2*Hold[D[hpm[xp, xm, yy], yy]]/l];
E0target["py"] := ReleaseHold[Hold[D[dd[xp, xm, yy], xp, yy]] - 1/4*Hold[D[hpp[xp, xm, yy], xm, yy]]];
E0target["mp"] := ReleaseHold[-1/4*Hold[D[hmp[xp, xm, yy], {yy, 2}]] - 1/2*Hold[D[hmp[xp, xm, yy], yy]]/l];
E0target["mm"] := ReleaseHold[-1/4*Hold[D[hmm[xp, xm, yy], {yy, 2}]] - 1/4*Hold[D[hmp[xp, xm, yy], {xm, 2}]] - 1/2*Hold[D[hmm[xp, xm, yy], yy]]/l];
E0target["my"] := ReleaseHold[-1/4*Hold[D[hmp[xp, xm, yy], xm, yy]]];
E0target["yp"] := ReleaseHold[-1/4*Hold[D[hmp[xp, xm, yy], xp, yy]]];
E0target["ym"] := ReleaseHold[Hold[D[dd[xp, xm, yy], xm, yy]] - 1/4*Hold[D[hmm[xp, xm, yy], xp, yy]]];
E0target["yy"] := ReleaseHold[Hold[D[dd[xp, xm, yy], {yy, 2}]]];
E0target["s0"] := ReleaseHold[4*Hold[D[dd[xp, xm, yy], {yy, 2}]] + Hold[D[hmp[xp, xm, yy], xm, xp]] + 8*Hold[D[dd[xp, xm, yy], yy]]/l];
E1NRtarget["pp"] := ReleaseHold[-1/2*Lp[xp]*Hold[D[hmp[xp, xm, yy], xm, xp]] - 1/4*Hold[D[Lp[xp], xp]]*Hold[D[hmp[xp, xm, yy], xm]] - 2*Lp[xp]*Hold[D[dd[xp, xm, yy], yy]]/l];
E1NRtarget["pm"] := ReleaseHold[Lm[xm]*Hold[D[dd[xp, xm, yy], {xp, 2}]] - 1/2*Lm[xm]*Hold[D[hpp[xp, xm, yy], xm, xp]] + Lp[xp]*Hold[D[dd[xp, xm, yy], {xm, 2}]] - 1/2*Lp[xp]*Hold[D[hmm[xp, xm, yy], xm, xp]] - 1/4*Hold[D[Lm[xm], xm]]*Hold[D[hpp[xp, xm, yy], xp]] - 1/4*Hold[D[Lp[xp], xp]]*Hold[D[hmm[xp, xm, yy], xm]] - 1/2*W1[xp, xm]*Hold[D[dd[xp, xm, yy], yy]]/l];
E1NRtarget["py"] := ReleaseHold[-1/4*Lm[xm]*Hold[D[hpp[xp, xm, yy], xp, yy]] + Lp[xp]*Hold[D[dd[xp, xm, yy], xm, yy]] - 2*Lp[xp]*Hold[D[dd[xp, xm, yy], xm]]/l + Lp[xp]*Hold[D[hmm[xp, xm, yy], xp]]/l + (1/4)*W1[xp, xm]*Hold[D[hmp[xp, xm, yy], xp]]/l + (1/2)*hmm[xp, xm, yy]*Hold[D[Lp[xp], xp]]/l + (1/8)*hmp[xp, xm, yy]*Hold[D[W1[xp, xm], xp]]/l];
E1NRtarget["mp"] := ReleaseHold[0];
E1NRtarget["mm"] := ReleaseHold[-1/2*Lm[xm]*Hold[D[hmp[xp, xm, yy], xm, xp]] - 1/4*Hold[D[Lm[xm], xm]]*Hold[D[hmp[xp, xm, yy], xp]] - 2*Lm[xm]*Hold[D[dd[xp, xm, yy], yy]]/l];
E1NRtarget["my"] := ReleaseHold[-1/4*Lm[xm]*Hold[D[hmp[xp, xm, yy], xp, yy]]];
E1NRtarget["yp"] := ReleaseHold[-1/4*Lp[xp]*Hold[D[hmp[xp, xm, yy], xm, yy]]];
E1NRtarget["ym"] := ReleaseHold[Lm[xm]*Hold[D[dd[xp, xm, yy], xp, yy]] - 1/4*Lp[xp]*Hold[D[hmm[xp, xm, yy], xm, yy]] - 2*Lm[xm]*Hold[D[dd[xp, xm, yy], xp]]/l + Lm[xm]*Hold[D[hpp[xp, xm, yy], xm]]/l + (1/4)*W1[xp, xm]*Hold[D[hmp[xp, xm, yy], xm]]/l + (1/8)*hmp[xp, xm, yy]*Hold[D[W1[xp, xm], xm]]/l + (1/2)*hpp[xp, xm, yy]*Hold[D[Lm[xm], xm]]/l];
E1NRtarget["yy"] := ReleaseHold[Lm[xm]*Hold[D[hpp[xp, xm, yy], yy]]/l + Lp[xp]*Hold[D[hmm[xp, xm, yy], yy]]/l + (1/4)*W1[xp, xm]*Hold[D[hmp[xp, xm, yy], yy]]/l];
E1NRtarget["s0"] := ReleaseHold[Lm[xm]*Hold[D[hmp[xp, xm, yy], {xp, 2}]] + Lp[xp]*Hold[D[hmp[xp, xm, yy], {xm, 2}]] + 2*Lm[xm]*Hold[D[hpp[xp, xm, yy], yy]]/l + 2*Lp[xp]*Hold[D[hmm[xp, xm, yy], yy]]/l + (1/2)*W1[xp, xm]*Hold[D[hmp[xp, xm, yy], yy]]/l];
DE1target["pp"] := ReleaseHold[Hold[D[dd[xp, xm, yy], {xp, 2}]]];
DE1target["pm"] := ReleaseHold[0];
DE1target["py"] := ReleaseHold[-1/4*Hold[D[hpm[xp, xm, yy], xp, yy]]];
DE1target["mp"] := ReleaseHold[-2*Hold[D[dd[xp, xm, yy], yy]]/l];
DE1target["mm"] := ReleaseHold[Hold[D[dd[xp, xm, yy], {xm, 2}]]];
DE1target["my"] := ReleaseHold[Hold[D[dd[xp, xm, yy], xm, yy]] - 1/4*Hold[D[hmm[xp, xm, yy], xp, yy]] - 2*Hold[D[dd[xp, xm, yy], xm]]/l];
DE1target["yp"] := ReleaseHold[Hold[D[dd[xp, xm, yy], xp, yy]] - 1/4*Hold[D[hpp[xp, xm, yy], xm, yy]] - 2*Hold[D[dd[xp, xm, yy], xp]]/l];
DE1target["ym"] := ReleaseHold[-1/4*Hold[D[hpm[xp, xm, yy], xm, yy]]];
DE1target["yy"] := ReleaseHold[Hold[D[hpm[xp, xm, yy], yy]]/l];
DE1target["s0"] := ReleaseHold[Hold[D[hmm[xp, xm, yy], {xp, 2}]] + Hold[D[hpp[xp, xm, yy], {xm, 2}]] - 8*Hold[D[dd[xp, xm, yy], xm, xp]] + 2*Hold[D[hpm[xp, xm, yy], yy]]/l];
h2Ltarget["R", "pp"] := ReleaseHold[-1/8*l^3*Hold[D[r0[xp, xm], xm, {xp, 3}]] + (1/2)*l*Lp[xp]*Hold[D[r0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], xm]]];
h2Ltarget["R", "pm"] := ReleaseHold[(1/16)*l^5*Hold[D[r0[xp, xm], {xm, 3}, {xp, 3}]] - 1/4*l^3*Lm[xm]*Hold[D[r0[xp, xm], xm, {xp, 3}]] - 1/4*l^3*Lp[xp]*Hold[D[r0[xp, xm], {xm, 3}, xp]] - 1/4*l^3*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], {xp, 3}]] - 1/4*l^3*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], {xm, 3}]] - 1/4*l^3*Hold[D[a0[xp, xm], xm, {xp, 3}]] - 1/4*l^3*Hold[D[b0[xp, xm], {xm, 3}, xp]] + l^3*Hold[D[v0[xp, xm], {xm, 2}, {xp, 2}]] - 1/2*l*Lm[xm]*Lp[xp]*Hold[D[r0[xp, xm], xm, xp]] - 2*l*Lm[xm]*Hold[D[v0[xp, xm], {xp, 2}]] + l*Lm[xm]*Hold[D[b0[xp, xm], xm, xp]] - 2*l*Lp[xp]*Hold[D[v0[xp, xm], {xm, 2}]] + l*Lp[xp]*Hold[D[a0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lm[xm], xm]]*Hold[D[b0[xp, xm], xp]] + (1/2)*l*Hold[D[Lp[xp], xp]]*Hold[D[a0[xp, xm], xm]] + (1/2)*l*Hold[D[Rm[xp, xm], {xp, 2}]] + (1/2)*l*Hold[D[Rp[xp, xm], {xm, 2}]]];
h2Ltarget["R", "mp"] := ReleaseHold[-1/2*l*Hold[D[r0[xp, xm], xm, xp]]];
h2Ltarget["R", "mm"] := ReleaseHold[-1/8*l^3*Hold[D[r0[xp, xm], {xm, 3}, xp]] + (1/2)*l*Lm[xm]*Hold[D[r0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], xp]]];
dd2Ltarget["R"] := ReleaseHold[-1/16*l^3*Hold[D[r0[xp, xm], {xm, 2}, {xp, 2}]]];
dd2target["R"] := ReleaseHold[-1/32*l^4*Hold[D[r0[xp, xm], {xm, 2}, {xp, 2}]] + (1/8)*l^2*Lm[xm]*Hold[D[r0[xp, xm], {xp, 2}]] + (1/8)*l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 2}]] + (1/8)*l^2*Hold[D[a0[xp, xm], {xp, 2}]] + (1/8)*l^2*Hold[D[b0[xp, xm], {xm, 2}]] - 1/2*l^2*Hold[D[v0[xp, xm], xm, xp]]];
Fptarget["R"] := ReleaseHold[-1/8*l^4*Hold[D[r0[xp, xm], {xm, 2}, {xp, 3}]] + (1/4)*l^2*Lm[xm]*Hold[D[r0[xp, xm], {xp, 3}]] + l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 2}, xp]] + (3/4)*l^2*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], {xm, 2}]] + (1/4)*l^2*Hold[D[a0[xp, xm], {xp, 3}]] + (1/4)*l^2*Hold[D[b0[xp, xm], {xm, 2}, xp]] - l^2*Hold[D[v0[xp, xm], xm, {xp, 2}]] - 2*Lm[xm]*Lp[xp]*Hold[D[r0[xp, xm], xp]] - Lm[xm]*r0[xp, xm]*Hold[D[Lp[xp], xp]] - 2*Lp[xp]*Hold[D[a0[xp, xm], xp]] + 4*Lp[xp]*Hold[D[v0[xp, xm], xm]] - a0[xp, xm]*Hold[D[Lp[xp], xp]]];
Fmtarget["R"] := ReleaseHold[-1/8*l^4*Hold[D[r0[xp, xm], {xm, 3}, {xp, 2}]] + l^2*Lm[xm]*Hold[D[r0[xp, xm], xm, {xp, 2}]] + (1/4)*l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 3}]] + (3/4)*l^2*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], {xp, 2}]] + (1/4)*l^2*Hold[D[b0[xp, xm], {xm, 3}]] + (1/4)*l^2*Hold[D[a0[xp, xm], xm, {xp, 2}]] - l^2*Hold[D[v0[xp, xm], {xm, 2}, xp]] - 2*Lm[xm]*Lp[xp]*Hold[D[r0[xp, xm], xm]] - 2*Lm[xm]*Hold[D[b0[xp, xm], xm]] + 4*Lm[xm]*Hold[D[v0[xp, xm], xp]] - Lp[xp]*r0[xp, xm]*Hold[D[Lm[xm], xm]] - b0[xp, xm]*Hold[D[Lm[xm], xm]]];
h2Ltarget["NR", "pp"] := ReleaseHold[(1/2)*l*Lp[xp]*Hold[D[r0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], xm]]];
h2Ltarget["NR", "pm"] := ReleaseHold[-1/4*l^3*Lm[xm]*Hold[D[r0[xp, xm], xm, {xp, 3}]] - 1/4*l^3*Lp[xp]*Hold[D[r0[xp, xm], {xm, 3}, xp]] - 1/4*l^3*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], {xp, 3}]] - 1/4*l^3*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], {xm, 3}]] - 2*l*Lm[xm]*Hold[D[v0[xp, xm], {xp, 2}]] + l*Lm[xm]*Hold[D[b0[xp, xm], xm, xp]] - 2*l*Lp[xp]*Hold[D[v0[xp, xm], {xm, 2}]] + l*Lp[xp]*Hold[D[a0[xp, xm], xm, xp]] - 1/8*l*W1[xp, xm]*Hold[D[r0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lm[xm], xm]]*Hold[D[b0[xp, xm], xp]] + (1/2)*l*Hold[D[Lp[xp], xp]]*Hold[D[a0[xp, xm], xm]] + (1/2)*l*Hold[D[Rm[xp, xm], {xp, 2}]] + (1/2)*l*Hold[D[Rp[xp, xm], {xm, 2}]]];
h2Ltarget["NR", "mp"] := ReleaseHold[0];
h2Ltarget["NR", "mm"] := ReleaseHold[(1/2)*l*Lm[xm]*Hold[D[r0[xp, xm], xm, xp]] + (1/2)*l*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], xp]]];
dd2Ltarget["NR"] := ReleaseHold[0];
dd2target["NR"] := ReleaseHold[(1/8)*l^2*Lm[xm]*Hold[D[r0[xp, xm], {xp, 2}]] + (1/8)*l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 2}]]];
Fptarget["NR"] := ReleaseHold[(1/4)*l^2*Lm[xm]*Hold[D[r0[xp, xm], {xp, 3}]] + l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 2}, xp]] + (3/4)*l^2*Hold[D[Lp[xp], xp]]*Hold[D[r0[xp, xm], {xm, 2}]] - 2*Lp[xp]*Hold[D[a0[xp, xm], xp]] + 4*Lp[xp]*Hold[D[v0[xp, xm], xm]] - 1/2*W1[xp, xm]*Hold[D[r0[xp, xm], xp]] - a0[xp, xm]*Hold[D[Lp[xp], xp]] - 1/4*r0[xp, xm]*Hold[D[W1[xp, xm], xp]]];
Fmtarget["NR"] := ReleaseHold[l^2*Lm[xm]*Hold[D[r0[xp, xm], xm, {xp, 2}]] + (1/4)*l^2*Lp[xp]*Hold[D[r0[xp, xm], {xm, 3}]] + (3/4)*l^2*Hold[D[Lm[xm], xm]]*Hold[D[r0[xp, xm], {xp, 2}]] - 2*Lm[xm]*Hold[D[b0[xp, xm], xm]] + 4*Lm[xm]*Hold[D[v0[xp, xm], xp]] - 1/2*W1[xp, xm]*Hold[D[r0[xp, xm], xm]] - b0[xp, xm]*Hold[D[Lm[xm], xm]] - 1/4*r0[xp, xm]*Hold[D[W1[xp, xm], xm]]];
slotRules[k_] := {a0 -> Symbol["a" <> ToString[k]], b0 -> Symbol["b" <> ToString[k]], r0 -> Symbol["r" <> ToString[k]],
   c0 -> Symbol["c" <> ToString[k]], v0 -> Symbol["v" <> ToString[k]], Rp -> Symbol["Rp" <> ToString[k]], Rm -> Symbol["Rm" <> ToString[k]]};
srcK[k_] := <|"pp" -> Symbol["b" <> ToString[k]][xp, xm], "pm" -> Symbol["c" <> ToString[k]][xp, xm],
   "mp" -> Symbol["r" <> ToString[k]][xp, xm], "mm" -> Symbol["a" <> ToString[k]][xp, xm]|>;
respK[k_] := <|"pp" -> Symbol["Rp" <> ToString[k]][xp, xm], "mm" -> Symbol["Rm" <> ToString[k]][xp, xm],
   "pm" -> Symbol["H" <> ToString[k]][xp, xm]|>;
epsS[s_] := If[s == "R", 1, 0];
Jsrc[s_] := If[s == "R", Lp[xp] Lm[xm], W1[xp, xm]/4];
solution[s_, k_] := Module[{a, b, r, c, v, h1, h1b, C0, h2, h2L, dd1, f},
   {a, b, r, c, v} = Table[Symbol[n <> ToString[k]][xp, xm], {n, {"a", "b", "r", "c", "v"}}];
   dd1 = -(l/8) D[r, xp, xm];
   h1 = <|"mp" -> 0, "pp" -> -(l/2) D[r, {xp, 2}], "mm" -> -(l/2) D[r, {xm, 2}],
      "pm" -> (l/2) (4 D[v, xp, xm] - D[a, {xp, 2}] - D[b, {xm, 2}] - (l^2/4) D[r, {xp, 2}, {xm, 2}])|>;
   h1b = <|"mp" -> 0, "pp" -> 0, "mm" -> 0, "pm" -> (l^2/8) D[r, {xp, 2}, {xm, 2}]|>;
   C0 = Symbol["cs" <> ToString[k]] + If[s == "R", 4 v - (l^2/4) D[r, xp, xm], 0];
   h2 = <|"pp" -> respK[k]["pp"], "mm" -> respK[k]["mm"], "pm" -> respK[k]["pm"], "mp" -> C0|>;
   h2L = Association[Table[ch -> (h2Ltarget[s, ch] /. slotRules[k]), {ch, {"pp", "pm", "mp", "mm"}}]];
   f = Association[Table[ch -> srcK[k][ch] + yy h1[ch] + yy^2 h1b[ch] + z (h2[ch] + yy h2L[ch]), {ch, {"pp", "pm", "mp", "mm"}}]];
   f["dd"] = v + yy dd1 + z ((dd2target[s] /. slotRules[k]) + yy (dd2Ltarget[s] /. slotRules[k]));
   f];
(* SMfinitepartexplicit: u^{-1} delta A^y_{a bbar} = u^{-1} (1/2) d_y h_{a bbar} and u^{-1} delta B^y = 4 u^{-1} d_y delta d on the radial ansatz *)
Do[Module[{f = solution[s, 1], coefZ, coefY, lhsA, lhsB, chs = {"pp", "pm", "mp", "mm"}},
      coefZ[e_, n_] := Coefficient[Expand[e], z, n]; coefY[e_, j_] := Coefficient[Expand[e], yy, j];
      lhsA = Association[Table[ch -> Expand[(1/z) DyZY[f[ch]]/2], {ch, chs}]];
      lhsB = Expand[(1/z) 4 DyZY[f["dd"]]];
      NRH`CheckZero["SMfinitepartexplicit on " <> s <> ": u^{-1} delta A^y_{a bbar} = u^{-1}(h^(1)/2 + y h^(1b)) - h^(2)/l + h^(2L)/2 - (y/l) h^(2L), u^{-1} delta B^y = 4 u^{-1} delta d^(1) - (8/l) delta d^(2) + 4 delta d^(2L) - (8y/l) delta d^(2L)",
         Flatten[{Table[lhsA[ch] - ((1/z) (coefY[coefZ[f[ch], 0], 1]/2 + yy coefY[coefZ[f[ch], 0], 2]) - coefY[coefZ[f[ch], 1], 0]/l
               + coefY[coefZ[f[ch], 1], 1]/2 - (yy/l) coefY[coefZ[f[ch], 1], 1]), {ch, chs}],
            lhsB - (4 (1/z) coefY[coefZ[f["dd"], 0], 1] - (8/l) coefY[coefZ[f["dd"], 1], 0] + 4 coefY[coefZ[f["dd"], 1], 1] - (8 yy/l) coefY[coefZ[f["dd"], 1], 1])}]]],
   {s, {"R", "NR"}}];
constraintRules[s_] := Flatten[Table[With[{Fp = Fptarget[s] /. slotRules[k] /. NRLpsi, Fm = Fmtarget[s] /. slotRules[k] /. NRLpsi,
      Rpk = Symbol["Rp" <> ToString[k]], Rmk = Symbol["Rm" <> ToString[k]]},
   {Derivative[i_, j_][Rpk][xp, xm] /; j >= 1 :> D[Fp, {xp, i}, {xm, j - 1}],
    Derivative[i_, j_][Rmk][xp, xm] /; i >= 1 :> D[Fm, {xp, i - 1}, {xm, j}]}], {k, 2}]];
pairing[X_, Y_] := -(X["mm"] Y["pp"] + X["pp"] Y["mm"] + X["mp"] Y["pm"] + X["pm"] Y["mp"]);   (* X^{a bbar} Y_{a bbar} *)

(* SMdirectcutoffbilinear: assembled from the general variation of the boundary term with the computed momenta,
   e^{-2d}[h1^{p qbar}(delta_2 A^y - 2 delta_2 d A^y) + 2 delta_1 d (delta_2 B^y - 2 delta_2 d (B^y + 4/l))]/(16 pi G) *)
orderedBilinear[s_] := Module[{f1 = solution[s, 1], f2 = solution[s, 2], m = If[s == "R", momR, momNR], dA2, dB2, A0, B0, sub2, e},
   sub2 = {hpp -> Function[{a, b, c}, Evaluate[f2["pp"] /. {xp -> a, xm -> b, yy -> c}]],
      hpm -> Function[{a, b, c}, Evaluate[f2["pm"] /. {xp -> a, xm -> b, yy -> c}]],
      hmp -> Function[{a, b, c}, Evaluate[f2["mp"] /. {xp -> a, xm -> b, yy -> c}]],
      hmm -> Function[{a, b, c}, Evaluate[f2["mm"] /. {xp -> a, xm -> b, yy -> c}]],
      dd -> Function[{a, b, c}, Evaluate[f2["dd"] /. {xp -> a, xm -> b, yy -> c}]]};
   (* the momentum variations are linear operators on the generic fluctuation; apply them to the second solution with
      the radial derivative d_y = d_yy - (2z/l) d_z acting on the z-dependent fields *)
   applyTo[op_] := Expand[op /. Derivative[i_, j_, k_][f_][xp, xm, yy] :> Nest[DyZY, D[f[xp, xm, yy] /. sub2, {xp, i}, {xm, j}], k]
      /. f_[xp, xm, yy] /; MemberQ[{hpp, hpm, hmp, hmm, dd}, f] :> (f[xp, xm, yy] /. sub2)];
   dA2 = Association[Table[ch -> applyTo[m["dA"][[Sequence @@ (ch /. {"pp" -> {1, 1}, "pm" -> {1, 2}, "mp" -> {2, 1}, "mm" -> {2, 2}})]]],
      {ch, {"pp", "pm", "mp", "mm"}}]];
   dB2 = applyTo[m["dB"]];
   A0 = Association[Table[ch -> m["A"][[Sequence @@ (ch /. {"pp" -> {1, 1}, "pm" -> {1, 2}, "mp" -> {2, 1}, "mm" -> {2, 2}})]], {ch, {"pp", "pm", "mp", "mm"}}]];
   B0 = m["B"];
   e = (1/z) (pairing[f1, dA2] - 2 f2["dd"] pairing[f1, A0] + 2 f1["dd"] (dB2 - 2 f2["dd"] (B0 + 4/l)));
   e = Expand[e /. NRLpsi];
   e = Sum[Coefficient[e, z, k] z^k, {k, -1, 0}];      (* terms vanishing as Y -> infinity dropped *)
   e];
displayedBilinear[s_] := Module[{f1 = solution[s, 1], f2 = solution[s, 2], dyf2, e},
   dyf2 = Association[Table[ch -> DyZY[f2[ch]], {ch, {"pp", "pm", "mp", "mm"}}]];
   e = (1/z) (pairing[f1, dyf2]/2 + 8 f1["dd"] DyZY[f2["dd"]])
      - (4/l) f2["dd"] ((Lp[xp] /. NRLpsi[[1]]) f1["mm"] + (Lm[xm] /. NRLpsi[[2]]) f1["pp"] + (Jsrc[s] /. NRLpsi) f1["mp"] + epsS[s] f1["pm"]);
   e = Expand[e /. NRLpsi];
   Sum[Coefficient[e, z, k] z^k, {k, -1, 0}]];
Print["  assembling the ordered cutoff bilinear on both saddles ..."];
Do[
   NRH`CheckZero["SMdirectcutoffbilinear on " <> s <> ": the general second variation with the computed momenta equals the displayed ordered bilinear through Y^0",
      Together[orderedBilinear[s] - displayedBilinear[s]]],
   {s, {"R", "NR"}}];

(* SMderivativectNR / SMderivativectR: quadratic counterterm on cutoff fields (Y = yy, e^{2Y/l} = 1/z) *)
ctQ[s_, A_, B_] := Module[{Dp, Dm, Vp, Vm, Q, LpS = Lp[xp] /. NRLpsi[[1]], LmS = Lm[xm] /. NRLpsi[[2]]},
   Dp[e_] := D[e, xp]; Dm[e_] := D[e, xm];
   Vp[e_] := LpS Dp[e] + 1/2 D[LpS, xp] e; Vm[e_] := LmS Dm[e] + 1/2 D[LmS, xm] e;
   Q = (1/z) ((l/2) A["mm"] Dp[Dp[B["mp"]]] + (l/2) A["pp"] Dm[Dm[B["mp"]]] - 2 l A["dd"] Dp[Dm[B["mp"]]] + (l^3/16) A["mp"] Dp[Dp[Dm[Dm[B["mp"]]]]])
      + 2 yy A["mm"] Vp[Dm[B["mp"]]] + 2 yy A["pp"] Vm[Dp[B["mp"]]]
      + (l yy^2/2) A["mp"] (Vp[Dm[Dm[Dm[B["mp"]]]]] + Vm[Dp[Dp[Dp[B["mp"]]]]]);
   If[s == "R", Q += l^2 yy A["dd"] Dp[Dp[Dm[Dm[B["mp"]]]]] + (l^3 yy^2/16) A["mp"] Dp[Dp[Dp[Dm[Dm[Dm[B["mp"]]]]]]]];
   Q];
ctBilinear[s_] := Module[{f1 = solution[s, 1], f2 = solution[s, 2], e},
   e = -1/2 (ctQ[s, f1, f2] + ctQ[s, f2, f1]);          (* delta_2 delta_1 S_ct^(2) = (1/16 pi G) (-1/2)(Q_12 + Q_21) *)
   e = Expand[e /. NRLpsi];
   Sum[Coefficient[e, z, k] z^k, {k, -1, 0}]];
divergentPart[e_] := Module[{ex = Expand[e]}, Coefficient[ex, z, -1]/z + Sum[Coefficient[Coefficient[ex, z, 0], yy, j] yy^j, {j, 1, Exponent[Coefficient[ex, z, 0], yy]}]];
finitePart[e_] := Coefficient[Coefficient[Expand[e], z, 0], yy, 0];
(* Integration-by-parts canonicalizer: every bilinear term is brought to the form j_{1,I}(x) q_I[j_2](x).  Slot-1 derivatives
   are moved onto the slot-2 factor; the nonlocal responses R_pm enter only through the stress constraints SMctlocality
   (d_- h^(2)_{op opb} = F_+, d_+ h^(2)_{om omb} = F_-), which are imposed in both slots.  The quotient by total tangential
   derivatives is thereby made explicit; the type-changing response H and the zero mode c_s are kept as leftovers. *)
slot1Fns = {a1, b1, r1, c1, v1, Rp1, Rm1, H1}; slot2Fns = {a2, b2, r2, c2, v2, Rp2, Rm2, H2};
termList[e_] := Which[e === 0, {}, Head[e] === Plus, List @@ e, True, {e}];
slotFactor[t_, fns_] := Module[{cs = Cases[t, e : (Derivative[_, _][f_][xp, xm] | f_[xp, xm]) /; MemberQ[fns, f] :> e, {0, Infinity}]},
   If[cs === {}, None, First[cs]]];
derivOrders[e_] := Replace[e, {Derivative[i_, j_][f_][xp, xm] :> {f, i, j}, f_[xp, xm] :> {f, 0, 0}}];
constraintRulesK[s_, k_] := With[{Fp = Fptarget[s] /. slotRules[k] /. NRLpsi, Fm = Fmtarget[s] /. slotRules[k] /. NRLpsi,
      Rpk = Symbol["Rp" <> ToString[k]], Rmk = Symbol["Rm" <> ToString[k]]},
   {Derivative[i_, j_][Rpk][xp, xm] /; j >= 1 :> D[Fp, {xp, i}, {xm, j - 1}],
    Derivative[i_, j_][Rmk][xp, xm] /; i >= 1 :> D[Fm, {xp, i - 1}, {xm, j}]}];
canonicalize[expr_, s_] := Module[{queue, q, left = 0, t, f1, X1, i, j, rest, f2, X2, i2, j2, f2p, c, Fp1, Fm1, need, var, rules2},
   Fp1 = Fptarget[s] /. slotRules[1] /. NRLpsi; Fm1 = Fmtarget[s] /. slotRules[1] /. NRLpsi;
   rules2 = constraintRulesK[s, 2];
   queue = termList[Expand[expr]];
   q = Association[# -> 0 & /@ {a1, b1, r1, c1, v1}];
   While[queue =!= {},
      t = Last[queue]; queue = Most[queue];
      f1 = slotFactor[t, slot1Fns];
      If[f1 === None,
         If[! FreeQ[t, cs1],
            f2 = slotFactor[t, slot2Fns];
            If[f2 === None, left += t; Continue[]];
            {X2, i2, j2} = derivOrders[f2];
            If[i2 + j2 == 0, left += t; Continue[]];
            If[i2 >= 1, f2p = D[X2[xp, xm], {xp, i2 - 1}, {xm, j2}]; var = xp, f2p = D[X2[xp, xm], {xp, i2}, {xm, j2 - 1}]; var = xm];
            c = t/f2; queue = Join[queue, termList[Expand[-D[c, var] f2p]]]; Continue[]];
         left += t; Continue[]];
      {X1, i, j} = derivOrders[f1];
      Which[
         X1 === Rp1 && j >= 1, queue = Join[queue, termList[Expand[(t/f1) D[Fp1, {xp, i}, {xm, j - 1}]]]],
         X1 === Rm1 && i >= 1, queue = Join[queue, termList[Expand[(t/f1) D[Fm1, {xp, i - 1}, {xm, j}]]]],
         MemberQ[{Rp1, Rm1}, X1],
            rest = (-1)^(i + j) D[t/f1, {xp, i}, {xm, j}];
            need = If[X1 === Rp1, xm, xp];
            Do[f2 = slotFactor[tt, slot2Fns];
               If[f2 === None, left += X1[xp, xm] tt; Continue[]];
               {X2, i2, j2} = derivOrders[f2];
               If[(need === xm && j2 == 0) || (need === xp && i2 == 0), left += X1[xp, xm] tt; Continue[]];
               f2p = If[need === xm, D[X2[xp, xm], {xp, i2}, {xm, j2 - 1}], D[X2[xp, xm], {xp, i2 - 1}, {xm, j2}]];
               c = tt/f2;
               queue = Join[queue, termList[Expand[-If[X1 === Rp1, Fp1, Fm1] c f2p - X1[xp, xm] D[c, need] f2p]]],
               {tt, termList[Expand[rest]]}],
         X1 === H1, left += t,
         True, q[X1] += (-1)^(i + j) D[t/f1, {xp, i}, {xm, j}]]];
   q = Map[Expand[# //. rules2] &, q]; left = Expand[left //. rules2];
   {q, left}];
allZeroQ[{q_, left_}] := And @@ (NRH`ZeroQ[Together[#]] & /@ Append[Values[q], left]);
Do[
   Module[{bil, ct, qc, qn},
      bil = orderedBilinear[s]; ct = ctBilinear[s];
      qc = canonicalize[divergentPart[bil] + divergentPart[ct], s];
      NRH`Check["SMderivativect on " <> s <> ": delta_2 delta_1 S_ct^(2) cancels every e^{2Y/l} poly(Y) and Y^k divergence of the ordered bilinear (mod total tangential derivatives, with the stress constraints imposed)",
         allZeroQ[qc]];
      qn = canonicalize[divergentPart[bil], s];
      NRH`Check["negative control on " <> s <> ": without the counterterm the divergences remain", ! allZeroQ[qn]];
      NRH`Check["SMctlocality on " <> s <> ": with the stress constraints the divergent rows of the bare bilinear are local in both slots (no undetermined responses" <> If[s == "NR", ", no W_1)", ")"],
         FreeQ[{Values[qn[[1]]], qn[[2]]}, Rp1 | Rm1 | H1 | Rp2 | Rm2 | H2] && (s == "R" || FreeQ[{Values[qn[[1]]], qn[[2]]}, W1])]],
   {s, {"R", "NR"}}];
NRH`CheckZero["SMctadjoint: (L+ d+ + (1/2) d+ L+)^dagger = -(L+ d+ + (1/2) d+ L+) on test functions",
   Module[{op, f = ff[xp], g = gg[xp]}, op[e_] := Lp[xp] D[e, xp] + 1/2 D[Lp[xp], xp] e;
      Together[g op[f] + f op[g] - D[Lp[xp] f g, xp]]]];

(* SMfullfinitehessian, SMdirectstressrows, SMdirectremainingrows: historical finite-source-scheme Hessian and contact rows, retained as ancillary checks. *)
finiteDisplay[s_] := Module[{f1 = solution[s, 1], f2 = solution[s, 2], p1, p2, coefY, coefZ},
   coefZ[e_, n_] := Coefficient[Expand[e], z, n]; coefY[e_, j_] := Coefficient[Expand[e], yy, j];
   p1 = Association[Table[ch -> <|"0" -> coefY[coefZ[f1[ch], 0], 0], "1" -> coefY[coefZ[f1[ch], 0], 1], "2" -> coefY[coefZ[f1[ch], 1], 0], "2L" -> coefY[coefZ[f1[ch], 1], 1]|>,
      {ch, {"pp", "pm", "mp", "mm", "dd"}}]];
   p2 = Association[Table[ch -> <|"0" -> coefY[coefZ[f2[ch], 0], 0], "1" -> coefY[coefZ[f2[ch], 0], 1], "2" -> coefY[coefZ[f2[ch], 1], 0], "2L" -> coefY[coefZ[f2[ch], 1], 1]|>,
      {ch, {"pp", "pm", "mp", "mm", "dd"}}]];
   Expand[pairing[Association[Table[ch -> p1[ch]["0"], {ch, {"pp", "pm", "mp", "mm"}}]], Association[Table[ch -> -(1/l) p2[ch]["2"] + p2[ch]["2L"]/2, {ch, {"pp", "pm", "mp", "mm"}}]]]
      + 1/2 pairing[Association[Table[ch -> p1[ch]["2"], {ch, {"pp", "pm", "mp", "mm"}}]], Association[Table[ch -> p2[ch]["1"], {ch, {"pp", "pm", "mp", "mm"}}]]]
      + p1["dd"]["0"] (-(16/l) p2["dd"]["2"] + 8 p2["dd"]["2L"]) + 8 p1["dd"]["2"] p2["dd"]["1"]
      - (4/l) p2["dd"]["0"] ((Lp[xp] /. NRLpsi[[1]]) a1[xp, xm] + (Lm[xm] /. NRLpsi[[2]]) b1[xp, xm] + (Jsrc[s] /. NRLpsi) r1[xp, xm] + epsS[s] c1[xp, xm]) /. NRLpsi]];
(* SMcontactconstraintreduction: integrate the slot-1 response terms by parts and use the constraints *)
localize[e_, s_] := Module[{ex = Expand[e], rules},
   rules = {Rp1[xp, xm] Derivative[0, n_][r2][xp, xm] /; n >= 1 :> (-1)^n r2[xp, xm] D[Rp1[xp, xm], {xm, n}],
      Rm1[xp, xm] Derivative[n_, 0][r2][xp, xm] /; n >= 1 :> (-1)^n r2[xp, xm] D[Rm1[xp, xm], {xp, n}]};
   Expand[(ex /. rules) //. constraintRules[s]]];
(* Historical rows q_{s,I} with slot-2 arguments (SMdirectstressrows, SMdirectremainingrows); no longer displayed in the compact SM. *)
qDisplayed[s_] := Module[{eps = epsS[s], J = Jsrc[s], D2 = Function[{f, v}, D[f, Sequence @@ v]], a = a2[xp, xm], b = b2[xp, xm], r = r2[xp, xm], c = c2[xp, xm], v = v2[xp, xm], LpS, LmS},
   LpS = Lp[xp]; LmS = Lm[xm];
   <|"a" -> Rp2[xp, xm]/l + eps l^3/16 D[r, {xp, 3}, xm] - 3 l/4 LpS D[r, xp, xm] - l/2 D[LpS, xp] D[r, xm] - 4/l LpS v,
     "b" -> Rm2[xp, xm]/l + eps l^3/16 D[r, xp, {xm, 3}] - 3 l/4 LmS D[r, xp, xm] - l/2 D[LmS, xm] D[r, xp] - 4/l LmS v,
     "r" -> H2[xp, xm]/l - eps l^3/16 (D[a, {xp, 3}, xm] + D[b, xp, {xm, 3}]) + eps l^3/4 D[v, {xp, 2}, {xm, 2}] - eps l^5/64 D[r, {xp, 3}, {xm, 3}]
        - l^3/16 (D[LmS, xm] D[r, {xp, 3}] + D[LpS, xp] D[r, {xm, 3}]) + l/4 J D[r, xp, xm]
        + l/2 (D[J, xp] D[r, xm] + D[J, xm] D[r, xp] + D[J, xp, xm] r) - 4/l J v,
     "c" -> cs2/l,
     "v" -> eps (-l (D[a, {xp, 2}] + D[b, {xm, 2}] - 4 D[v, xp, xm]) + l^3/4 D[r, {xp, 2}, {xm, 2}]) - l (LpS D[r, {xm, 2}] + LmS D[r, {xp, 2}])|> /. NRLpsi];
Do[
   Module[{bil, fin, loc, rows, target, ok},
      bil = orderedBilinear[s];
      fin = finitePart[bil];
      NRH`CheckZero["SMfullfinitehessian on " <> s <> ": the Y^0 coefficient of the ordered bilinear equals the displayed finite integrand",
         Together[fin - finiteDisplay[s]]];
      loc = localize[fin, s];
      NRH`Check["SMcontactconstraintreduction on " <> s <> ": after integrating the slot-1 response terms by parts the finite part is local in the first source",
         FreeQ[loc, Rp1] && FreeQ[loc, Rm1] && FreeQ[loc, H1]];
      rows = Association[Table[n -> Together[VariationalD[loc, Symbol[n <> "1"][xp, xm], {xp, xm}] //. constraintRules[s]], {n, {"a", "b", "r", "c", "v"}}]];
      target = qDisplayed[s];
      ok = AllTrue[{"a", "b", "r", "c", "v"}, Together[rows[#] - target[#]] === 0 &];
      If[! ok, Do[Print["   row ", n, " on ", s, ": ", Together[rows[n] - target[n]]], {n, {"a", "b", "r", "c", "v"}}]];
      NRH`Check["SMdirectstressrows, SMdirectremainingrows on " <> s <> ": the Euler-Lagrange derivatives of the finite bilinear with respect to the five first-slot sources are the displayed rows q_{s,I}", ok];
      NRH`CheckZero["SMdirectkernel, SM" <> s <> "scalarcontactrow: (G_s)_{5J}|action response = q_5/(4 * 16 pi G) reproduces the displayed scalar contact row",
         Together[rows["v"]/(4 16 Pi G) - If[s == "R",
            l/(64 Pi G) (-D[a2[xp, xm], {xp, 2}] - D[b2[xp, xm], {xm, 2}] + (l^2/4) D[r2[xp, xm], {xp, 2}, {xm, 2}] - Lp[xp] D[r2[xp, xm], {xm, 2}] - Lm[xm] D[r2[xp, xm], {xp, 2}] + 4 D[v2[xp, xm], xp, xm]),
            -l/(64 Pi G) (Lp[xp] D[r2[xp, xm], {xm, 2}] + Lm[xm] D[r2[xp, xm], {xp, 2}])] /. NRLpsi]];
      NRH`Check["fourth row on " <> s <> ": the fourth h_pm/W0-source row couples only to c_s/l (the 4 delta d^(0)/l terms cancel against the measure term)",
         Together[rows["c"] - cs2/l] === 0]],
   {s, {"R", "NR"}}];

(* SMLorentziandeltaconvention: normalization of the joint Lorentzian distribution *)
NRH`CheckZero["SMLorentziandeltaconvention: the box integral of (1/(2 pi i)) d+ d- ln(-x+ x- + i eps) over [-R,R]^2 is 1",
   Limit[1/(2 Pi I) (Log[-R^2 + I ep] - Log[R^2 + I ep] - Log[R^2 + I ep] + Log[-R^2 + I ep]), ep -> 0, Direction -> "FromAbove", Assumptions -> R > 0] - 1];
NRH`CheckZero["SMPBHkernel: d_-^{-1} delta^2 = (1/(2 pi i)) d_+ ln(-D+ D- + i eps) -> 1/(2 pi i D+) away from the light cone",
   Together[1/(2 Pi I) D[Log[-dp dm], dp] - 1/(2 Pi I dp)]];

(* SMstresssourceexample, SMgeneralWardoperators, SMstressparticular: kernels from the rows *)
Module[{Fp, Vnr, Vr, kernelA, aNR, aR, mNR, target},
   Fp = Fptarget["R"] /. {b0 -> Function[{x, y}, 0], r0 -> Function[{x, y}, 0], c0 -> Function[{x, y}, 0], v0 -> Function[{x, y}, 0]} /. a0 -> aa;     (* only h^(0)_{om omb} = a turned on *)
   Vr[e_] := 2 Lp[xp] D[e, xp] + D[Lp[xp], xp] e - l^2/4 D[e, {xp, 3}];
   Vnr[e_] := 2 Lp[xp] D[e, xp] + D[Lp[xp], xp] e;
   NRH`CheckZero["SMstresssourceexample: with only h^(0)_{om omb} the R stress constraint is d_- h^(2)_{op opb} = -(2L+ d+ + d+L+ - (l^2/4) d+^3) h^(0)_{om omb}",
      Together[Fp + Vr[aa[xp, xm]]]];
   NRH`CheckZero["SMstresssourceexample: on NR the third-derivative term is absent",
      Together[(Fptarget["NR"] /. {b0 -> Function[{x, y}, 0], r0 -> Function[{x, y}, 0], c0 -> Function[{x, y}, 0], v0 -> Function[{x, y}, 0]} /. a0 -> aa) + Vnr[aa[xp, xm]]]];
   (* q_{s,1} = R+/l + local; R+ = -V+ d_-^{-1} a + f+  =>  Q_11 = -(1/l) V+ d_-^{-1} delta, G = Q/(64 pi G), time ordered: /i *)
   NRH`CheckZero["SMstressparticular: -(1/l)/(64 pi G)/i = -kappa/(4 i) with kappa = 1/(16 pi G l)",
      Together[-(1/l)/(64 Pi G)/I + 1/(16 Pi G l)/(4 I)]];
   (* position space with the Lorentzian d_-^{-1} delta -> 1/(2 pi i D+) *)
   (* on kernels of the separation D+ = x+ - x'+, d_+ acts as d/dD+ while L+ is evaluated at the first point *)
   VrK[f_] := 2 Lp[xp] D[f, dp] + D[Lp[xp], xp] f - l^2/4 D[f, {dp, 3}];
   VnrK[f_] := 2 Lp[xp] D[f, dp] + D[Lp[xp], xp] f;
   kernelA[V_] := -1/(16 Pi G l)/(4 I) V[1/(2 Pi I dp)];
   aNR = Together[kernelA[VnrK]]; aR = Together[kernelA[VrK]];
   target = (D[Lp[xp], xp]/dp - 2 Lp[xp]/dp^2)/(128 Pi^2 G l);
   NRH`CheckZero["SMgeneralpositionkernels, Rcorrelators: A_+^NR = [L+'/D+ - 2L+/D+^2]/(128 pi^2 G l)", aNR - target];
   NRH`CheckZero["SMgeneralpositionkernels, Rcorrelators: A_+^R = A_+^NR + 3l/(256 pi^2 G D+^4)", aR - target - 3 l/(256 Pi^2 G dp^4)];
   NRH`CheckZero["SMRtwopt: at L+ = 0, (8 pi)^2 A_+^R = (c/2)/D+^4 with c = 3l/(2G)",
      Together[(8 Pi)^2 (aR /. Lp -> (0 &)) - (3 l/(2 G))/2/dp^4]];
   NRH`Check["the Lorentzian route is real: the 1/i of the action variation and the 1/i of d_-^{-1} delta combine to -1/(2 pi)",
      Together[1/I 1/(2 Pi I) + 1/(2 Pi)] === 0]];

(* SMindependentsourcegenerator, SMhairparticularsolution, SMgeneralparticularkernels: the nonchiral completion *)
Module[{alpha = al[xp, xm], xiy, xiR, xiNR, rEx = RiemannianSaddleExact[], nrEx, sR, sNR, dHR, dHNR, hR, hNR, ddR, ddNR, Vp, Wp, Hr, Hnr},
   xiy = -(l/2) D[alpha, xp];
   (* xi = (xi~_+, xi~_-, xi~_y; xi^+, xi^-, xi^y): only alpha^+ turned on, radial-gauge completion through u *)
   xiR = {-yy D[xiy, xp] + l z Lp[xp] D[xiy, xm], yy D[xiy, xm] - l z Lm[xm] D[xiy, xp], 0, alpha - (l/2) z D[xiy, xm], -(l/2) z D[xiy, xp], xiy};   (* xi^-_R = alpha^- - (l/2) u d_+ xi^y with alpha^- = 0 *)
   xiNR = {-yy D[xiy, xp] + l z (1/psip[xp]^2) D[xiy, xm], yy D[xiy, xm] - l z (1/psim[xm]^2) D[xiy, xp], 0, alpha, 0, xiy};
   sR = SaddleSeries[rEx, 1]; nrEx = NonRiemannianSaddleExact[W1[xp, xm] z + NRW2 z^2]; sNR = SaddleSeries[nrEx, 1];
   dHR = SeriesZM[GenLieH[xiR, sR["H"], xsZY], 1]; dHNR = SeriesZM[GenLieH[xiNR, sNR["H"], xsZY], 1];
   hR = SeriesZM[Transpose[JJ . sR["V"]] . dHR . (JJ . sR["Vb"]), 1];
   hNR = SeriesZM[Transpose[JJ . sNR["V"]] . dHNR . (JJ . sNR["Vb"]), 1];
   ddR = SeriesZ[GenLieD[xiR, sR["d"], xsZY], 1]; ddNR = SeriesZ[GenLieD[xiNR, sNR["d"], xsZY], 1];
   NRH`CheckZero["SMindependentsourcegenerator: radial gauge h_{p ybar} = h_{y qbar} = 0 through u on both saddles, and delta d^(0) = 0",
      {hR[[All, 3]], hR[[3, All]], hNR[[All, 3]], hNR[[3, All]], Coefficient[ddR, z, 0], Coefficient[ddNR, z, 0]}];
   NRH`CheckZero["SMindependentsourcegenerator: the only leading source is h^(0)_{om omb} = -2 d_- alpha^+ (both saddles)",
      {(Coefficient[hR[[1 ;; 2, 1 ;; 2]], z, 0] /. yy -> 0) - {{0, 0}, {0, -2 D[alpha, xm]}}, (Coefficient[hNR[[1 ;; 2, 1 ;; 2]], z, 0] /. yy -> 0) - {{0, 0}, {0, -2 D[alpha, xm]}}}];
   NRH`CheckZero["SMindependentsourcegenerator: the y-linear piece of the leading block is the current SM3.4-3.5 response h^(1)_{op omb} = -(l/2) d_+^2 h^(0)_{om omb} = l d_+^2 d_- alpha^+ (no y^2 term)",
      {Coefficient[Coefficient[hR[[1 ;; 2, 1 ;; 2]], z, 0], yy, 1] - {{0, l D[alpha, {xp, 2}, xm]}, {0, 0}}, Coefficient[Coefficient[hNR[[1 ;; 2, 1 ;; 2]], z, 0], yy, 1] - {{0, l D[alpha, {xp, 2}, xm]}, {0, 0}},
       Coefficient[Coefficient[hR[[1 ;; 2, 1 ;; 2]], z, 0], yy, 2], Coefficient[Coefficient[hNR[[1 ;; 2, 1 ;; 2]], z, 0], yy, 2]}];
   Vp[e_, s_] := 2 Lp[xp] D[e, xp] + D[Lp[xp], xp] e - If[s == "R", l^2/4 D[e, {xp, 3}], 0];
   Wp[e_] := 2 W1[xp, xm] D[e, xp] + D[W1[xp, xm], xp] e - l^2 Lm[xm] D[e, {xp, 3}];
   NRH`CheckZero["SMhairparticularsolution: h^(2)_{op opb} = 2 V_+^s alpha^+ on both saddles",
      {Coefficient[hR[[1, 1]], z, 1] - 2 Vp[alpha, "R"], Coefficient[hNR[[1, 1]], z, 1] - 2 (Vp[alpha, "NR"] /. NRLpsi)}];
   NRH`CheckZero["SMhairparticularsolution: h^(2)_{op omb}|R = 2 L- V_+^R alpha - (l^2/2) L+ d+ d-^2 alpha, |NR = (1/2) W_+ alpha - (l^2/2) L+ d+ d-^2 alpha",
      {Coefficient[hR[[1, 2]], z, 1] - (2 Lm[xm] Vp[alpha, "R"] - l^2/2 Lp[xp] D[alpha, xp, {xm, 2}]),
       Coefficient[hNR[[1, 2]], z, 1] - ((Wp[alpha]/2 - l^2/2 Lp[xp] D[alpha, xp, {xm, 2}]) /. NRLpsi)}];
   (* mixed kernels: q_{s,3} contains H_s/l; alpha = -(1/2) d_-^{-1} h^(0)_{om omb}; time ordering /i *)
   (* with alpha^+ = -(1/2) d_-^{-1} h^(0)_{om omb} represented by a generic function gF *)
   Hr = 2 Lm[xm] Vp[-1/2 gF[xp, xm], "R"]; Hnr = Wp[-1/2 gF[xp, xm]]/2;
   NRH`CheckZero["SMgeneralparticularkernels: M_+^R = L_- A_+^R and M_+^NR = -(kappa/(16 i)) W_+ d_-^{-1} delta (normalization chain)",
      {Together[(Hr/l)/(64 Pi G)/I - Lm[xm] (-(1/(16 Pi G l))/(4 I) Vp[gF[xp, xm], "R"])], Together[(Hnr/l)/(64 Pi G)/I + 1/(16 Pi G l)/(16 I) Wp[gF[xp, xm]]]}];
   WpK[f_] := 2 W1[xp, xm] D[f, dp] + D[W1[xp, xm], xp] f - l^2 Lm[xm] D[f, {dp, 3}];
   NRH`CheckZero["SMgeneralpositionkernels, Rcorrelators: M_+^NR = [W1'/D+ - 2W1/D+^2 + 6 l^2 L-/D+^4]/(512 pi^2 G l) with d_-^{-1} delta -> 1/(2 pi i D+)",
      Together[-(1/(16 Pi G l))/(16 I) WpK[1/(2 Pi I dp)] - (D[W1[xp, xm], xp]/dp - 2 W1[xp, xm]/dp^2 + 6 l^2 Lm[xm]/dp^4)/(512 Pi^2 G l)]];
   NRH`CheckZero["SMNRhairhessian: the hair-hair coefficient (1/l)/(64 pi G)/i = 1/(64 pi i G l)",
      Together[(1/l)/(64 Pi G)/I - 1/(64 Pi I G l)]]];

(* SMRhairuniformization, SMRhairconditional, SMRhairconstantL: historical conditional
   uniformization identities, retained as ancillary algebra. These factors do not establish
   the current SM3.3 K3 interior prescription or derive its H_R response and hair coefficient. *)
Module[{f = Exp[2 alp xx], schw, Bpm, BpmC},
   schw = D[f, {xx, 3}]/D[f, xx] - 3/2 (D[f, {xx, 2}]/D[f, xx])^2;
   NRH`CheckZero["SMRhairuniformization: {e^{2 alpha x}; x} = -2 alpha^2 = -4 L/l^2 for alpha = Sqrt[2L]/l",
      Simplify[schw + 4 (alp^2 l^2/2)/l^2]];
   Bpm = D[f, xx]^2 (D[f, xx] /. xx -> xx2)^2/(f - (f /. xx -> xx2))^4;
   NRH`CheckZero["SMRhairconstantL: B(x,x') = [alpha/sinh(alpha D)]^4 for f = e^{2 alpha x}",
      Simplify[Bpm - (alp/Sinh[alp (xx - xx2)])^4]];
   NRH`CheckZero["SMRhairconstantL: the L -> 0 limit of the thermal factor is 1/D^4",
      Limit[(alp/Sinh[alp dp])^4, alp -> 0] - 1/dp^4];
   NRH`CheckZero["SMRhairconditional: the stress coefficient 3l/(256 pi^2 G) equals (8 pi)^{-2} c/2",
      Together[3 l/(256 Pi^2 G) - 1/(8 Pi)^2 (3 l/(2 G))/2]];
   NRH`CheckZero["SMRhairconditional: algebraic rewriting only: 9 l^5/(64 pi^2 G r^8) with r^2 = -2 D+ D- is 9 l^5/(1024 pi^2 G) (-D+ D-)^{-4}",
      Together[9 l^5/(64 Pi^2 G (-2 dp dm)^4) - 9 l^5/(1024 Pi^2 G) (-dp dm)^-4]]];
NRH`CheckZero["SMonepointequivalence: (1/2) delta O^_I/delta j_J with O^_I = (1/2) delta S/delta j_I is (1/4) delta^2 S/delta j_I delta j_J",
   Together[1/2 D[1/2 D[SS[j1, j2], j1], j2] - 1/4 D[SS[j1, j2], j1, j2]]];

NRH`FileSummary[];
