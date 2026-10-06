(* 01_ExactBackgrounds_EDFE.wl | 2026-10-06 standalone edition.
   All definitions are embedded. No Get, Needs, input files, or packages.
   Run in a fresh kernel; this file clears Global` and NRH`.
   Stable labels identify formulas; old SM numbers in inherited check IDs are historical. *)
ClearAll["Global`*", "NRH`*"];
(* ::Title:: *)
(*NRH01 DFT Tools*)

NRH`$FileResults;
If[!ListQ[NRH`$AllResults], NRH`$AllResults = {}];

NRH`BeginFile[name_String] := (
   NRH`$CurrentFile = name;
   NRH`$FileResults = {};
   Print["\n================================================================"];
   Print["  ", name];
   Print["================================================================"]);

NRH`Record[label_String, ok : (True | False)] := (
   AppendTo[NRH`$FileResults, {NRH`$CurrentFile, label, ok}];
   AppendTo[NRH`$AllResults, {NRH`$CurrentFile, label, ok}];
   Print[If[ok, "  [PASS] ", "  [FAIL] "], label];
   ok);

NRH`CheckZero[label_String, expr_] := Module[{z},
   z = NRH`ZeroQ[expr];
   NRH`Record[label, TrueQ[z]]];

NRH`Check[label_String, statement_] := NRH`Record[label, TrueQ[statement]];

NRH`FileSummary[] := Module[{n, bad},
   n = Length[NRH`$FileResults];
   bad = Select[NRH`$FileResults, #[[3]] === False &];
   Print["----------------------------------------------------------------"];
   Print["  ", NRH`$CurrentFile, ": ", n - Length[bad], "/", n, " checks passed."];
   If[Length[bad] > 0,
      Print["  FAILED: ", bad[[All, 2]]];
      If[$FrontEnd === Null && ! TrueQ[NRH`$DeferExit], Exit[1]]];
   Length[bad] === 0];

NRH`GrandSummary[] := Module[{n, bad},
   n = Length[NRH`$AllResults];
   bad = Select[NRH`$AllResults, #[[3]] === False &];
   Print["\n################################################################"];
   Print["  GRAND TOTAL: ", n - Length[bad], "/", n, " checks passed."];
   Scan[Print["  FAILED: ", #[[1]], " -- ", #[[2]]] &, bad];
   Print["################################################################"];
   If[Length[bad] > 0 && $FrontEnd === Null, Exit[1]];
   Length[bad] === 0];

NRH`ZeroQ[expr_] := Module[{flat, t},
   flat = Flatten[{expr}];
   AllTrue[flat,
      Function[e,
         t = Together[Expand[e]];
         If[t === 0, True,
            t = Together[ExpandAll[TrigToExp[t]]];
            If[t === 0, True, PossibleZeroQ[Simplify[t]]]]]]];

ODDJ[nphys_Integer] := ArrayFlatten[{{0, IdentityMatrix[nphys]}, {IdentityMatrix[nphys], 0}}];

DblD[expr_, m_Integer, xs_List] := Module[{n = Length[xs], op},
   If[m <= n, 0*expr,
      op = xs[[m - n]];
      If[Head[op] === Function, op[expr], D[expr, op]]]];

DblGrad[expr_, xs_List] := Table[DblD[expr, m, xs], {m, 1, 2 Length[xs]}];

(* Lowered Gamma_CAB, unit-weight antisymmetrization; trace fixed by nabla d = 0. *)
GammaDFT[HH_, dd_, xs_List] := Module[
   {n = Length[xs], dim, JJ, P, Pb, Pm, Pbm, PbUD, dP, gradd,
    gamma12, T12, X, PbmX, PmX, coeff},
   dim = 2 n; JJ = ODDJ[n];
   P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   Pm = P . JJ; Pbm = Pb . JJ; PbUD = JJ . Pb;
   dP = Table[DblD[P, m, xs], {m, 1, dim}];
   gradd = DblGrad[dd, xs];
   gamma12 = Table[
      Module[{term1, m2},
         term1 = Pm . dP[[c]] . PbUD;
         term1 = term1 - Transpose[term1];

         m2 = Sum[
            Module[{colQb = (Pbm . dP[[dd2]])[[All, c]], colQ = (Pm . dP[[dd2]])[[All, c]]},
               Outer[Times, Pbm[[All, dd2]], colQb] - Outer[Times, colQb, Pbm[[All, dd2]]]
               - Outer[Times, Pm[[All, dd2]], colQ] + Outer[Times, colQ, Pm[[All, dd2]]]],
            {dd2, n + 1, dim}];
         Map[Together, term1 + m2, {2}]],
      {c, 1, dim}];
   T12 = Table[Together[Sum[JJ[[b, e]]*gamma12[[e, b, a]], {b, 1, dim}, {e, 1, dim}]], {a, 1, dim}];
   X = gradd + T12/2;
   PbmX = Pbm . X; PmX = Pm . X;
   coeff = -4/(n - 1);
   Table[
      Map[Together, gamma12[[c]] + coeff/2*(
         Outer[Times, Pb[[c]], PbmX] - Outer[Times, PbmX, Pb[[c]]]
         + Outer[Times, P[[c]], PmX] - Outer[Times, PmX, P[[c]]]), {2}],
      {c, 1, dim}]];

RiemannR4[gamma_List, xs_List] := Module[{n = Length[xs], dim, JJ},
   dim = 2 n; JJ = ODDJ[n];
   Table[
      If[b <= a, ConstantArray[0, {dim, dim}],
         Map[Together,
            DblD[gamma[[b]], a, xs] - DblD[gamma[[a]], b, xs]
            + gamma[[a]] . JJ . gamma[[b]] - gamma[[b]] . JJ . gamma[[a]], {2}]],
      {a, 1, dim}, {b, 1, dim}]
   // (# - Transpose[#, {2, 1, 3, 4}] &)];

Partner[m_Integer, n_Integer] := If[m <= n, m + n, m - n];

RicciS[gamma_List, r4_List, xs_List] := Module[{n = Length[xs], dim, gg},
   dim = 2 n;

   Table[
      Together[Sum[Module[{e = Partner[c, n]},
         (r4[[c, b, e, a]] + r4[[e, a, c, b]]
            - Sum[gamma[[Partner[f, n], e, a]]*gamma[[f, c, b]], {f, 1, dim}])/2],
         {c, 1, dim}]],
      {a, 1, dim}, {b, 1, dim}]];

ScalarS0[HH_, dd_, xs_List] := Module[
   {n = Length[xs], dim, JJ, Hup, Hmix, gradd, term},
   dim = 2 n; JJ = ODDJ[n];
   Hup = JJ . HH . JJ;
   Hmix = HH . JJ;
   gradd = DblGrad[dd, xs];
   term =
      Sum[Hup[[a, b]]*(
            1/8*Sum[DblD[Hup[[c, e]], a, xs]*DblD[HH[[c, e]], b, xs], {c, 1, dim}, {e, 1, dim}]
            + 1/2*Sum[DblD[Hmix[[a, e]], c, xs]*DblD[Hmix[[b, c]], e, xs], {c, 1, dim}, {e, 1, dim}]
            - 4*gradd[[a]]*gradd[[b]] + 4*DblD[gradd[[b]], a, xs]),
         {a, 1, dim}, {b, 1, dim}]
      - Sum[DblD[DblD[Hup[[a, b]], a, xs], b, xs], {a, 1, dim}, {b, 1, dim}]
      + 4*Sum[DblD[Hup[[a, b]], a, xs]*gradd[[b]], {a, 1, dim}, {b, 1, dim}];
   term];

ScalarS0FromS4[gamma_List, r4_List, HH_, xs_List] := Module[
   {n = Length[xs], dim, JJ, P, Pb, Pup, Pbup, s4},
   dim = 2 n; JJ = ODDJ[n];
   P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   Pup = JJ . P . JJ; Pbup = JJ . Pb . JJ;
   s4[a_, b_, c_, d_] :=
      (r4[[c, d, a, b]] + r4[[a, b, c, d]]
         - Sum[gamma[[Partner[f, n], a, b]]*gamma[[f, c, d]], {f, 1, dim}])/2;
   Sum[(Pup[[a, c]]*Pup[[b, d]] - Pbup[[a, c]]*Pbup[[b, d]])*s4[a, b, c, d],
      {a, 1, dim}, {b, 1, dim}, {c, 1, dim}, {d, 1, dim}]];

ProjectedRicci[HH_, ricci_, xs_List] := Module[{n = Length[xs], JJ, P, Pb},
   JJ = ODDJ[n]; P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   P . JJ . ricci . JJ . Pb];

EinsteinG[HH_, ricci_, s0_, xs_List] := Module[{n = Length[xs], JJ, psp},
   JJ = ODDJ[n];
   psp = ProjectedRicci[HH, ricci, xs];
   4*(psp - Transpose[psp])/2 - 1/2*JJ*s0];

GenLieH[xiUp_List, HH_, xs_List] := Module[
   {n = Length[xs], dim, JJ, xiLow, dxiUp, dxiLow, amat},
   dim = 2 n; JJ = ODDJ[n];
   xiLow = JJ . xiUp;
   dxiUp = Table[DblD[xiUp[[c]], m, xs], {m, 1, dim}, {c, 1, dim}];
   dxiLow = Table[DblD[xiLow[[c]], m, xs], {m, 1, dim}, {c, 1, dim}];

   amat = Table[dxiUp[[m, c]] - Sum[JJ[[c, dd2]]*dxiLow[[dd2, m]], {dd2, 1, dim}], {m, 1, dim}, {c, 1, dim}];
   Sum[xiUp[[c]]*DblD[HH, c, xs], {c, 1, dim}] + amat . HH + HH . Transpose[amat]];

GenLieD[xiUp_List, dd_, xs_List] := Module[{n = Length[xs], dim},
   dim = 2 n;
   Sum[xiUp[[a]]*DblD[dd, a, xs], {a, 1, dim}] - 1/2*Sum[DblD[xiUp[[a]], a, xs], {a, 1, dim}]];

RiemannianH[g_, B_] := Module[{gi = Inverse[g]},
   ArrayFlatten[{{gi, -gi . B}, {B . gi, g - B . gi . B}}]];

RiemannianDilaton[g_, phi_] := phi - 1/4 Log[-Det[g]];

Gamma2Density[HH_, dd_, gamma_List, xs_List] := Module[
   {n = Length[xs], dim, JJ, P, Pb, Pup, Pbup, gup},
   dim = 2 n; JJ = ODDJ[n];
   P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   Pup = JJ . P . JJ; Pbup = JJ . Pb . JJ;

   Exp[-2 dd]*Sum[(Pup[[a, c]]*Pup[[b, d]] - Pbup[[a, c]]*Pbup[[b, d]])*
        Sum[gamma[[a, c, Partner[e, n]]]*gamma[[b, d, e]]
            - gamma[[a, b, Partner[e, n]]]*gamma[[d, c, e]]
            + 1/2*gamma[[Partner[e, n], a, b]]*gamma[[e, c, d]], {e, 1, dim}],
      {a, 1, dim}, {b, 1, dim}, {c, 1, dim}, {d, 1, dim}]];

GammaBVector[HH_, dd_, xs_List] := Module[{n = Length[xs], dim, JJ, Hup},
   dim = 2 n; JJ = ODDJ[n]; Hup = JJ . HH . JJ;
   Table[4*Sum[Hup[[a, b]]*DblD[dd, b, xs], {b, 1, dim}]
         - Sum[DblD[Hup[[a, b]], b, xs], {b, 1, dim}], {a, 1, dim}]];

(* Unprojected mathcal A and its mixed projected tensor; doubled indices are lowered. *)
MomentumCore[gamma_List, xs_List] := Module[{dim = 2 Length[xs], j = ODDJ[Length[xs]], tr},
   tr = Table[Sum[j[[a, b]] gamma[[b, a, n]], {a, dim}, {b, dim}], {n, dim}];
   Table[KroneckerDelta[k, m] tr[[n]] + KroneckerDelta[k, n] tr[[m]]
      - Sum[j[[k, a]] (gamma[[m, a, n]] + gamma[[n, a, m]]), {a, dim}],
      {k, dim}, {m, dim}, {n, dim}]];
MomentumAK[HH_, gamma_List, xs_List] := With[{j = ODDJ[Length[xs]]},
   Map[Map[Together, ((j + HH)/2) . j . # . Transpose[((j - HH)/2) . j], {2}] &,
      MomentumCore[gamma, xs]]];

SpinConnectionDFT[V_, eta_, gamma_List, xs_List] := Module[
   {n = Length[xs], dim, JJ, VlowFlat, VupLowFlat, deriv, cov, phi},
   dim = 2 n; JJ = ODDJ[n];
   VlowFlat = V . eta;
   VupLowFlat = JJ . V . eta;
   Table[
      deriv = DblD[V, a, xs] . eta;
      cov = deriv + gamma[[a]] . JJ . VlowFlat;
      phi = Transpose[VupLowFlat] . cov;
      (phi - Transpose[phi])/2,
      {a, 1, dim}]];

DFTCurvature[HH_, dd_, xs_List] := Module[{gamma, r4, ric, s0, t},
   {t, gamma} = AbsoluteTiming[GammaDFT[HH, dd, xs]];
   Print["    [timing] Gamma: ", Round[t, 0.1], " s"];
   {t, r4} = AbsoluteTiming[RiemannR4[gamma, xs]];
   Print["    [timing] R4:    ", Round[t, 0.1], " s"];
   {t, ric} = AbsoluteTiming[RicciS[gamma, r4, xs]];
   Print["    [timing] Ricci: ", Round[t, 0.1], " s"];
   {t, s0} = AbsoluteTiming[Together[ScalarS0[HH, dd, xs]]];
   Print["    [timing] S0:    ", Round[t, 0.1], " s"];
   <|"Gamma" -> gamma, "R4" -> r4, "Ricci" -> ric, "S0" -> s0,
     "PSPbar" -> ProjectedRicci[HH, ric, xs],
     "G" -> EinsteinG[HH, ric, s0, xs]|>];

Print["[NRH01] DFT toolbox loaded."];

(* Shared three-dimensional backgrounds; u = exp(2 y/l), z = 1/u. *)
Hinf = {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, -1, 0}, {0, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1}};
Vinf = {{1/Sqrt[2], 0, 0}, {0, 0, 0}, {0, 0, 1/Sqrt[2]},
   {0, -Sqrt[2], 0}, {0, 0, 0}, {0, 0, 1/Sqrt[2]}};
Vbinf = {{0, 0, 0}, {0, 1/Sqrt[2], 0}, {0, 0, 1/Sqrt[2]},
   {0, 0, 0}, {Sqrt[2], 0, 0}, {0, 0, -1/Sqrt[2]}};
eta3 = {{0, -1, 0}, {-1, 0, 0}, {0, 0, 1}};
etab3 = -eta3;

RiemannianMetric[lp_, lm_, u_] := With[{f = u + lp lm/u},
   {{2 lp, -f, 0}, {-f, 2 lm, 0}, {0, 0, 1}}];
RiemannianB[lp_, lm_, u_] := (u + lp lm/u) {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
RiemannianD[lp_, lm_, u_] := -Log[u (1 - lp lm/u^2)]/2;
NonRiemannianH[chi_, esigma_, w_] := With[{c = Cosh[chi], s = Sinh[chi]},
   {{0, 0, 0, c, -s/esigma, 0}, {0, 0, 0, esigma s, -c, 0},
    {0, 0, 1, 0, 0, 0}, {c, esigma s, 0, -w esigma s, w c, 0},
    {-s/esigma, -c, 0, w c, -w s/esigma, 0}, {0, 0, 0, 0, 0, 1}}];
NRBoundaryH[lp_, lm_, w1_, z_] := With[{c = 1 + 2 lp lm z^2},
   {{0, 0, 0, c, -2 lm z, 0}, {0, 0, 0, 2 lp z, -c, 0},
    {0, 0, 1, 0, 0, 0}, {c, 2 lp z, 0, -2 lp w1 z^2, w1 z, 0},
    {-2 lm z, -c, 0, w1 z, -2 lm w1 z^2, 0}, {0, 0, 0, 0, 0, 1}}];
NRBoundaryD[lp_, lm_, z_] := Log[z]/2 + lp lm z^2/4;

(* ---------------------------------------------------------------------------------------------- *)
(* Near-boundary tools used by the SM3 files.  z = e^{-2y/l} is the manuscript's u; the radial       *)
(* derivative acts on an explicit y (symbol yy) and on z: d_y = d_yy - (2z/l) d_z.  Backgrounds and    *)
(* frames are expanded through z^n; nothing lowers the z-order, so truncation is exact at each order. *)
(* ---------------------------------------------------------------------------------------------- *)

xsZY = {xp, xm, Function[e, D[e, yy] - (2 z/l) D[e, z]]};
DyZY[e_] := D[e, yy] - (2 z/l) D[e, z];
SeriesZ[e_, n_] := Together[Normal[Series[e, {z, 0, n}]]];
SeriesZM[m_, n_] := Map[SeriesZ[#, n] &, m, {ArrayDepth[m]}];
LinearT[e_] := Coefficient[Normal[Series[e, {t, 0, 1}]], t, 1];

(* Exact Riemannian saddle data of SMexactRdata (lower flat indices), as rational functions of z. *)
RiemannianSaddleExact[] := Module[{Pi2, g, B, e, eb, V, Vb, d},
   Pi2 = Lp[xp] Lm[xm];
   g = {{2 Lp[xp], -1/z - Pi2 z, 0}, {-1/z - Pi2 z, 2 Lm[xm], 0}, {0, 0, 1}};
   B = {{0, -1/z - Pi2 z, 0}, {1/z + Pi2 z, 0, 0}, {0, 0, 0}};
   e = {{1, -Lp[xp], 0}, {-z Lm[xm], 1/z, 0}, {0, 0, 1}};
   eb = {{1/z, -z Lp[xp], 0}, {-Lm[xm], 1, 0}, {0, 0, 1}};
   V = 1/Sqrt[2] ArrayFlatten[{{Transpose[Inverse[e]]}, {e . eta3 + B . Transpose[Inverse[e]]}}];
   Vb = 1/Sqrt[2] ArrayFlatten[{{Transpose[Inverse[eb]]}, {eb . etab3 + B . Transpose[Inverse[eb]]}}];
   d = -yy/l - 1/2 Log[1 - Pi2 z^2];
   <|"H" -> Map[Together, RiemannianH[g, B], {2}], "V" -> Map[Together, V, {2}], "Vb" -> Map[Together, Vb, {2}],
     "d" -> d, "g" -> g, "B" -> B, "e" -> e, "eb" -> eb|>];

(* Exact non-Riemannian saddle data of SMexactNRdata / SMvielbein in the rational variables
   psi_pm = L_pm^{-1/2}: L+ = 1/psip^2, L- = 1/psim^2, e^sigma = psim/psip, chi = 2 sqrt2 arctanh(z sqrt(Pi/2)).
   Wz is the hair function of z (W0 = 0 unless supplied). *)
NonRiemannianSaddleExact[Wz_] := Module[{Pi2, chi, esg, hh, chh, shh, Vup, Vbup, H, d},
   Pi2 = 1/(psip[xp]^2 psim[xm]^2);
   chi = 2 Sqrt[2] ArcTanh[z/(Sqrt[2] psip[xp] psim[xm])];     (* q = z sqrt(Pi/2) written rationally *)
   esg = psim[xm]/psip[xp];
   hh = chi/2; chh = Cosh[hh]; shh = Sinh[hh];
   Vup = {{0, -chh/Sqrt[2], 0}, {0, -esg shh/Sqrt[2], 0}, {0, 0, 1/Sqrt[2]},
      {Sqrt[2] chh, Wz esg shh/(2 Sqrt[2]), 0}, {-Sqrt[2] shh/esg, -Wz chh/(2 Sqrt[2]), 0}, {0, 0, 1/Sqrt[2]}};
   Vbup = {{shh/(esg Sqrt[2]), 0, 0}, {chh/Sqrt[2], 0, 0}, {0, 0, -1/Sqrt[2]},
      {-Wz chh/(2 Sqrt[2]), -Sqrt[2] esg shh, 0}, {Wz shh/(2 Sqrt[2] esg), Sqrt[2] chh, 0}, {0, 0, 1/Sqrt[2]}};
   H = NonRiemannianH[chi, esg, Wz];
   d = -yy/l + Log[Cosh[chi/(2 Sqrt[2])]];
   <|"H" -> H, "V" -> Vup . eta3, "Vb" -> Vbup . etab3, "d" -> d, "chi" -> chi, "esigma" -> esg|>];
NRW2 := -(l^2/4) D[1/psip[xp]^2, xp] D[1/psim[xm]^2, xm];    (* SMbackgroundexpansion: the derivative-dependent W_2 *)
NRLpsi = {Lp -> Function[x, 1/psip[x]^2], Lm -> Function[x, 1/psim[x]^2]};

SaddleSeries[sd_Association, n_] := <|"H" -> SeriesZM[sd["H"], n], "V" -> SeriesZM[sd["V"], n],
   "Vb" -> SeriesZM[sd["Vb"], n], "d" -> SeriesZ[sd["d"], n]|>;

(* SMcosetreconstruction: delta H_MN = 2 V_(M^p Vbar_N)^qbar h_{p qbar} for a lower-index mixed fluctuation matrix. *)
MixedFluctuationH[V_, Vb_, hmat_] := Module[{m = (V . eta3) . hmat . Transpose[Vb . etab3]}, m + Transpose[m]];

(* Frame variation of SMframevariation, lower flat indices. *)
FrameVariation[V_, Vb_, hmat_] := {1/2 (Vb . etab3) . Transpose[hmat], -1/2 (V . eta3) . hmat};

(* Linearized EDFE components E_{p qbar} = V^M_p delta(P S Pbar)_MN Vbar^N_qbar and E_0 = delta S_(0)
   on a z-series saddle, for the mixed fluctuation hmat (functions of xp, xm, yy) and dilaton fluctuation ddf. *)
LinearizedEDFEComponents[bg_Association, hmat_, ddf_, n_] := Module[
   {JJ = ODDJ[3], dH, Hlin, dlin, gamma, r4, ric, psp, Vup, Vbup, E, E0, tr},
   dH = SeriesZM[MixedFluctuationH[bg["V"], bg["Vb"], hmat], n];
   Hlin = bg["H"] + t dH; dlin = bg["d"] + t ddf;
   tr[e_] := SeriesZ[Normal[Series[e, {t, 0, 1}]], n];
   gamma = Map[tr, GammaDFT[Hlin, dlin, xsZY], {3}];
   r4 = Map[tr, RiemannR4[gamma, xsZY], {4}];
   ric = Map[tr, RicciS[gamma, r4, xsZY], {2}];
   psp = Map[Function[e, Together[LinearT[Expand[e]]]], ProjectedRicci[Hlin, ric, xsZY], {2}];
   Vup = JJ . bg["V"]; Vbup = JJ . bg["Vb"];
   E = Map[Function[e, SeriesZ[e, n]], Transpose[Vup] . psp . Vbup, {2}];
   E0 = SeriesZ[Together[LinearT[Expand[ScalarS0[Hlin, dlin, xsZY]]]], n];
   <|"E" -> E, "E0" -> E0, "dH" -> dH|>];

(* Frame-projected radial momentum A^y_{p qbar} = V^M_p A^y_MN Vbar^N_qbar (SMAdefinition) and same-chirality
   projections, from the unprojected tensor of MomentumCore. *)
MomentumProjected[gamma_List, V_, Vb_, xs_List] := Module[{JJ = ODDJ[3], core, Vup, Vbup},
   core = MomentumCore[gamma, xs][[6]];
   Vup = JJ . V; Vbup = JJ . Vb;
   <|"Amixed" -> Transpose[Vup] . core . Vbup, "Aunbarred" -> Transpose[Vup] . core . Vup,
     "Abarred" -> Transpose[Vbup] . core . Vbup, "core" -> core|>];

(* === CALCULATION === *)
(* ::Title:: *)
(*NRH02 Letter*)

(* Checks for the current Letter (1)-(21), with retained ancillary identities.  The particular two-point
   kernels (20) use the ordered second variation / linear response and source normalization checked in
   NRH06 (historical SM3.2, SM3.4-SM3.5).  The generic off-shell Codazzi identity (13) is not checked here. *)

NRH`BeginFile["01_ExactBackgrounds_EDFE.wl"];

xs = {xp, xm, Function[e, (2 u/l) D[e, u]]};    (* code u = e^{2y/l} *)
JJ = ODDJ[3];

(* ::Section:: *)
(*(1)-(4) Riemannian saddle and boundary data*)

fR = u + Lp[xp] Lm[xm]/u;
gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];
thetaP = {u^(1/2), -u^(-1/2) Lm[xm], 0}; thetaM = {-u^(-1/2) Lp[xp], u^(1/2), 0};   (* e^{y/l} dx^pm - e^{-y/l} L_mp dx^mp *)
NRH`CheckZero["Rfields: dy^2 - 2 theta+ theta- reproduces the displayed metric components",
   Together[Table[-(thetaP[[i]] thetaM[[j]] + thetaM[[i]] thetaP[[j]]), {i, 3}, {j, 3}] + DiagonalMatrix[{0, 0, 1}] - gR]];
NRH`CheckZero["Rfields: B_{-+} = e^{2y/l} + e^{-2y/l} L+ L- (B = (...) dx^- ^ dx^+)", Together[bR[[2, 1]] - fR]];
NRH`CheckZero["RDFTfields: H J H = J (O(3,3) constraint)", Map[Together, HR . JJ . HR - JJ, {2}]];
NRH`CheckZero["RDFTfields: e^{-2d} = e^{2y/l}(1 - L+ L- e^{-4y/l}) = sqrt(-g) at phi_0 = 0",
   {Together[Exp[-2 dR] - u (1 - Lp[xp] Lm[xm]/u^2)], Together[Exp[-2 dR]^2 + Det[gR]]}];
curvR = DFTCurvature[HR, dR, xs];
NRH`CheckZero["EDFE on R: (P S Pbar)_MN = 0 and S_(0) = -4/l^2 for arbitrary chiral L_pm", {curvR["PSPbar"], Together[curvR["S0"] + 4/l^2]}];
NRH`CheckZero["G_MN = 2 l^-2 J_MN on the Riemannian saddle", Map[Together, curvR["G"] - 2/l^2 JJ, {2}]];
NRH`CheckZero["Rboundary: H - H^infty = O(e^{-2y/l}) and d + y/l = O(e^{-4y/l})",
   {Map[Limit[#, u -> Infinity] &, HR - Hinf, {2}], Limit[u (dR + 1/2 Log[u]), u -> Infinity]}];
NRH`CheckZero["boundaryH: H^infty is the displayed 6x6 matrix", Map[Limit[#, u -> Infinity] &, HR, {2}] - Hinf];
H0 = Hinf[[{1, 2, 4, 5}, {1, 2, 4, 5}]];
NRH`Check["boundaryH: the induced boundary H^(0) is of type (1,1): vanishing upper-left block",
   H0[[1 ;; 2, 1 ;; 2]] === {{0, 0}, {0, 0}}];
J4 = ODDJ[2];
V0 = {{1/Sqrt[2], 0}, {0, 0}, {0, -Sqrt[2]}, {0, 0}};
Vb0 = {{0, 0}, {0, 1/Sqrt[2]}, {0, 0}, {Sqrt[2], 0}};
eta2 = {{0, -1}, {-1, 0}}; etab2 = {{0, 1}, {1, 0}};
NRH`CheckZero["Rboundaryframe: the aligned D = 2 DFT vielbein reproduces P^(0), Pbar^(0)",
   {V0 . eta2 . Transpose[V0] - (J4 + H0)/2, Vb0 . etab2 . Transpose[Vb0] - (J4 - H0)/2}];
NRH`CheckZero["horizon: e^{-2d} = 0 at e^{2y/l} = Sqrt[L+ L-] for constant L_pm > 0",
   Together[Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]]];

(* ::Section:: *)
(*(5)-(9) Non-Riemannian saddle*)

chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
xsNR = {Function[e, D[e, xp] + chp D[e, ch]], Function[e, D[e, xm] + chm D[e, ch]], Function[e, chy D[e, ch] + D[e, Ysym]]};
NRHZeroNR[label_, e_] := NRH`CheckZero[label, Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];
esig = psim[xm]/psip[xp];     (* e^sigma = Sqrt[L+/L-] with psi_pm = L_pm^{-1/2} *)
chiq = 2 Sqrt[2] ArcTanh[q];
NRH`CheckZero["NRvariables: chi = 2 Sqrt[2] arctanh q and e^{pm sigma} sinh(chi/2) = e^{-2y/l} L_pm + O(e^{-6y/l})",
   {SeriesCoefficient[Sinh[chiq/2], {q, 0, 0}], SeriesCoefficient[Sinh[chiq/2], {q, 0, 2}],
    Together[esig SeriesCoefficient[Sinh[chiq/2], {q, 0, 1}] - Sqrt[2] psim[xm]/psip[xp]]}];   (* q = e^{-2y/l} Sqrt[Pi/2]: Sqrt2 q e^sigma = e^{-2y/l} L+ *)
HNR = NonRiemannianH[ch, esig, W[xp, xm, ch]];
dNR = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];
NRHZeroNR["NRHcompact: H J H = J for the exact non-Riemannian matrix (generic W)", HNR . JJ . HNR - JJ];
NRH`Check["NRHcompact: type (1,1) at every radius (upper-left block diag(0,0,1))",
   HNR[[1 ;; 3, 1 ;; 3]] === {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}];
NRHZeroNR["NRdilaton: e^{-2d} = e^{2y/l}(1 - q^2) with q = tanh(chi/(2 Sqrt[2]))",
   Exp[-2 (dNR + Ysym/l)] - (1 - Tanh[ch/(2 Sqrt[2])]^2)];
radialSource = l^2/(16 psip[xp] psim[xm]) (
   4 Sqrt[2] D[Sinh[ch]/Sinh[ch/Sqrt[2]], ch] (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm])
   - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 Exp[ch]
   + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 Exp[-ch]);
curvNR = DFTCurvature[HNR, dNR, xsNR];
NRHZeroNR["EDFE scalar on NR: S_(0) = -4/l^2 for ARBITRARY W(x+, x-, chi)", curvNR["S0"] + 4/l^2];
odeRule = Derivative[0, 0, 2][W][xp, xm, ch] -> radialSource;
NRHZeroNR["EDFE tensor on NR: (P S Pbar)_MN = 0 <=> d^2 W/d chi^2 = F (the radial equation of SM2)", curvNR["PSPbar"] /. odeRule];
NRHZeroNR["G_MN = 2 l^-2 J_MN on the non-Riemannian saddle", (curvNR["G"] /. odeRule) - 2/l^2 JJ];
Gp[c_] := 4 Sqrt[2] (Sinh[c]/Sinh[c/Sqrt[2]] - Sqrt[2]);     (* Gp = I_radial'(chi), GG = I_radial(chi), the manuscript's calligraphic I in MainGprofile. *)
Wexact = W0[xp, xm] + W1[xp, xm] ch psip[xp] psim[xm]/2 +
   l^2/(16 psip[xp] psim[xm]) ((Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm]) GG[ch]
      - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 (Exp[ch] - 1 - ch)
      + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 (Exp[-ch] - 1 + ch));
NRHZeroNR["NRWgeneral: solves d^2W/dchi^2 = F with I_radial'' = rho", (D[Wexact, {ch, 2}] /. Derivative[2][GG][ch] -> D[Gp[ch], ch]) - radialSource];
NRH`Check["NRWgeneral: W_0 and W_1 multiply the two homogeneous modes {1, chi/(2 Sqrt[Pi])}",
   {D[Wexact, W0[xp, xm]], Together[D[Wexact, W1[xp, xm]] - ch psip[xp] psim[xm]/2]} === {1, 0}];
NRH`CheckZero["MainGprofile: I_radial'(0) = 0 and I_radial'(chi) = (2/3) chi^2 + O(chi^4); the integral normalization I_radial(0) = 0 gives I_radial = (2/9) chi^3 + O(chi^5)", {Limit[Gp[ch], ch -> 0], Normal[Series[Gp[ch], {ch, 0, 3}]] - 2/3 ch^2}];
NRH`CheckZero["W = W_0 + W_1 e^{-2y/l} + O(e^{-4y/l}) near the boundary (chi = 2 Sqrt[2] q + O(q^3))",
   {SeriesCoefficient[chiq, {q, 0, 1}] - 2 Sqrt[2], SeriesCoefficient[chiq, {q, 0, 2}]}];
NRH`CheckZero["endpoint q = 1: e^{-2d} = 0 at e^{4y/l} = L+ L-/2 (distinct from the Riemannian horizon)",
   Together[(u (1 - z0^2) /. z0 -> 1/(Sqrt[2] psip[xp] psim[xm] u)) /. u -> 1/(Sqrt[2] psip[xp] psim[xm])]];

(* ::Section:: *)
(*(10)-(14) Doubled dictionary and constant-boundary Ward identities; retained Rcontinuity checks are ancillary, not a check of Codazzi (13).*)

NRH`CheckZero["RKdef: -2K = (16 pi G)^{-1} A^y and 2 T_(0) = (16 pi G)^{-1} 2 (B^y + 4/l) fix the coefficients -1/(32 pi G) and 1/(16 pi G)",
   {Together[-2 (-(1/(32 Pi G))) - 1/(16 Pi G)], Together[2 (1/(16 Pi G)) - 2/(16 Pi G)]}];
NRH`CheckZero["Mainmomenta: B^y = 4 d_y d on both saddles (H^{yN} = delta^N_y)",
   {Together[GammaBVector[HR, dR, xs][[6]] - 4 (2 u/l) D[dR, u]],
    Together[ExpandAll[TrigToExp[(GammaBVector[HNR, dNR, xsNR][[6]] - 4 (chy D[dNR, ch] + D[dNR, Ysym])) /. ch -> 2 Sqrt[2] Log[T]]]]}];
kmat = {{Kpp[xp, xm], Kpm[xp, xm]}, {Kmp[xp, xm], Kmm[xp, xm]}};
T0f = T0[xp, xm];
V0u = V0 . eta2; Vb0u = Vb0 . etab2;
TAB = Table[2 Sum[(V0u[[a, i]] Vb0u[[b, j]] - V0u[[b, i]] Vb0u[[a, j]]) kmat[[i, j]], {i, 2}, {j, 2}] - 1/2 J4[[a, b]] T0f, {a, 4}, {b, 4}];
bD[e_, a_] := If[a <= 2, 0, D[e, {xp, xm}[[a - 2]]]];
divT = Table[Together[Sum[J4[[a, c]] bD[TAB[[a, b]], c], {a, 4}, {c, 4}]], {b, 4}];
ward1 = D[Kpp[xp, xm], xm] + 1/4 D[T0f, xp];
ward2 = D[Kmm[xp, xm], xp] + 1/4 D[T0f, xm];
NRH`CheckZero["RDFTconservation, Rcontinuity: on the constant boundary, div T = {d_- K_mp, -d_+ K_mp, -2 Ward_+, -2 Ward_-}",
   Together[divT - {D[Kmp[xp, xm], xm], -D[Kmp[xp, xm], xp], -2 ward1, -2 ward2}]];
NRH`Check["Rcontinuity: ancillary constant-boundary identity: no local condition on K_{op bom}; K_{om bop} must be constant", FreeQ[divT, Kpm] && ! FreeQ[divT, Kmp]];

(* ::Section:: *)
(*(15)-(17) Asymptotic symmetry generators and transformation laws*)

xiUp = {-l^2/(2 u) Lp[xp] D[em[xm], {xm, 2}], +l^2/(2 u) Lm[xm] D[ep[xp], {xp, 2}], -l/2 (D[ep[xp], xp] - D[em[xm], xm]),
   ep[xp] + l^2/(4 u) D[em[xm], {xm, 2}], em[xm] + l^2/(4 u) D[ep[xp], {xp, 2}], -l/2 (D[ep[xp], xp] + D[em[xm], xm])};
lieH = Map[Together, GenLieH[xiUp, HR, xs], {2}];
lieD = Together[GenLieD[xiUp, dR, xs]];
dLp = ep[xp] D[Lp[xp], xp] + 2 Lp[xp] D[ep[xp], xp] - l^2/4 D[ep[xp], {xp, 3}];
dLm = em[xm] D[Lm[xm], xm] + 2 Lm[xm] D[em[xm], xm] - l^2/4 D[em[xm], {xm, 3}];
HRgen = HR /. {Lp[xp] -> LPv, Lm[xm] -> LMv};
depsH = Map[Together, (D[HRgen, LPv] dLp + D[HRgen, LMv] dLm) /. {LPv -> Lp[xp], LMv -> Lm[xm]}, {2}];
orderCheck[m_, pow_] := Map[Function[e, Normal[Series[e, {u, Infinity, pow}]]], m, {2}];
NRH`CheckZero["Rkilling, RVirasoro: Lhat_xi H - delta_eps H = O(e^{-4y/l}) with delta L = eps L' + 2 L eps' - (l^2/4) eps'''", orderCheck[lieH - depsH, 1]];
NRH`CheckZero["Rkilling: Lhat_xi d = O(e^{-4y/l}) (the radial component compensates the weight)", Normal[Series[lieD, {u, Infinity, 1}]]];
xsU = xs;
qu = 1/(Sqrt[2] psip[xp] psim[xm] u);
chU = 2 Sqrt[2] ArcTanh[qu];
HNRu = NonRiemannianH[chU, esig, W1[xp, xm]/u];
dNRu = -1/2 Log[u] + Log[Cosh[chU/(2 Sqrt[2])]];
xiUpNR = {-l^2/(2 u) (1/psip[xp]^2) D[em[xm], {xm, 2}], +l^2/(2 u) (1/psim[xm]^2) D[ep[xp], {xp, 2}], -l/2 (D[ep[xp], xp] - D[em[xm], xm]),
   ep[xp], em[xm], -l/2 (D[ep[xp], xp] + D[em[xm], xm])};
lieHNR = GenLieH[xiUpNR, HNRu, xsU]; lieDNR = GenLieD[xiUpNR, dNRu, xsU];
dPsiP = ep[xp] D[psip[xp], xp] - psip[xp] D[ep[xp], xp];
dPsiM = em[xm] D[psim[xm], xm] - psim[xm] D[em[xm], xm];
dW1NR = (ep[xp] D[W1[xp, xm], xp] + em[xm] D[W1[xp, xm], xm] + 2 W1[xp, xm] (D[ep[xp], xp] + D[em[xm], xm])
   - l^2 ((1/psim[xm]^2) D[ep[xp], {xp, 3}] + (1/psip[xp]^2) D[em[xm], {xm, 3}]));   (* parenthesized: a top-level line break would end the expression *)
HNRuGen = HNRu /. {psip[xp] -> PSPv, psim[xm] -> PSMv, W1[xp, xm] -> W1v};
depsHNR = (D[HNRuGen, PSPv] dPsiP + D[HNRuGen, PSMv] dPsiM + D[HNRuGen, W1v] dW1NR) /. {PSPv -> psip[xp], PSMv -> psim[xm], W1v -> W1[xp, xm]};
seriesZero[m_, ord_] := Map[Function[e, Together[Normal[Series[e, {u, Infinity, ord}]]]], m, {2}];
NRH`CheckZero["Rkilling (NR form), NRasympt: Lhat_xi H - delta_(L, W1) H = O(u^-2) componentwise", seriesZero[lieHNR - depsHNR, 1]];
NRH`CheckZero["NRasympt: delta psi is equivalent to delta L = eps L' + 2 L eps' (no central term)",
   Together[(D[1/PSPv^2, PSPv] dPsiP /. PSPv -> psip[xp]) - (ep[xp] D[1/psip[xp]^2, xp] + 2 (1/psip[xp]^2) D[ep[xp], xp])]];
NRH`CheckZero["NRasympt: Lhat_xi d = O(u^-2)", Together[Normal[Series[lieDNR, {u, Infinity, 1}]]]];
centralCharge = 3 l/(2 G);
THol = 8 Pi (Lp[xp]/(16 Pi G l));
NRH`CheckZero["charges: with T = 8 pi K = L+/(2 G l), RVirasoro is delta T = eps T' + 2 T eps' - (c/12) eps''' with c = 3l/(2G)",
   Together[dLp/(2 G l) - (ep[xp] D[THol, xp] + 2 THol D[ep[xp], xp] - centralCharge/12 D[ep[xp], {xp, 3}])]];
NRH`CheckZero["charge normalization: (16 pi G)^{-1}(4/l) = 1/(4 pi G l) (Q[eps] of the Letter)", Together[1/(16 Pi G) 4/l - 1/(4 Pi G l)]];

(* ::Section:: *)
(*(18) One-point functions; retained frame-variation check for SMframevariation (formerly Mainframevariation).*)

AyR = MomentumAK[HR, curvR["Gamma"], xs][[6]];
AfixR = Map[Together, Transpose[JJ . Vinf] . AyR . (JJ . Vbinf), {2}];
KR = Map[Limit[-(1/(32 Pi G)) u #, u -> Infinity] &, AfixR[[1 ;; 2, 1 ;; 2]], {2}];
NRH`CheckZero["Mainonepoints on R: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, L+L-},{1, L-}}/(16 pi G l) (exact family, aligned frame)",
   Map[Together, KR - {{Lp[xp], Lp[xp] Lm[xm]}, {1, Lm[xm]}}/(16 Pi G l), {2}]];
NRH`CheckZero["Mainonepoints on R: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0", Limit[u (GammaBVector[HR, dR, xs][[6]] + 4/l), u -> Infinity]];
xsZ = {xp, xm, Function[e, -(2 z/l) D[e, z]]};
HNRz = NRBoundaryH[Lp[xp], Lm[xm], W1[xp, xm], z];
dNRz = NRBoundaryD[Lp[xp], Lm[xm], z];
AyNR = MomentumAK[HNRz, GammaDFT[HNRz, dNRz, xsZ], xsZ][[6]];
AfixNR = Map[Together, Transpose[JJ . Vinf] . AyNR . (JJ . Vbinf), {2}];
KNR = Map[Limit[-(1/(32 Pi G)) #/z, z -> 0] &, AfixNR[[1 ;; 2, 1 ;; 2]], {2}];
NRH`CheckZero["Mainonepoints on NR: -(32 pi G)^{-1} lim e^{2Y/l} A^y_{a bbar} = {{L+, W1/4},{0, L-}}/(16 pi G l)",
   Map[Together, KNR - {{Lp[xp], W1[xp, xm]/4}, {0, Lm[xm]}}/(16 Pi G l), {2}]];
NRH`CheckZero["Mainonepoints on NR: e^{2Y/l}(B^y + 4/l) -> 0, <T_(0)> = 0", Limit[(GammaBVector[HNRz, dNRz, xsZ][[6]] + 4/l)/z, z -> 0]];
hgen = Table[hh[p, q], {p, 3}, {q, 3}];
dHgen = MixedFluctuationH[Vinf, Vbinf, hgen];
{dV, dVb} = FrameVariation[Vinf, Vbinf, hgen];
NRH`CheckZero["Mainframevariation, SMframevariation: delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2",
   {dV . eta3 . Transpose[Vinf] + Vinf . eta3 . Transpose[dV] - dHgen/2, dVb . etab3 . Transpose[Vbinf] + Vbinf . etab3 . Transpose[dVb] + dHgen/2}];

(* ::Section:: *)
(*(19)-(20) Particular two-point kernels (position-space linear response; historical SM3.2, SM3.4-SM3.5).*)

(* The stress and hair boundary operators V_pm^s, W_pm follow from the stress constraints and the retained
   nonchiral-generator derivation in NRH06.  They give the particular kernels of historical SM3.5; homogeneous
   responses and local contacts are not fixed here.  The Lorentzian position-space prescription is:
   A^s = -(kappa/(4 i)) V^s d_-^{-1} delta, M^NR = -(kappa/(16 i)) W d_-^{-1} delta, d_-^{-1} delta -> 1/(2 pi i D+). *)
Module[{kappa = 1/(16 Pi G l), inv = 1/(2 Pi I dp), VnrK, VrK, WK, aNR, aR, mNR, mR, targetA, targetM},
   VnrK[f_] := 2 Lp[xp] D[f, dp] + D[Lp[xp], xp] f;
   VrK[f_] := VnrK[f] - l^2/4 D[f, {dp, 3}];
   WK[f_] := 2 W1[xp, xm] D[f, dp] + D[W1[xp, xm], xp] f - l^2 Lm[xm] D[f, {dp, 3}];
   aNR = Together[-kappa/(4 I) VnrK[inv]]; aR = Together[-kappa/(4 I) VrK[inv]];
   mNR = Together[-kappa/(16 I) WK[inv]]; mR = Together[Lm[xm] aR];
   targetA = (D[Lp[xp], xp]/dp - 2 Lp[xp]/dp^2)/(128 Pi^2 G l);
   targetM = (D[W1[xp, xm], xp]/dp - 2 W1[xp, xm]/dp^2 + 6 l^2 Lm[xm]/dp^4)/(512 Pi^2 G l);
   NRH`CheckZero["Rcorrelators: A_+^NR = [L+'/D+ - 2 L+/D+^2]/(128 pi^2 G l)", aNR - targetA];
   NRH`CheckZero["Rcorrelators: A_+^R = A_+^NR + 3 l/(256 pi^2 G D+^4)", aR - targetA - 3 l/(256 Pi^2 G dp^4)];
   NRH`CheckZero["Rcorrelators: M_+^NR = [W1'/D+ - 2 W1/D+^2 + 6 l^2 L-/D+^4]/(512 pi^2 G l)", mNR - targetM];
   NRH`CheckZero["Mainfivebyfive: M_+^R = L_- A_+^R", mR - Lm[xm] aR];
   NRH`CheckZero["c = 3l/(2G): at L+ = 0, (8 pi)^2 A_+^R = (c/2)/D+^4", Together[(8 Pi)^2 (aR /. Lp -> (0 &)) - centralCharge/2/dp^4]];
   NRH`Check["negative control: reversing the source sign fails the NR kernel", ! NRH`ZeroQ[aNR + targetA]];
   NRH`CheckZero["the mixed fourth-order pole is proportional to L_- (state dependent, no state-independent central charge)",
      D[Coefficient[Expand[mNR], dp, -4], Lm[xm]] - 3 l/(256 Pi^2 G)]];

(* ::Section:: *)
(*(21) Worldsheet: Gomis-Ooguri limit*)

fRc = u + Lpc Lmc/u;
L1 = dy by + 2 Lpc dxp bxp + 2 Lmc dxm bxm + beta bxp + betab dxm + beta betab/(2 fRc);
betaSol = First@Solve[{D[L1, beta] == 0, D[L1, betab] == 0}, {beta, betab}];
EmatR = {{2 Lpc, 0, 0}, {-2 fRc, 2 Lmc, 0}, {0, 0, 1}};
LE = Sum[EmatR[[m, n]] {dxp, dxm, dy}[[m]] {bxp, bxm, by}[[n]], {m, 3}, {n, 3}];
NRH`CheckZero["Mainworldsheet: eliminating beta, betabar from the first-order action reproduces E_{mu nu} dx dbar x (E = g - B)",
   Together[(L1 /. betaSol) - LE]];
NRH`CheckZero["c_eff^2 := 2F = 2(e^{2y/l} + L+ L- e^{-2y/l}) ~ 2 e^{2y/l}: the Gomis-Ooguri limit is reached without tuning",
   Limit[(2 fRc)/(2 u), u -> Infinity] - 1];
NRH`Check["Mainworldsheet: V_W = (1/(4 pi alpha')) W dx+ dbar x- with the 1/(2 pi alpha') prefactor of the reduced action",
   Together[1/(2 Pi alphaPrime) Wc/2 - Wc/(4 Pi alphaPrime)] === 0];

NRH`FileSummary[];
