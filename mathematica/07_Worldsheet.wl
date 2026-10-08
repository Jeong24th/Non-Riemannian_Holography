(* 07_Worldsheet.wl | 2026-10-08 standalone edition.
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
(*SM3 Lorentzian Worldsheet Reduction and Radial Vertex Operators*)

NRH`BeginFile["07_Worldsheet.wl"];

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
NRH`CheckZero["SMRfirstorder: multiplier shifts retain the full L_pm-dependent finite-radius interaction",
   Together[(L1 /. {beta -> betaNew - 2 Lp dxp, betab -> betabNew - 2 Lm bxm})
    - (dy by + betaNew bxp + betabNew dxm
       + (betaNew - 2 Lp dxp) (betabNew - 2 Lm bxm)/(2 fR))]];
NRH`CheckZero["-det g_par = F^2 - 4 L+ L- = e^{-4 d_R}",
   Together[-Det[{{2 Lp, -fR}, {-fR, 2 Lm}}] - (fR^2 - 4 Lp Lm)]];
NRH`CheckZero["c_eff^2 = 2F -> 4 Sqrt[L+L-] at the horizon u = Sqrt[L+L-]",
   Together[(2 fR /. u -> Sqrt[Lp Lm]) - 4 Sqrt[Lp Lm]]];

(* Internal legacy identifiers tauP/tauM denote the manuscript's X/bar-X. *)
tauP = {Cosh[ch/2], -esig^-1 Sinh[ch/2], 0};
tauM = {-esig Sinh[ch/2], Cosh[ch/2], 0};
Hupper = {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}};
NRH`CheckZero["SNCtau: H^{mu nu} X_nu = H^{mu nu} bar-X_nu = 0",
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
NRH`CheckZero["symmetric-block route - antisymmetric-clock route = W (bar-X . dx)(X . dbar x)",
   Together[symRoute - antisymRoute - Wc mdotd tdotb]];

NRH`CheckZero["at chi -> 0 the W coupling reduces to (W/2) dx^+ dbar x^- (+ constraint terms)",
   Together[(symRoute /. ch -> 0 /. esig -> 1) - Wc/2 (dxp + 0) (bxm + 0) - Wc/2 dxm bxp]];
NRH`Check["with the 1/(2 pi alpha') prefactor this is V_W = (1/(4 pi alpha')) W dx+ dbar x-",
   Together[1/(2 Pi alphaPrime) Wc/2 - Wc/(4 Pi alphaPrime)] === 0];

(* Ancillary historical calculation: the classical wound-string energy passage
   and SMlongstringE were deleted from the manuscript on 2026-10-08.
   The following seven algebraic checks are retained as historical evidence,
   not coverage of a current manuscript claim. *)
Print["Ancillary historical checks: deleted classical wound-string energy passage."];
gR = RiemannianMetric[Lp, Lm, u];
bR = RiemannianB[Lp, Lm, u];
et = {1/Sqrt[2], 1/Sqrt[2], 0};
ephi = {l wN/Sqrt[2], -l wN/Sqrt[2], 0};
g2 = {{et . gR . et, et . gR . ephi}, {ephi . gR . et, ephi . gR . ephi}};

ephi1 = {l/Sqrt[2], -l/Sqrt[2], 0};
gtphi = {{et . gR . et, et . gR . ephi1}, {ephi1 . gR . et, ephi1 . gR . ephi1}};
Btphi = et . bR . ephi1;
NRH`Check["Ancillary deleted-energy passage: constant L_pm gives t-independent (t,phi) metric and B components",
   FreeQ[{gtphi, Btphi}, t] && FreeQ[{gR, bR}, xp] && FreeQ[{gR, bR}, xm]];
NRH`CheckZero["Ancillary deleted-energy passage: gamma^{tau a} g_{t nu} d_a X^nu = 1 on the static embedding",
   Together[Sum[Inverse[g2][[1, a]] (et . gR . {et, ephi}[[a]]), {a, 2}] - 1]];
NRH`CheckZero["Ancillary deleted-energy passage: B_{t phi} = l F and B_{t nu} X'^nu = w B_{t phi}",
   {Together[Btphi - l (u + Lp Lm/u)], Together[et . bR . ephi - wN Btphi]}];
sqrtMinusGamma = l wN (u - Lp Lm/u);
Pt = 1/(2 Pi alphaPrime) (-sqrtMinusGamma Sum[Inverse[g2][[1, a]] (et . gR . {et, ephi}[[a]]), {a, 2}] + et . bR . ephi);
NRH`CheckZero["Ancillary deleted-energy passage: E = -Int P_t reproduces E(y) = -(2 w l/alpha') L+L- e^{-2y/l}",
   Together[-2 Pi Pt + 2 wN l/alphaPrime Lp Lm/u]];
NRH`CheckZero["Ancillary deleted-energy passage: -det g_(t,phi) equals the squared Nambu-Goto area density",
   Together[-Det[g2] - (l wN (u - Lp Lm/u))^2]];

BtphiPerW = l (u + Lp Lm/u);
NRH`CheckZero["Ancillary deleted-energy passage: area minus B coupling gives E(y)=-(2 w l/alpha') L+L- e^{-2y/l}",
   Together[wN l/alphaPrime ((u - Lp Lm/u) - (u + Lp Lm/u)) + 2 wN l/alphaPrime Lp Lm/u]];
NRH`CheckZero["Ancillary deleted-energy passage: with phi_0=0 the area density equals l e^{-2d}",
   Together[l (u - Lp Lm/u) - l u (1 - Lp Lm/u^2)]];
Print["End ancillary historical checks; resume current SM3 calculations."];

fluxH = Integrate[D[l^2 Cos[th]^2, th], {th, 0, Pi/2}] (2 Pi) (2 Pi);
NRH`CheckZero["k = |Int_{S^3} H| / (4 pi^2 alpha') = l^2/alpha'",
   Together[(-fluxH)/(4 Pi^2 alphaPrime) - l^2/alphaPrime]];

(* One real Lorentzian null separation; this is one chiral component of
   <y(s)y(u)> = -(alpha'/2) Log[(s+ - u+)(s- - u-)]. The causal boundary
   value is an input and is not established by this separated-point algebra. *)
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
NRH`Check["V_W: the undressed W_0 longitudinal operator has no self-contraction pair inside {x^+, x^-}^2; radial W_1 dressing is checked separately",
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

NRH`CheckZero["SMBRSTfusion: each real null separation has exponent -2 alpha'/l^2 for a = -2/l",
   Together[-alphaPrime/2 (-2/l)^2 + 2 alphaPrime/l^2]];
hn = nn + qq nn (1 - nn);
NRH`CheckZero["the n-fold fused weight h_n = n + q n(1-n) equals n + h_y(-2n/l) with q = alpha'/l^2",
   Together[(hn - (nn + (-alphaPrime/4 (-2 nn/l) (-2 nn/l + 2/l)))) /. qq -> alphaPrime/l^2]];
NRH`Check["h_n = 1 exactly at n = 1 or at the resonant value q = 1/n",
   Solve[hn == 1, qq] === {{qq -> 1/nn}} && Together[(hn /. nn -> 1) - 1] === 0];

(* SNCdualB and SNCreconstruction with current X/bar-X notation; exact frames. *)
Module[{cc=chiLocal,ee=sigmaExp,ww=wLocal,tp,tm,yv,yb,bmat,lower,hh,omega,v,vb},
 tp={Cosh[cc/2],-Sinh[cc/2]/ee,0}; tm={-ee Sinh[cc/2],Cosh[cc/2],0};
 yv={Cosh[cc/2],ee Sinh[cc/2],0}; yb={Sinh[cc/2]/ee,Cosh[cc/2],0};
 bmat=-ww/2(Outer[Times,tp,tm]-Outer[Times,tm,tp]);
 NRH`CheckZero["SNCdualB: dual vectors and X wedge bar-X = dx+ wedge dx-",
 {tp.yv-1,tm.yb-1,tp.yb,tm.yv,bmat-{{0,-ww/2,0},{ww/2,0,0},{0,0,0}}}];
 lower=DiagonalMatrix[{0,0,1}]+Outer[Times,tp,bmat.yv]+Outer[Times,bmat.yv,tp]
 -Outer[Times,tm,bmat.yb]-Outer[Times,bmat.yb,tm];
 hh=NonRiemannianH[cc,ee,ww];
 NRH`CheckZero["SNCreconstruction: all lower-right generalized-metric components",lower-hh[[4;;6,4;;6]]];
 omega=ArrayFlatten[{{IdentityMatrix[3],ConstantArray[0,{3,3}]},{bmat,IdentityMatrix[3]}}];
 v[w_] := {{0,-Cosh[cc/2]/Sqrt[2],0},{0,-ee Sinh[cc/2]/Sqrt[2],0},{0,0,1/Sqrt[2]},
 {Sqrt[2] Cosh[cc/2],w ee Sinh[cc/2]/(2 Sqrt[2]),0},{-Sqrt[2] Sinh[cc/2]/ee,-w Cosh[cc/2]/(2 Sqrt[2]),0},{0,0,1/Sqrt[2]}};
 vb[w_] := {{Sinh[cc/2]/(ee Sqrt[2]),0,0},{Cosh[cc/2]/Sqrt[2],0,0},{0,0,-1/Sqrt[2]},
 {-w Cosh[cc/2]/(2 Sqrt[2]),-Sqrt[2] ee Sinh[cc/2],0},{w Sinh[cc/2]/(2 Sqrt[2] ee),Sqrt[2] Cosh[cc/2],0},{0,0,1/Sqrt[2]}};
 NRH`CheckZero["SMvielbein: both frames are the B transforms of W=0 frames",{omega.v[0]-v[ww],omega.vb[0]-vb[ww]}];
];

(* Direct Lorentzian reduction of SMdyg before longitudinal constraints.
   This block omits the overall 1/(4 pi alpha'). Set rt=exp(chi/2), so the
   rational ct,st expressions implement cosh(chi/2),sinh(chi/2) exactly. *)
Module[{rt, ee, wc, ct, st, xcov, xbar, yvec, ybar, transverse, mixed, lower,
   hh, up, um, uy, vp, vm, vy, ap, am, ay, abp, abm, aby, du, dv, au, av,
   ddu, ddv, raw, transverseRules, betaOld, betabarOld, expected, be, beb,
   newbe, multiplierLag, shifted, nullJac, eps01, flatH, stressP, stressM},
 ct = (rt + 1/rt)/2; st = (rt - 1/rt)/2;
 xcov = {ct, -st/ee, 0}; xbar = {-ee st, ct, 0};
 yvec = {ct, ee st, 0}; ybar = {st/ee, ct, 0};
 transverse = DiagonalMatrix[{0, 0, 1}];
 mixed = Outer[Times, yvec, xcov] - Outer[Times, ybar, xbar];
 lower = transverse + wc (Outer[Times, xcov, xbar] + Outer[Times, xbar, xcov]);
 hh = ArrayFlatten[{{transverse, mixed}, {Transpose[mixed], lower}}];
 NRH`CheckZero["SNCdualB: completeness Y X + bar-Y bar-X + transverse = identity",
   Outer[Times, yvec, xcov] + Outer[Times, ybar, xbar] + transverse - IdentityMatrix[3]];
 NRH`CheckZero["SNCreconstruction: complete X/bar-X generalized metric obeys H J H = J",
   hh . ODDJ[3] . hh - ODDJ[3]];
 nullJac = {{1, 1}, {1, -1}}/Sqrt[2]; eps01 = {{0, 1}, {-1, 0}};
 NRH`CheckZero["SMdyg Lorentzian orientation: epsilon^{+-} = -1",
   nullJac . eps01 . Transpose[nullJac] - {{0, -1}, {1, 0}}];
 du = {up, um, uy}; dv = {vp, vm, vy}; au = {ap, am, ay}; av = {abp, abm, aby};
 ddu = Join[-au, du]; ddv = Join[-av, dv];
 raw = ddu . hh . ddv + du . av - dv . au;
 transverseRules = {ay -> -uy, aby -> vy};
 NRH`CheckZero["SMphysicalsectionA: transverse auxiliary equations a_y=-partial y and bar-a_y=bar-partial y",
   {D[raw, ay], D[raw, aby]} /. transverseRules];
 betaOld = -yvec . au; betabarOld = ybar . av;
 expected = 2 uy vy + 2 betaOld (xcov . dv) + 2 betabarOld (xbar . du)
   + wc ((xcov . du) (xbar . dv) + (xbar . du) (xcov . dv));
 NRH`CheckZero["SMdygGO: exact Lorentzian doubled-action reduction at finite chi",
   (raw /. transverseRules) - expected];
 NRH`CheckZero["SMdygGO: no longitudinal auxiliary bilinear at any radius",
   Table[D[raw /. transverseRules, aa, bb], {aa, {ap, am}}, {bb, {abp, abm}}]];
 multiplierLag = 2 uy vy + 2 be (xcov . dv) + 2 beb (xbar . du)
   + wc ((xcov . du) (xbar . dv) + (xbar . du) (xcov . dv));
 NRH`CheckZero["SMdygconstraints: longitudinal multipliers impose X dot bar-partial x and bar-X dot partial x",
   {D[multiplierLag, be] - 2 xcov . dv, D[multiplierLag, beb] - 2 xbar . du}];
 shifted = multiplierLag /. be -> newbe - wc/2 (xbar . du);
 NRH`CheckZero["SMdygGO: unit-Jacobian beta shift leaves exactly W (X dot partial x)(bar-X dot bar-partial x)",
   shifted - (2 uy vy + 2 newbe (xcov . dv) + 2 beb (xbar . du)
     + wc (xcov . du) (xbar . dv))];
 NRH`CheckZero["SMdygGO: beta redefinition has unit derivative",D[newbe - wc/2 (xbar . du), newbe] - 1];
 flatH = hh /. {rt -> 1, wc -> 0};
 stressP = (ddu . flatH . ddu /. transverseRules) /. um -> 0;
 stressM = (ddv . flatH . ddv /. transverseRules) /. vp -> 0;
 NRH`CheckZero["SMGO: doubled metric variation gives both normalized Lorentzian chiral stresses",
   {stressP - 2 (uy^2 - ap up), stressM - 2 (vy^2 + abm vm)}];
];

Module[{pp, pm, zz, halfchi, ee, xx, xb},
 halfchi = Sqrt[2] ArcTanh[zz/(Sqrt[2] pp pm)]; ee = pm/pp;
 xx = {Cosh[halfchi], -Sinh[halfchi]/ee, 0};
 xb = {-ee Sinh[halfchi], Cosh[halfchi], 0};
 NRH`CheckZero["SNCtau: first X/bar-X corrections are -z L_- dx^- and -z L_+ dx^+",
   Normal[Series[{xx, xb}, {zz, 0, 1}]] - {{1, -zz/pm^2, 0}, {-zz/pp^2, 1, 0}}];
];

(* W0 has no radial dressing. Contract beta(s) with W0(x(u)) and partial x(u),
   then expand partial x(s). The right-moving calculation is identical. *)
Module[{ss, uu, alpha, wc, wx, gamma1, gamma2, barred1, delta, contraction,
   operator, derivative},
 delta = ss - uu;
 contraction = (-1/alpha) (gamma1 + delta gamma2) barred1
   (wx (-alpha/delta) gamma1 + wc D[-alpha/delta, uu]);
 operator = wc gamma1 barred1;
 derivative = (wx gamma1^2 + wc gamma2) barred1;
 NRH`CheckZero["SMvertex: Lorentzian T V_W0 OPE has unit double-pole weight and the derivative simple pole",
   Normal[Series[contraction - operator/delta^2 - derivative/delta, {ss, uu, -1}]]];
];

NRH`CheckZero["SMBRSTcentral: bosonic level k_bos=k+2 gives the same central charge",
   3 (kk + 2)/((kk + 2) - 2) - (3 + 6/kk)];

NRH`FileSummary[];
