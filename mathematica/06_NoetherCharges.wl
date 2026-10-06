(* 06_NoetherCharges.wl | 2026-10-06 standalone edition.
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
(*NRH07 SM3 Covariant Charges and Asymptotic Algebras*)

(* Historical SM3.6: surface-charge one-form, Riemannian normalization, non-Riemannian charges,
   C-bracket and cocycles. Assertion identifiers are retained from the expanded archive.
   In the R calculation below, u = exp(+2y/l), reciprocal to the manuscript radial variable;
   in the NR calculation z = exp(-2y/l), equal to the manuscript u. Their derivatives and
   boundary limits accordingly differ. Antisymmetric K^(y+) is minus the printed K^(+y).
   Mixed parameter-bracket leftovers are tested by their vanishing surface potentials,
   not by claiming that the raw parameter bracket vanishes. The final L+T/4 identity is
   ancillary algebra and does not alter the physical NR charge central term. *)

NRH`BeginFile["06_NoetherCharges.wl"];

JJ = ODDJ[3];

NoetherK[HH_, xUp_, a_, b_, xs_] := Module[
   {n = Length[xs], dim, Hup, HfirstUp, HsecondUp, xDown, val},
   dim = 2 n;
   Hup = JJ . HH . JJ; HfirstUp = JJ . HH; HsecondUp = HH . JJ;
   xDown = JJ . xUp;
   val = 0;
   Do[
      val -= Hup[[c, a]] (DblD[xUp[[b]], c, xs] + Sum[JJ[[b, f]] DblD[xDown[[c]], f, xs], {f, dim}]);
      val += Hup[[c, b]] (DblD[xUp[[a]], c, xs] + Sum[JJ[[a, f]] DblD[xDown[[c]], f, xs], {f, dim}]);
      Do[
         val -= (Hup[[c, a]] Hup[[b, dd]] - Hup[[c, b]] Hup[[a, dd]]) DblD[HH[[dd, e]], c, xs] xUp[[e]];
         val -= 1/2 HsecondUp[[e, c]] (HfirstUp[[a, dd]] DblD[Hup[[b, dd]], c, xs]
              - HfirstUp[[b, dd]] DblD[Hup[[a, dd]], c, xs]) xUp[[e]],
         {dd, dim}, {e, dim}],
      {c, dim}];
   Do[
      val += Sum[JJ[[a, f]] DblD[HfirstUp[[b, e]], f, xs]
           - JJ[[b, f]] DblD[HfirstUp[[a, e]], f, xs], {f, dim}] xUp[[e]],
      {e, dim}];
   val];

KhatComp[HH_, dd_, xUp_, a_, b_, xs_] := Module[{bv = GammaBVector[HH, dd, xs]},
   NoetherK[HH, xUp, a, b, xs] + xUp[[a]] bv[[b]] - xUp[[b]] bv[[a]]];

xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};

gR = RiemannianMetric[Lp[xp], Lm[xm], u];
bR = RiemannianB[Lp[xp], Lm[xm], u];
HR = Map[Together, RiemannianH[gR, bR], {2}];
dR = RiemannianD[Lp[xp], Lm[xm], u];

xiPlus = {0, l^2/(2 u) Lm[xm] D[ep[xp], {xp, 2}], -l/2 D[ep[xp], xp],
   ep[xp], l^2/(4 u) D[ep[xp], {xp, 2}], -l/2 D[ep[xp], xp]};
xiMinus = {-l^2/(2 u) Lp[xp] D[em[xm], {xm, 2}], 0, +l/2 D[em[xm], xm],
   l^2/(4 u) D[em[xm], {xm, 2}], em[xm], -l/2 D[em[xm], xm]};

KfullP = Together[Exp[-2 dR] KhatComp[HR, dR, xiPlus, 5, 6, xsU]];
NRH`CheckZero["lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''",
   Together[Limit[KfullP, u -> Infinity]
      - (4/l ep[xp] Lp[xp] - 2 l D[ep[xp], {xp, 2}])]];

KfullM = Together[Exp[-2 dR] KhatComp[HR, dR, xiMinus, 6, 4, xsU]];
NRH`CheckZero["lim e^{-2d} Khat^{y+}[eps-] = -(4/l) eps- L- + 2 l eps-'' (mirror)",
   Together[Limit[KfullM, u -> Infinity]
      - (-(4/l) em[xm] Lm[xm] + 2 l D[em[xm], {xm, 2}])]];

deltaL[e_] := e D[Lp[xp], xp] + 2 Lp[xp] D[e, xp] - l^2/4 D[e, {xp, 3}];
alpha12 = e1[xp] D[e2[xp], xp] - e2[xp] D[e1[xp], xp];
cocycle = Together[e1[xp] deltaL[e2[xp]] - alpha12 Lp[xp]
   + l^2/4 e1[xp] D[e2[xp], {xp, 3}]];
ELx[f_, e_] := Together[D[e, f[xp]] - D[D[e, Derivative[1][f][xp]], xp]
   + D[D[e, Derivative[2][f][xp]], {xp, 2}] - D[D[e, Derivative[3][f][xp]], {xp, 3}]];
NRH`CheckZero["Virasoro cocycle: (charge bracket density) + (l^2/4) e1 e2''' is a total derivative",
   {ELx[Lp, cocycle], ELx[e1, cocycle], ELx[e2, cocycle]}];
NRH`Check["the central term itself is NOT a total derivative (the center is real)",
   ! NRH`ZeroQ[ELx[e1, e1[xp] D[e2[xp], {xp, 3}]]]];

xsZ = {xp, xm, Function[e, -(2 z/l) D[e, z]]};
HNRz = NRBoundaryH[Lp[xp], Lm[xm], W1[xp, xm], z];
dNRz = NRBoundaryD[Lp[xp], Lm[xm], z];
eDenz = (1 - Lp[xp] Lm[xm] z^2/2)/z;

NRH`Check["truncation obeys H J H = J through z^2",
   Module[{c = Expand[HNRz . JJ . HNRz - JJ]},
      AllTrue[Flatten[c], PossibleZeroQ[Coefficient[#, z, 0]] && PossibleZeroQ[Coefficient[#, z, 1]] && PossibleZeroQ[Coefficient[#, z, 2]] &]]];
NRH`CheckZero["state-dependent falloffs delta H^-_+ = 2 z dL+, delta H^+_- = -2 z dL-, delta H_{+-} = z dW1",
   {D[HNRz[[2, 4]], Lp[xp]] - 2 z, D[HNRz[[1, 5]], Lm[xm]] + 2 z, D[HNRz[[4, 5]], W1[xp, xm]] - z}];

varyRules = {Lp -> Function[x, Lp[x] + tt dLpF[x]], Lm -> Function[x, Lm[x] + tt dLmF[x]],
   W1 -> Function[{x, y2}, W1[x, y2] + tt dW1F[x, y2]]};
HNRzT = HNRz /. varyRules; dNRzT = dNRz /. varyRules; eDenzT = eDenz /. varyRules;
dH = D[HNRzT, tt] /. tt -> 0; dd0 = D[dNRzT, tt] /. tt -> 0;

gammaZ = GammaDFT[HNRz, dNRz, xsZ];
HupZ = JJ . HNRz . JJ; dHup = JJ . dH . JJ;
ThetaHat = Table[
   Module[{val},
      val = 4 Sum[HupZ[[a, b]] DblD[dd0, b, xsZ], {b, 6}];
      Do[
         val -= DblD[dHup[[a, b]], b, xsZ];
         Do[val += gammaZ[[b, c, f]] JJ[[f, a]] dHup[[c, b]]
             + gammaZ[[b, c, f]] JJ[[f, b]] dHup[[a, c]], {c, 6}, {f, 6}],
         {b, 6}];
      Together[eDenz val - (D[eDenzT (GammaBVector[HNRzT, dNRzT, xsZ][[a]]), tt] /. tt -> 0)]],
   {a, 4, 6}];
NRH`CheckZero["lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary",
   Map[Limit[#, z -> 0] &, ThetaHat]];

chargeOneForm[xiOf_, aa_, bb_] := Module[{xi, xiT, varied, fieldDep, thetaTerm},
   xi = xiOf[Lp[xp], Lm[xm]];
   xiT = xiOf[Lp[xp] + tt dLpF[xp], Lm[xm] + tt dLmF[xm]];
   varied = D[eDenzT KhatComp[HNRzT, dNRzT, xiT, aa, bb, xsZ], tt] /. tt -> 0;
   fieldDep = eDenz KhatComp[HNRz, dNRz, D[xiT, tt] /. tt -> 0, aa, bb, xsZ];
   thetaTerm = xi[[aa]] ThetaHat[[bb - 3]] - xi[[bb]] ThetaHat[[aa - 3]];
   Limit[Together[varied - fieldDep + thetaTerm], z -> 0]];

xiP = Function[{lp, lm}, {0, l^2 z lm D[ep[xp], {xp, 2}]/2, -l D[ep[xp], xp]/2, ep[xp], 0, -l D[ep[xp], xp]/2}];
xiM = Function[{lp, lm}, {-l^2 z lp D[em[xm], {xm, 2}]/2, 0, +l D[em[xm], xm]/2, 0, em[xm], -l D[em[xm], xm]/2}];

kPlus = chargeOneForm[xiP, 5, 6];
kMinus = chargeOneForm[xiM, 4, 6];
NRH`CheckZero["k^{-y}[eps+] = (4/l) eps+ dL+",
   Together[kPlus - 4/l ep[xp] dLpF[xp]]];
NRH`CheckZero["k^{+y}[eps-] = (4/l) eps- dL-",
   Together[kMinus - 4/l em[xm] dLmF[xm]]];
NRH`Check["W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise",
   FreeQ[{kPlus, kMinus}, W1] && FreeQ[{kPlus, kMinus}, dW1F] &&
   FreeQ[kPlus, dLmF] && FreeQ[kMinus, dLpF]];

CBracket[x_, y_, xs_] := Module[{xd = JJ . x, yd = JJ . y, dim = 6},
   Table[
      Sum[x[[b]] DblD[y[[a]], b, xs] - y[[b]] DblD[x[[a]], b, xs], {b, dim}]
      + 1/2 Sum[yd[[b]] Sum[JJ[[a, c]] DblD[x[[b]], c, xs], {c, dim}]
              - xd[[b]] Sum[JJ[[a, c]] DblD[y[[b]], c, xs], {c, dim}], {b, dim}],
      {a, dim}]];

xiPe = Function[{e}, {0, l^2 z Lm[xm] D[e, {xp, 2}]/2, -l D[e, xp]/2, e, 0, -l D[e, xp]/2}];
alphaP = e1[xp] D[e2[xp], xp] - e2[xp] D[e1[xp], xp];
bracketDiff = Together[CBracket[xiPe[e1[xp]], xiPe[e2[xp]], xsZ] - xiPe[alphaP]];
NRH`Check["the same-chirality C-bracket closes up to a closed B-gauge parameter (slot x~+ only)",
   Together[bracketDiff[[2 ;; 6]]] === {0, 0, 0, 0, 0} && ! PossibleZeroQ[bracketDiff[[1]]]];
NRH`CheckZero["the leftover reducibility parameter is chiral and closed: d_- and d_y of it vanish",
   {D[bracketDiff[[1]], xm], D[bracketDiff[[1]], z]}];
NRH`CheckZero["the closed B-gauge parameter carries no surface potential",
   Limit[Together[eDenz KhatComp[HNRz, dNRz, {zp[xp], 0, 0, 0, 0, 0}, 5, 6, xsZ]], z -> 0]];

NRH`CheckZero["NR cocycle e1 (e2 L' + 2 L e2') - alpha L = d/dx (e1 e2 L) => c_charge = 0",
   Together[e1[xp] (e2[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e2[xp], xp])
      - alphaP Lp[xp] - D[e1[xp] e2[xp] Lp[xp], xp]]];

NRH`CheckZero["(i) k^{-y}[eps+] = delta[(4/l) eps+ L+] (the charge exists and is integrable)",
   Together[kPlus - D[4/l ep[xp] LQ, LQ] dLpF[xp]]];
NRH`CheckZero["(i) mirror: k^{+y}[eps-] = delta[(4/l) eps- L-]",
   Together[kMinus - D[4/l em[xm] LQ, LQ] dLmF[xm]]];

adjNR[a_, b_, x_] := a[x] D[b[x], x] - b[x] D[a[x], x];
brPP = Together[(kPlus /. ep -> e1f) /.
   dLpF -> Function[x, e2f[x] Derivative[1][Lp][x] + 2 Lp[x] Derivative[1][e2f][x]]];
NRH`CheckZero["(ii) {Q[e1+], Q[e2+]} - Q[[e1,e2]] is a total derivative: centerless plus sector",
   {ELx[Lp, Together[brPP - 4/l adjNR[e1f, e2f, xp] Lp[xp]]],
    ELx[e1f, Together[brPP - 4/l adjNR[e1f, e2f, xp] Lp[xp]]],
    ELx[e2f, Together[brPP - 4/l adjNR[e1f, e2f, xp] Lp[xp]]]}];
brMM = Together[(kMinus /. em -> e1g) /.
   dLmF -> Function[x, e2g[x] Derivative[1][Lm][x] + 2 Lm[x] Derivative[1][e2g][x]]];
ELm[f_, e_] := Together[D[e, f[xm]] - D[D[e, Derivative[1][f][xm]], xm]
   + D[D[e, Derivative[2][f][xm]], {xm, 2}] - D[D[e, Derivative[3][f][xm]], {xm, 3}]];
NRH`CheckZero["(ii) minus-sector mirror: {Q[e1-], Q[e2-]} - Q[[e1,e2]] is a total derivative",
   {ELm[Lm, Together[brMM - 4/l adjNR[e1g, e2g, xm] Lm[xm]]],
    ELm[e1g, Together[brMM - 4/l adjNR[e1g, e2g, xm] Lm[xm]]],
    ELm[e2g, Together[brMM - 4/l adjNR[e1g, e2g, xm] Lm[xm]]]}];
NRH`CheckZero["(ii) opposite chiralities Poisson-commute: delta_{eps-} L+ = 0 kills the mixed bracket",
   {kPlus /. dLpF -> (0 &), kMinus /. dLmF -> (0 &)}];

brR = Together[4/l e1f[xp] (e2f[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e2f[xp], xp]
      - l^2/4 D[e2f[xp], {xp, 3}])];
centralDensity = -4/l l^2/4 e1f[xp] D[e2f[xp], {xp, 3}];
NRH`CheckZero["(iii) R bracket - adjoint - central = total derivative (Brown-Henneaux center isolated)",
   {ELx[Lp, Together[brR - 4/l adjNR[e1f, e2f, xp] Lp[xp] - centralDensity]],
    ELx[e1f, Together[brR - 4/l adjNR[e1f, e2f, xp] Lp[xp] - centralDensity]],
    ELx[e2f, Together[brR - 4/l adjNR[e1f, e2f, xp] Lp[xp] - centralDensity]]}];
NRH`CheckZero["(iii) normalization chain: (16 pi G)^{-1} (4/l) = 1/(4 pi G l), the Letter's charge normalization",
   Together[1/(16 Pi G) 4/l - 1/(4 Pi G l)]];

NRH`CheckZero["(iv) the central cocycle is antisymmetric modulo total derivatives",
   {ELx[e1f, Together[e1f[xp] D[e2f[xp], {xp, 3}] + e2f[xp] D[e1f[xp], {xp, 3}]]],
    ELx[e2f, Together[e1f[xp] D[e2f[xp], {xp, 3}] + e2f[xp] D[e1f[xp], {xp, 3}]]]}];
NRH`CheckZero["(iv) Witt Jacobi identity: [[e1,e2],e3] + cyclic = 0 exactly",
   Module[{br = Function[{a, b}, a D[b, xp] - b D[a, xp]]},
      Together[br[br[e1f[xp], e2f[xp]], e3f[xp]] + br[br[e2f[xp], e3f[xp]], e1f[xp]]
         + br[br[e3f[xp], e1f[xp]], e2f[xp]]]]];
NRH`CheckZero["(iv) Gelfand-Fuchs cocycle condition: c(e1,[e2,e3]) + cyclic = total derivative",
   Module[{cc = Function[{a, b}, a D[b, {xp, 3}]], br = Function[{a, b}, a D[b, xp] - b D[a, xp]],
      jj, el4},
      el4[f_, e_] := Together[D[e, f[xp]] - D[D[e, Derivative[1][f][xp]], xp]
         + D[D[e, Derivative[2][f][xp]], {xp, 2}] - D[D[e, Derivative[3][f][xp]], {xp, 3}]
         + D[D[e, Derivative[4][f][xp]], {xp, 4}]];
      jj = Together[cc[e1f[xp], br[e2f[xp], e3f[xp]]] + cc[e2f[xp], br[e3f[xp], e1f[xp]]]
         + cc[e3f[xp], br[e1f[xp], e2f[xp]]]];
      {el4[e1f, jj], el4[e2f, jj], el4[e3f, jj]}]];

xiMeOf[e_, LPval_] := {-l^2 z LPval D[e, {xm, 2}]/2, 0, +l D[e, xm]/2, 0, e, -l D[e, xm]/2};
dpLp = e1f[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e1f[xp], xp];
dmLm = e2g[xm] D[Lm[xm], xm] + 2 Lm[xm] D[e2g[xm], xm];
deltaPXm = D[xiMeOf[e2g[xm], Lp[xp] + tt dpLp], tt] /. tt -> 0;
deltaMXp = D[(xiPe[e1f[xp]] /. Lm[xm] -> Lm[xm] + tt dmLm), tt] /. tt -> 0;
mixedAdj = Together[CBracket[xiPe[e1f[xp]], xiMeOf[e2g[xm], Lp[xp]], xsZ] - deltaPXm + deltaMXp];
NRH`Check["(v) the adjusted mixed bracket has a leftover (the raw closure fails, as it should)",
   ! NRH`ZeroQ[mixedAdj]];
NRH`CheckZero["(v) but its Khat surface potentials vanish at the boundary: mixed charge bracket = 0",
   {Limit[Together[eDenz KhatComp[HNRz, dNRz, mixedAdj, 5, 6, xsZ]], z -> 0],
    Limit[Together[eDenz KhatComp[HNRz, dNRz, mixedAdj, 4, 6, xsZ]], z -> 0]}];

alphaM = e1g[xm] D[e2g[xm], xm] - e2g[xm] D[e1g[xm], xm];
bracketDiffM = Together[CBracket[xiMeOf[e1g[xm], Lp[xp]], xiMeOf[e2g[xm], Lp[xp]], xsZ]
   - xiMeOf[alphaM, Lp[xp]]];
NRH`Check["(vi) minus-sector C-bracket closes up to a closed B-gauge parameter (slot x~- only)",
   Together[bracketDiffM[[{1, 3, 4, 5, 6}]]] === {0, 0, 0, 0, 0} && ! PossibleZeroQ[bracketDiffM[[2]]]];
NRH`CheckZero["(vi) that leftover is chiral and closed, and carries no surface potential",
   {D[bracketDiffM[[2]], xp], D[bracketDiffM[[2]], z],
    Limit[Together[eDenz KhatComp[HNRz, dNRz, {0, zm[xm], 0, 0, 0, 0}, 4, 6, xsZ]], z -> 0]}];

NRH`CheckZero["footnote: historical auxiliary algebra only: with delta T = eps T' + 2 T eps' - l^2 eps''' the combination L + T/4 obeys the law with -(l^2/4) eps'''",
   Module[{dL = e1[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e1[xp], xp],
      dT = e1[xp] D[TT[xp], xp] + 2 TT[xp] D[e1[xp], xp] - l^2 D[e1[xp], {xp, 3}], comb},
      comb = Lp[xp] + TT[xp]/4;
      Together[dL + dT/4 - (e1[xp] D[comb, xp] + 2 comb D[e1[xp], xp] - l^2/4 D[e1[xp], {xp, 3}])]]];

(* Integrate the actual surface one-form along L(s)=s L, W1(s)=s W1.
   The periodic spatial circle has dx+=l/sqrt2 dphi, dx-=-l/sqrt2 dphi. *)
qPlus=Integrate[kPlus/.dLpF[xp]->Lp[xp],{path,0,1}]/(16 Pi G);
qMinus=Integrate[kMinus/.dLmF[xm]->Lm[xm],{path,0,1}]/(16 Pi G);
NRH`CheckZero["integrated Q coefficients from k, with background subtraction",
 {qPlus-ep[xp] Lp[xp]/(4 Pi G l),qMinus-em[xm] Lm[xm]/(4 Pi G l)}];
qDensity=FullSimplify[(l/Sqrt[2])(qPlus+qMinus)];
Print["Q[epsilon] = Integral(dphi,0,2pi) of ",InputForm[qDensity]," with x_pm=(t +/- l phi)/sqrt2."];
qConstant=Integrate[qDensity/.{ep[xp]->ep0,em[xm]->em0,Lp[xp]->lp0,Lm[xm]->lm0},{angle,0,2 Pi}];
NRH`CheckZero["constant-data circle charge",qConstant-(ep0 lp0+em0 lm0)/(2 Sqrt[2] G)];
Print["Constant-data Q = ",InputForm[qConstant]];

NRH`FileSummary[];
