(* 04_GeneralBackground_LEDFE.wl | 2026-10-06 standalone edition.
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
(*NRH05 SM3 Linearized Dynamics and Holographic Renormalization*)

(* Historical SM3: 3.1 constrained variations, 3.2 action variations, 3.3 linearized EDFE and interior
   condition, 3.4 Riemannian branch, 3.5 non-Riemannian branch, 3.6 covariant charges.
   This file retains expanded component algebra underlying historical SM3.1 and SM3.3-3.5;
   stable assertion identifiers also refer to ancillary displays removed in the compact rewrite.
   General linearized curvature equations are checked through u. The current integrated
   a/b/omega coefficient formulas require a separate equivalence bridge; the new fixed-flux
   K3 interior prescription and H_R matching are not tested here. The Einstein-tensor assembly
   assumes its defining formula and does not prove the universal Box or Codazzi identities.
   z = e^{-2y/l} (manuscript u); yy is the explicit y. Two-point action checks are in NRH06. *)

NRH`BeginFile["04_GeneralBackground_LEDFE.wl"];

JJ = ODDJ[3];

(* ::Section:: *)
(*SM3.1 Constrained field variations*)

hgen = Table[hh[p, q], {p, 3}, {q, 3}];
Pinf = (JJ + Hinf)/2; Pbinf = (JJ - Hinf)/2;
NRH`CheckZero["SMprojectors: P + Pbar = J, P J P = P, Pbar J Pbar = Pbar, P J Pbar = 0 on the limiting metric",
   {Pinf + Pbinf - JJ, Pinf . JJ . Pinf - Pinf, Pbinf . JJ . Pbinf - Pbinf, Pinf . JJ . Pbinf}];
dHgen = MixedFluctuationH[Vinf, Vbinf, hgen];
NRH`CheckZero["SMcosetreconstruction: P delta H P = 0 = Pbar delta H Pbar for the reconstructed mixed variation",
   {Pinf . JJ . dHgen . JJ . Pinf, Pbinf . JJ . dHgen . JJ . Pbinf}];
NRH`CheckZero["SMmixedfluctuation: projecting the reconstruction back gives h_{p qbar}",
   Transpose[JJ . Vinf] . dHgen . (JJ . Vbinf) - hgen];
NRH`CheckZero["constraint: delta(H J H) = 0 for the mixed variation (H J delta H + delta H J H = 0)",
   Hinf . JJ . dHgen + dHgen . JJ . Hinf];
{dV, dVb} = FrameVariation[Vinf, Vbinf, hgen];
NRH`CheckZero["SMframevariation: delta(V eta V^T) = delta H/2 and delta(Vbar etabar Vbar^T) = -delta H/2",
   {dV . eta3 . Transpose[Vinf] + Vinf . eta3 . Transpose[dV] - dHgen/2,
    dVb . etab3 . Transpose[Vbinf] + Vbinf . etab3 . Transpose[dVb] + dHgen/2}];
NRH`Check["coset count: 9 mixed components - (3 diffeomorphisms + 2 B-gauge) = 4 tangential fields",
   3*3 - (3 + 2) == 4];

(* SMbackgroundconnection: torsionless, compatible, dilaton-trace connection on both general saddles (through z^2) *)
sdR = SaddleSeries[RiemannianSaddleExact[], 2];
sdNR = SaddleSeries[NonRiemannianSaddleExact[W1[xp, xm] z + NRW2 z^2], 2];
Do[
   Module[{bg = data[[2]], gamma, P, compat, trv},
      gamma = Map[SeriesZ[#, 2] &, GammaDFT[bg["H"], bg["d"], xsZY], {3}];
      P = (JJ + bg["H"])/2;
      compat = Table[SeriesZ[DblD[P, c, xsZY][[a, b]]
          + Sum[(gamma[[c]] . JJ)[[a, dd]] P[[dd, b]], {dd, 6}] + Sum[(gamma[[c]] . JJ)[[b, dd]] P[[a, dd]], {dd, 6}], 2],
         {c, 6}, {a, 6}, {b, 6}];
      NRH`CheckZero["SMbackgroundconnection on " <> data[[1]] <> ": nabla_C P_AB = 0 through z^2", compat];
      trv = Table[Sum[JJ[[b, e]] gamma[[e, b, a]], {b, 6}, {e, 6}], {a, 6}];
      NRH`CheckZero["SMbackgroundconnection on " <> data[[1]] <> ": Gamma^B_BA = -2 d_A d through z^2",
         Map[SeriesZ[#, 2] &, trv + 2 DblGrad[bg["d"], xsZY]]];
      NRH`CheckZero["SMbackgroundconnection on " <> data[[1]] <> ": torsionless cyclic sum",
         Map[SeriesZ[#, 2] &, Table[gamma[[c, a, b]] + gamma[[a, b, c]] + gamma[[b, c, a]], {c, 6}, {a, 6}, {b, 6}], {3}]]],
   {data, {{"R", sdR}, {"NR", sdNR}}}];

(* ::Section:: *)
(* Historical SM3 source conventions and SM3.4-3.5 backgrounds; expanded ancillary checks *)

NRH`CheckZero["SMflatmetrics, SMinfinityvielbein: V eta V^T = P^infty, Vbar etabar Vbar^T = Pbar^infty, V^T J Vbar = 0",
   {Vinf . eta3 . Transpose[Vinf] - Pinf, Vbinf . etab3 . Transpose[Vbinf] - Pbinf, Transpose[JJ . Vinf] . Vbinf}];
eta2 = eta3[[1 ;; 2, 1 ;; 2]]; etab2 = etab3[[1 ;; 2, 1 ;; 2]];
hlow2 = {{h0pp, h0pm}, {h0mp, h0mm}};
Kmat = {{Kpp, Kpm}, {Kmp, Kmm}};
NRH`CheckZero["SMsources: -h^{(0)a bbar} K_{a bbar} = h_mm K_pp + h_pp K_mm + h_mp K_pm + h_pm K_mp",
   -Sum[(eta2 . hlow2 . etab2)[[a, b]] Kmat[[a, b]], {a, 2}, {b, 2}] - (h0mm Kpp + h0pp Kmm + h0mp Kpm + h0pm Kmp)];

(* SMW0completion, SMWshift: the closed B shift transports the saddle to arbitrary W_0 *)
bshift[bpm_] := {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
   {0, bpm, 0, 1, 0, 0}, {-bpm, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};
nrEx = NonRiemannianSaddleExact[Wf[xp, xm, z]];
NRH`CheckZero["SMWshift: Omega_b H(W) Omega_b^T = H(W - 2 b_{+-}) exactly",
   Simplify[bshift[bpm] . nrEx["H"] . Transpose[bshift[bpm]] - NonRiemannianSaddleExact[Wf[xp, xm, z] - 2 bpm]["H"]]];
NRH`CheckZero["SMW0completion: b_{+-} = -W_0/2 gives H_s[W_0] with W -> W + W_0, and Omega_b is O(3,3)",
   {Simplify[bshift[-W0[xp, xm]/2] . nrEx["H"] . Transpose[bshift[-W0[xp, xm]/2]] - NonRiemannianSaddleExact[Wf[xp, xm, z] + W0[xp, xm]]["H"]],
    bshift[bpm] . JJ . Transpose[bshift[bpm]] - JJ}];
Module[{Om = bshift[-W0[xp, xm]/2], Vt, Vbt, Ht},
   Vt = Om . nrEx["V"]; Vbt = Om . nrEx["Vb"]; Ht = Om . nrEx["H"] . Transpose[Om];
   NRH`CheckZero["SMW0completion: the transported frames V_s[W_0] = Omega_b V_s[0] reconstruct P, Pbar of H_s[W_0]",
      Simplify[{Vt . eta3 . Transpose[Vt] - (JJ + Ht)/2, Vbt . etab3 . Transpose[Vbt] - (JJ - Ht)/2}]];
   NRH`CheckZero["SMW0completion: flat components h_{p qbar} are unchanged by the transport",
      Simplify[Transpose[JJ . Vt] . (Om . dHgen . Transpose[Om]) . (JJ . Vbt) - Transpose[JJ . nrEx["V"]] . dHgen . (JJ . nrEx["Vb"])]]];
NRH`CheckZero["SMW0variation: varying W_0 gives h^{(0)}_{op bom} = delta W_0/2 at the reference boundary",
   Module[{dOm = D[bshift[-(W0[xp, xm] + s dW0[xp, xm])/2], s], dH},
      dH = dOm . Hinf . Transpose[bshift[-W0[xp, xm]/2]] + bshift[-W0[xp, xm]/2] . Hinf . Transpose[dOm];
      (Transpose[JJ . bshift[-W0[xp, xm]/2] . Vinf] . dH . (JJ . bshift[-W0[xp, xm]/2] . Vbinf))[[1, 2]] - dW0[xp, xm]/2]];

(* SMbackgroundexpansion *)
rEx = RiemannianSaddleExact[];
NRH`CheckZero["SMbackgroundexpansion: d_R = -y/l - (1/2) ln(1 - Pi u^2) and d_NR = -y/l - (1/2) ln(1 - Pi u^2/2)",
   {Simplify[rEx["d"] + yy/l + 1/2 Log[1 - Lp[xp] Lm[xm] z^2]],
    Simplify[TrigToExp[nrEx["d"]] + yy/l + 1/2 Log[1 - z^2/(2 psip[xp]^2 psim[xm]^2)], psip[xp] > 0 && psim[xm] > 0 && z > 0 && z < psip[xp] psim[xm]]}];
NRH`CheckZero["SMbackgroundexpansion: e^{-2 d_s} = u^{-1}[1 + O(u^2)] on both saddles",
   {SeriesCoefficient[Exp[-2 rEx["d"]], {z, 0, 0}] - Exp[2 yy/l] , SeriesCoefficient[Exp[-2 rEx["d"]], {z, 0, -1}],
    SeriesCoefficient[Exp[-2 nrEx["d"]], {z, 0, 0}] - Exp[2 yy/l], SeriesCoefficient[Exp[-2 nrEx["d"]], {z, 0, -1}]} /. Exp[2 yy/l] -> 0];
(* the derivative-dependent hair coefficient W_2 from the exact solution NRWgeneral *)
Module[{chiZ = 2 Sqrt[2] ArcTanh[z/(Sqrt[2] psip[xp] psim[xm])], GG, Wex, c2},
   GG[c_] := 4 Sqrt[2] (Sinh[c]/Sinh[c/Sqrt[2]] - Sqrt[2]);   (* derivative of the manuscript radial profile I(chi), now written calligraphically *)
   Wex = W1[xp, xm] chiZ psip[xp] psim[xm]/2 + l^2/(16 psip[xp] psim[xm]) (
      (Derivative[2][psip][xp] psip[xp] + Derivative[2][psim][xm] psim[xm]) (Integrate[Normal[Series[GG[c], {c, 0, 7}]], c] /. c -> chiZ)   (* calligraphic I(chi) through chi^8 suffices for O(z^2) *)
      - 2 (Derivative[1][psip][xp] + Derivative[1][psim][xm])^2 (Exp[chiZ] - 1 - chiZ)
      + 2 (Derivative[1][psip][xp] - Derivative[1][psim][xm])^2 (Exp[-chiZ] - 1 + chiZ));
   c2 = SeriesCoefficient[Wex, {z, 0, 2}];
   NRH`CheckZero["SMbackgroundexpansion: W_NR = W_1 u - (l^2/4) L+' L-' u^2 + O(u^3) from the exact NRWgeneral",
      {SeriesCoefficient[Wex, {z, 0, 1}] - W1[xp, xm], Simplify[c2 - NRW2]}]];

(* ::Section:: *)
(* Historical SM3.3: direct linearized curvature through u; ancillary radial constraints *)

NRH`CheckZero["SMexactRdata: e eta e^T = g and ebar etabar ebar^T = -g; V, Vbar reconstruct P, Pbar and are orthogonal",
   Simplify[{rEx["e"] . eta3 . Transpose[rEx["e"]] - rEx["g"], rEx["eb"] . etab3 . Transpose[rEx["eb"]] + rEx["g"],
     rEx["V"] . eta3 . Transpose[rEx["V"]] - (JJ + rEx["H"])/2, rEx["Vb"] . etab3 . Transpose[rEx["Vb"]] - (JJ - rEx["H"])/2,
     Transpose[JJ . rEx["V"]] . rEx["Vb"]}]];
NRH`CheckZero["SMexactNRdata: the lowered SMvielbein frames reconstruct P, Pbar of the exact hairy metric and are orthogonal",
   Simplify[{nrEx["V"] . eta3 . Transpose[nrEx["V"]] - (JJ + nrEx["H"])/2, nrEx["Vb"] . etab3 . Transpose[nrEx["Vb"]] - (JJ - nrEx["H"])/2,
     Transpose[JJ . nrEx["V"]] . nrEx["Vb"]}]];
sdR1 = SaddleSeries[rEx, 1];
sdNR1 = SaddleSeries[NonRiemannianSaddleExact[W1[xp, xm] z + NRW2 z^2], 1];
NRH`CheckZero["both saddle frames tend to the limiting frame SMinfinityvielbein at u = 0",
   {(sdR1["V"] /. z -> 0) - Vinf, (sdR1["Vb"] /. z -> 0) - Vbinf, (sdNR1["V"] /. z -> 0) - Vinf, (sdNR1["Vb"] /. z -> 0) - Vbinf}];

(* SMFG, SMexactcomponentrecipe: radial-gauge fluctuation and direct linearized curvature, expanded through u; legacy labels retained. *)
hmat = {{hpp[xp, xm, yy], hpm[xp, xm, yy], 0}, {hmp[xp, xm, yy], hmm[xp, xm, yy], 0}, {0, 0, 0}};
ddf = dd[xp, xm, yy];
Print["  computing the linearized EDFE on the general Riemannian saddle through u ..."];
{tR, linR} = AbsoluteTiming[LinearizedEDFEComponents[sdR1, hmat, ddf, 1]];
Print["  [timing] R: ", Round[tR, 0.1], " s"];
Print["  computing the linearized EDFE on the general non-Riemannian saddle through u ..."];
{tNR, linNR} = AbsoluteTiming[LinearizedEDFEComponents[sdNR1, hmat, ddf, 1]];
Print["  [timing] NR: ", Round[tNR, 0.1], " s"];
compMap = <|"pp" -> {1, 1}, "pm" -> {1, 2}, "py" -> {1, 3}, "mp" -> {2, 1}, "mm" -> {2, 2}, "my" -> {2, 3},
   "yp" -> {3, 1}, "ym" -> {3, 2}, "yy" -> {3, 3}|>;
comps = {"pp", "pm", "py", "mp", "mm", "my", "yp", "ym", "yy", "s0"};
getE[res_, c_] := If[c == "s0", res["E0"], res["E"][[Sequence @@ compMap[c]]]];
E0R[c_] := Coefficient[getE[linR, c], z, 0]; E1R[c_] := Coefficient[getE[linR, c], z, 1];
E0NR[c_] := Coefficient[getE[linNR, c], z, 0]; E1NR[c_] := Coefficient[getE[linNR, c], z, 1];

(* SMexactnormalform: principal radial coefficient of the five evolution components *)
prin[e_, f_] := Coefficient[e, Derivative[0, 0, 2][f][xp, xm, yy]];
NRH`CheckZero["SMexactnormalform: M = diag(-1/4,-1/4,-1/4,-1/4,1) multiplies d_y^2 (h_pp,h_pm,h_mp,h_mm,delta d) in (E_pp,E_pm,E_mp,E_mm,E_yy)",
   {prin[E0R["pp"], hpp] + 1/4, prin[E0R["pm"], hpm] + 1/4, prin[E0R["mp"], hmp] + 1/4, prin[E0R["mm"], hmm] + 1/4, prin[E0R["yy"], dd] - 1}];
NRH`Check["SMexactCauchyconstraints: E_{p ybar}, E_{y qbar} and E_0 - 4 E_{y ybar} contain no second radial derivative (both saddles, through u)",
   FreeQ[Expand[{getE[linR, "py"], getE[linR, "my"], getE[linR, "yp"], getE[linR, "ym"], getE[linR, "s0"] - 4 getE[linR, "yy"],
      getE[linNR, "py"], getE[linNR, "my"], getE[linNR, "yp"], getE[linNR, "ym"], getE[linNR, "s0"] - 4 getE[linNR, "yy"]}],
      Derivative[_, _, 2][_][__]]];
(* SMexactconstraintpropagation: algebraic radial projection using the Einstein-tensor definition; this does not independently prove the Bianchi or Codazzi identity. *)
Module[{Ecomp, dPSP, dG, Es, tang},
   Ecomp = Table[ee[p, q], {p, 3}, {q, 3}];
   dPSP = (sdR1["V"] . eta3) . Ecomp . Transpose[sdR1["Vb"] . etab3];   (* delta(P S Pbar)_MN = V_M^p E_{p qbar} Vbar_N^qbar, lower-index frames *)
   dG = 4 (dPSP - Transpose[dPSP])/2 - 1/2 JJ e0s;
   Es = SeriesZM[dG[[3]], 1];                                       (* row tilde y *)
   tang = Table[Sqrt[2] Sum[(sdR1["Vb"] . etab3)[[N, a]] Ecomp[[3, a]] + (sdR1["V"] . eta3)[[N, a]] Ecomp[[a, 3]], {a, 3}], {N, 6}];
   NRH`CheckZero["SMexactconstraintpropagation: delta G_{y~ N}: tangential = sqrt2(Vbar_N^a E_{y a} + V_N^a E_{a y}), N=y gives -(E_0 - 4E_{yy})/2, N=y~ vanishes",
      SeriesZM[{Es[[{1, 2, 4, 5}]] - tang[[{1, 2, 4, 5}]], Es[[6]] + 1/2 (e0s - 4 ee[3, 3]), Es[[3]]}, 1]]];
(* SMexactTaylorSolution: historical toy recursion with two fields and one boundary variable, not an exact DFT radial solution. *)
Module[{Mm, A0, A1, B0, B1, h0v, h1v, hser, Eop, coeffs, recur},
   Mm = DiagonalMatrix[{-1/4, 1}];
   A0 = Table[a0f[i, j][x], {i, 2}, {j, 2}] + tt Table[a1f[i, j][x], {i, 2}, {j, 2}];       (* A_alpha(t) through t^1 *)
   A1 = Table[c0f[i, j][x], {i, 2}, {j, 2}] + tt Table[c1f[i, j][x], {i, 2}, {j, 2}];       (* coefficient of d_x *)
   B0 = Table[b0f[i, j][x], {i, 2}, {j, 2}] + tt Table[b1f[i, j][x], {i, 2}, {j, 2}];       (* coefficient of d_y *)
   h0v = {f0[x], g0[x]}; h1v = {f1[x], g1[x]};
   recur[n_, hs_] := -(Inverse[Mm] . Sum[
      Coefficient[A0, tt, m] . hs[[n - m + 1]] + Coefficient[A1, tt, m] . D[hs[[n - m + 1]], x]
      + (n - m + 1) Coefficient[B0, tt, m] . hs[[n - m + 2]], {m, 0, Min[n, 1]}])/((n + 2) (n + 1));
   hser = {h0v, h1v}; hser = Append[hser, recur[0, hser]]; hser = Append[hser, recur[1, hser]];
   Eop = Mm . D[Sum[hser[[k + 1]] tt^k, {k, 0, 3}], {tt, 2}] + B0 . D[Sum[hser[[k + 1]] tt^k, {k, 0, 3}], tt]
      + A0 . Sum[hser[[k + 1]] tt^k, {k, 0, 3}] + A1 . D[Sum[hser[[k + 1]] tt^k, {k, 0, 3}], x];
   NRH`CheckZero["SMexactTaylorSolution: the historical toy recursion makes the t^0 and t^1 coefficients of its normal form vanish",
      {Coefficient[Expand[Eop], tt, 0], Coefficient[Expand[Eop], tt, 1]}]];

(* ::Section:: *)
(* Historical SM3.4-3.5: expanded near-boundary component solutions with general sources *)

(* Historical expanded component operators, retained as ancillary checks underlying historical SM3.4-3.5. *)
(* Historical machine-transcribed targets; these are not a direct transcription of the current compact formulas. *)
(* Fields: hpp,hpm,hmp,hmm,dd of [xp,xm,yy]; sources a0,b0,r0,c0,v0 of [xp,xm].
   r0 = h_mp^(0) is the type-changing source; c0 = h_pm^(0) is the W0/B-source channel.
   Rp,Rm are the full h_pp^(2,0),h_mm^(2,0), not just the compact b coefficients.
   H_s is the full h_pm^(2,0) response; the R integration constant c_R is -4 f_0. *)
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

NRH`CheckZero["SMNRradialintegration: (d_y^2 + 2 d_y/l) annihilates 1 and u, gives 2/l on y and -2u/l on u y",
   {DyZY[DyZY[1]] + 2/l DyZY[1], DyZY[DyZY[z]] + 2/l DyZY[z], DyZY[DyZY[yy]] + 2/l DyZY[yy] - 2/l,
    DyZY[DyZY[z yy]] + 2/l DyZY[z yy] + 2 z/l, DyZY[DyZY[yy^2]] + 2/l DyZY[yy^2] - 2 - 4 yy/l}];
fieldsH = {hpp, hpm, hmp, hmm, dd};
actOn[op_, fac_, dyOp_] := Expand[op /. {Derivative[i_, j_, k_][f_][xp, xm, yy] /; MemberQ[fieldsH, f] :> Nest[dyOp, D[fac f[xp, xm, yy], {xp, i}, {xm, j}], k],
   f_[xp, xm, yy] /; MemberQ[fieldsH, f] :> fac f[xp, xm, yy]}];
NRH`CheckZero["SMtotalorderhierarchy: E_I[u^n f] = u^n E_I(d_y - 2n/l, d_+, d_-) f for the ten E^{(0)}_I and all five fields (n = 1, 2)",
   Flatten[Table[actOn[E0R[c], z^n, DyZY] - z^n actOn[E0R[c], 1, (D[#, yy] - 2 n/l #) &], {c, comps}, {n, 1, 2}]]];

W1toR = {W1 -> Function[{a, b}, 4 Lp[a] Lm[b]]};      (* "replace W_1 by 4 L+ L-", including its derivatives *)
okE0 = AllTrue[comps, Together[E0R[#] - E0target[#]] === 0 && Together[E0NR[#] - E0target[#]] === 0 &];
NRH`Check["SMEzerocomponents: the ten historical expanded leading operators E^{(0)} equal the u^0 part of the linearized EDFE on both saddles", okE0];
okE1NR = AllTrue[comps, Together[E1NR[#] - (E1NRtarget[#] /. NRLpsi)] === 0 &];
NRH`Check["SMEoneNRcomponents: the ten historical expanded E^{(1),NR} equal the u^1 part of the linearized EDFE on the general NR saddle", okE1NR];
okDE1 = AllTrue[comps, Together[E1R[#] - (E1NRtarget[#] /. W1toR) - DE1target[#]] === 0 &];
NRH`Check["SMRoperatorrule, SMDeltaEonecomponents: E^{(1),R} = E^{(1),NR}|_{W1 -> 4 L+ L-} + Delta E^{(1)} with the ten displayed Delta E^{(1)}", okDE1];
If[! okE1NR, Do[Print["   NR residual ", c, ": ", Together[E1NR[c] - (E1NRtarget[c] /. NRLpsi)]], {c, comps}]];
If[! okDE1, Do[Print["   R residual ", c, ": ", Together[E1R[c] - (E1NRtarget[c] /. W1toR) - DE1target[c]]], {c, comps}]];
NRH`Check["E^{(1),NR}_{om bop} = 0 while Delta E^{(1)}_{om bop} = -2 d_y delta d/l (the branch dependence of SMlogexample)",
   Together[E1NR["mp"]] === 0 && Together[E1R["mp"] + 2 D[dd[xp, xm, yy], yy]/l] === 0];

(* ---- solve the combined hierarchy SMcoupledorders order by order (SMradialansatz) ---- *)
W1toR = {W1 -> Function[{a, b}, 4 Lp[a] Lm[b]]};      (* "replace W_1 by 4 L+ L-", including its derivatives *)
E1op[s_, c_] := If[s == "NR", E1NRtarget[c], (E1NRtarget[c] /. W1toR) + DE1target[c]];
src = <|"pp" -> b0[xp, xm], "pm" -> c0[xp, xm], "mp" -> r0[xp, xm], "mm" -> a0[xp, xm]|>;
chans = {"pp", "pm", "mp", "mm"};
unk1 = Flatten[{Table[{h1s[c], h1bs[c]}, {c, chans}], dd1s}];
ansatz0[c_] := src[c] + yy h1s[c] + yy^2 h1bs[c];
ans0dd = v0[xp, xm] + yy dd1s;
(* the unknown coefficients are functions of x; represent them by symbols whose x-derivatives are new symbols *)
unkFun[c_] := {h1s[c] -> h1f[c][xp, xm], h1bs[c] -> h1bf[c][xp, xm]};
funRules = Join[Flatten[unkFun /@ chans], {dd1s -> dd1f[xp, xm]}];
subs0[s_] := {hpp -> Function[{xa, xb, xc}, Evaluate[ansatz0["pp"] /. funRules /. {xp -> xa, xm -> xb, yy -> xc}]],
   hpm -> Function[{xa, xb, xc}, Evaluate[ansatz0["pm"] /. funRules /. {xp -> xa, xm -> xb, yy -> xc}]],
   hmp -> Function[{xa, xb, xc}, Evaluate[ansatz0["mp"] /. funRules /. {xp -> xa, xm -> xb, yy -> xc}]],
   hmm -> Function[{xa, xb, xc}, Evaluate[ansatz0["mm"] /. funRules /. {xp -> xa, xm -> xb, yy -> xc}]],
   dd -> Function[{xa, xb, xc}, Evaluate[ans0dd /. funRules /. {xp -> xa, xm -> xb, yy -> xc}]]};
applyOp[op_, rules_] := Expand[op /. Derivative[i_, j_, k_][f_][xp, xm, yy] :> Nest[DyZY, D[f[xp, xm, yy] /. rules, {xp, i}, {xm, j}], k]
   /. f_[xp, xm, yy] /; MemberQ[{hpp, hpm, hmp, hmm, dd}, f] :> (f[xp, xm, yy] /. rules)];
eqs0 = Flatten[Table[CoefficientList[applyOp[E0target[c], subs0["R"]], yy], {c, comps}]];
eqs0 = DeleteCases[Together /@ eqs0, 0];
unkSyms0 = Join[Flatten[Table[{h1f[c][xp, xm], h1bf[c][xp, xm]}, {c, chans}]], {dd1f[xp, xm]}];
(* sequential solution: equations algebraic in an unknown are solved first, then substituted *)
solveSeq[eqsIn_, unks_] := Module[{eqs = eqsIn, sol = {}, remaining = unks, progress = True, alg, s, rule},
   While[progress && remaining =!= {},
      progress = False;
      Do[
         (* an equation is algebraic in uu if it contains uu, no derivative of uu, and no other unsolved unknown (or derivative thereof) *)
         alg = Select[eqs, FreeQ[#, Derivative[__][Head[uu]][__]] && ! FreeQ[#, uu] && FreeQ[#, Alternatives @@ (Head /@ DeleteCases[remaining, uu])] &];
         If[alg =!= {},
            s = Solve[First[alg] == 0, uu];
            If[s =!= {} && s =!= {{}},
               rule = Head[uu] -> Function[{a, b}, Evaluate[(uu /. First[s]) /. {xp -> a, xm -> b}]];
               sol = Join[sol /. rule, {rule}];
               eqs = DeleteCases[Together /@ (eqs /. rule), 0];
               remaining = DeleteCases[remaining, uu];
               progress = True]],
         {uu, remaining}]];
   <|"solution" -> sol, "remaining" -> eqs, "unsolved" -> remaining|>];
sol0 = solveSeq[eqs0, unkSyms0];
NRH`Check["u^0 order: the log branches h^{(1)}, h^{(1b)}, delta d^{(1)} are all determined and no equation restricts the sources",
   sol0["unsolved"] === {} && sol0["remaining"] === {}];
targ0 = {dd1f[xp, xm] -> -(l/8) D[r0[xp, xm], xp, xm], h1f["mp"][xp, xm] -> 0, h1bf["mp"][xp, xm] -> 0,
   h1f["pp"][xp, xm] -> -(l/2) D[r0[xp, xm], {xp, 2}], h1f["mm"][xp, xm] -> -(l/2) D[r0[xp, xm], {xm, 2}],
   h1bf["pp"][xp, xm] -> 0, h1bf["mm"][xp, xm] -> 0,
   h1bf["pm"][xp, xm] -> (l^2/8) D[r0[xp, xm], {xp, 2}, {xm, 2}],
   h1f["pm"][xp, xm] -> (l/2) (4 D[v0[xp, xm], xp, xm] - D[a0[xp, xm], {xp, 2}] - D[b0[xp, xm], {xm, 2}] - (l^2/4) D[r0[xp, xm], {xp, 2}, {xm, 2}])};
NRH`CheckZero["SMleadingdilaton, SMleadingstress, SMNRsol: the solved u^0 coefficients equal the displayed ones",
   Table[Together[(uu /. sol0["solution"]) - (uu /. targ0)], {uu, unkSyms0}]];

(* order u: full second line of SMcoupledorders with the u^0 solution inserted *)
ansatz1[s_, c_] := (ansatz0[c] /. funRules /. targ0) + z (h2f[c][xp, xm] + yy h2Lf[c][xp, xm] + yy^2 h2LLf[c][xp, xm]);
ans1dd[s_] := (ans0dd /. funRules /. targ0) + z (dd2f[xp, xm] + yy dd2Lf[xp, xm]);
subs1[s_] := {hpp -> Function[{xa, xb, xc}, Evaluate[ansatz1[s, "pp"] /. {xp -> xa, xm -> xb, yy -> xc}]],
   hpm -> Function[{xa, xb, xc}, Evaluate[ansatz1[s, "pm"] /. {xp -> xa, xm -> xb, yy -> xc}]],
   hmp -> Function[{xa, xb, xc}, Evaluate[ansatz1[s, "mp"] /. {xp -> xa, xm -> xb, yy -> xc}]],
   hmm -> Function[{xa, xb, xc}, Evaluate[ansatz1[s, "mm"] /. {xp -> xa, xm -> xb, yy -> xc}]],
   dd -> Function[{xa, xb, xc}, Evaluate[ans1dd[s] /. {xp -> xa, xm -> xb, yy -> xc}]]};
orderUeqs[s_] := Module[{eq, all = {}},
   Do[
      eq = Expand[applyOp[E0target[c], subs1[s]] + z applyOp[E1op[s, c], subs1[s]]];
      eq = Coefficient[eq, z, 1];                       (* the total order-u equation (SMcoupledorders, second line) *)
      all = Join[all, DeleteCases[Together /@ CoefficientList[eq, yy], 0]],
      {c, comps}];
   all];
unkSyms1 = Join[Flatten[Table[{h2LLf[c][xp, xm], h2Lf[c][xp, xm]}, {c, chans}]], {dd2Lf[xp, xm], dd2f[xp, xm]}];
Do[
   Module[{eqs1, sol1, cons, cons2, tcEqs, typeTarget, typeSub, targ1, fp, fm, eqp, eqm, stressRules},
      eqs1 = orderUeqs[s];
      sol1 = solveSeq[eqs1, unkSyms1];
      NRH`Check["u^1 order on " <> s <> ": all h^{(2LL)}, h^{(2L)}, delta d^{(2L)}, delta d^{(2)} are determined algebraically by the log-free responses h^{(2)}",
         sol1["unsolved"] === {}];
      cons = DeleteCases[Together /@ sol1["remaining"], 0];
      NRH`Check["u^1 order on " <> s <> ": the type-changing response h^{(2)}_{op bom} = H_s is unconstrained (SMgeneralradialsolution)", FreeQ[cons, h2f["pm"]]];
      (* SMtypeconstraint: the remaining equations that contain no stress response are d_pm[h^{(2)}_{om bop} - target] and their derivatives *)
      typeTarget = If[s == "R", 4 v0[xp, xm] - (l^2/4) D[r0[xp, xm], xp, xm], 0];
      typeSub = h2f["mp"] -> Function[{xa, xb}, Evaluate[cst + typeTarget /. {xp -> xa, xm -> xb}]];
      tcEqs = Select[cons, FreeQ[#, h2f["pp"]] && FreeQ[#, h2f["mm"]] &];
      NRH`Check["SMtypeconstraint on " <> s <> ": the type-changing channel is constrained (equations in h^{(2)}_{om bop} alone remain, and they are nontrivial)",
         tcEqs =!= {} && ! AllTrue[Together[tcEqs /. h2f["mp"] -> Function[{xa, xb}, gen[xa, xb]]], # === 0 &]];
      NRH`CheckZero["SMtypeconstraint on " <> s <> ": they are solved exactly by h^{(2)}_{om bop} = c_s + (" <> ToString[typeTarget, InputForm] <> ")",
         Together[tcEqs /. typeSub]];
      If[s == "R", NRH`Check["SMtypeconstraint on R: the shift 4 delta d^{(0)} - (l^2/4) d_+ d_- r is required (negative control with h^{(2)}_{om bop} = c_s)",
         ! AllTrue[Together[tcEqs /. h2f["mp"] -> Function[{xa, xb}, cst]], # === 0 &]]];
      (* impose the type constraint; the rest are the two stress constraints SM{R,NR}constraint0/1 and their derivatives *)
      cons2 = DeleteCases[Together /@ (cons /. typeSub), 0];
      eqp = Select[cons2, ! FreeQ[#, Derivative[0, 1][h2f["pp"]][xp, xm]] && FreeQ[#, Derivative[i_, j_][h2f["pp"]][__] /; i + j > 1] && FreeQ[#, h2f["mm"]] &];
      eqm = Select[cons2, ! FreeQ[#, Derivative[1, 0][h2f["mm"]][xp, xm]] && FreeQ[#, Derivative[i_, j_][h2f["mm"]][__] /; i + j > 1] && FreeQ[#, h2f["pp"]] &];
      NRH`Check["u^1 order on " <> s <> ": one equation is linear in d_- h^{(2)}_{op bop} alone and one in d_+ h^{(2)}_{om bom} alone", eqp =!= {} && eqm =!= {}];
      fp = Together[D[h2f["pp"][xp, xm], xm] /. First[Solve[First[eqp] == 0, Derivative[0, 1][h2f["pp"]][xp, xm]]]];
      fm = Together[D[h2f["mm"][xp, xm], xp] /. First[Solve[First[eqm] == 0, Derivative[1, 0][h2f["mm"]][xp, xm]]]];
      NRH`CheckZero["SM" <> s <> "constraint0, SM" <> s <> "constraint1: d_- h^{(2)}_{op bop} = F_+^" <> s <> " and d_+ h^{(2)}_{om bom} = F_-^" <> s <> " with the displayed F_pm",
         {Together[fp - Fptarget[s]], Together[fm - Fmtarget[s]]}];
      stressRules = {Derivative[i_, j_][h2f["pp"]][xp, xm] /; j >= 1 :> D[fp, {xp, i}, {xm, j - 1}],
         Derivative[i_, j_][h2f["mm"]][xp, xm] /; i >= 1 :> D[fm, {xp, i - 1}, {xm, j}]};
      NRH`CheckZero["u^1 order on " <> s <> ": with the type and stress constraints imposed every remaining order-u equation is satisfied (SMgeneralradialsolution: R_pm, H_s and c_s are the free data)",
         Together[cons2 //. stressRules]];
      (* the solved logarithmic and finite coefficients, with the type constraint imposed *)
      NRH`CheckZero["u^1 order on " <> s <> ": all four h^{(2LL)} vanish (no y^2 e^{-2y/l} terms)",
         Together[Table[h2LLf[c][xp, xm] /. sol1["solution"], {c, chans}] /. typeSub]];
      targ1 = Join[Table[h2Lf[c][xp, xm] -> (h2Ltarget[s, c] /. {Rp -> h2f["pp"], Rm -> h2f["mm"]}), {c, chans}],
         {dd2Lf[xp, xm] -> dd2Ltarget[s], dd2f[xp, xm] -> dd2target[s]}];
      NRH`CheckZero["SM" <> s <> "h2L_*, SM" <> s <> "dd2L, SM" <> s <> "dd2: the solved u^1 coefficients equal the displayed ones (R_pm = h^{(2)}_{pm bpm})",
         Table[Together[((uu /. sol1["solution"]) /. typeSub) - (uu /. targ1)], {uu, DeleteCases[unkSyms1, h2LLf[_][__]]}]]],
   {s, {"R", "NR"}}];

(* SMlogfreeconditions: historical sufficient restrictions annihilating the leading logarithmic coefficients; neither an iff nor an all-orders log-free proof. *)
NRH`CheckZero["SMlogfreeconditions: the leading logs vanish if d_+ d_- r = d_+^2 r = d_-^2 r = 0, d_+ d_- delta d^{(0)} = 0 and d_+^2 a + d_-^2 b = 0 (sufficient conditions)",
   Module[{c = {D[r0[xp, xm], xp, xm], D[r0[xp, xm], {xp, 2}], D[r0[xp, xm], {xm, 2}], v0[xp, xm]}},
      Together[({dd1f[xp, xm], h1f["pp"][xp, xm], h1f["mm"][xp, xm], h1bf["pm"][xp, xm], h1f["pm"][xp, xm]} /. targ0)
         /. {Derivative[1, 1][r0][xp, xm] -> 0, Derivative[2, 0][r0][xp, xm] -> 0, Derivative[0, 2][r0][xp, xm] -> 0,
             Derivative[2, 2][r0][xp, xm] -> 0, Derivative[1, 1][v0][xp, xm] -> 0,
             Derivative[2, 0][a0][xp, xm] -> -Derivative[0, 2][b0][xp, xm]}]]];

(* === Direct bridge to the current compact SM2.3 and SM2.4 formulas ===
   All five sources, arbitrary chiral L_pm, arbitrary NR W1; no constant-profile test.
   Inverse derivatives are tested by differentiating their defining equations.
   Chiral zero modes and H_s remain free. The remainder is O(z^2 poly(y)). *)
Module[{r=r0[xp,xm],a=a0[xp,xm],b=b0[xp,xm],v=v0[xp,xm],lp=Lp[xp],lm=Lm[xm],
 app,amm,apm,bpm,dbpp,dbmm,dwp,dwm,corrp,corrm,rlog,nru,solutions,sub,ops,res,stressRules},
 app=b-lp r; amm=a-lm r; apm=2 v-2 f0-l^2/8 D[r,xp,xm];
 bpm=-lm app-lp amm+l^2/4(D[amm,{xp,2}]+D[app,{xm,2}]-2 D[apm,xp,xm]);
 dbpp=D[bpm,xp]+lm D[app,xp]-lp D[amm,xp];
 dbmm=D[bpm,xm]-lm D[app,xm]+lp D[amm,xm];
 dwp=lm/2(l^2 D[r,{xp,3}]-8 lp D[r,xp]-4 D[lp,xp] r);
 dwm=lp/2(l^2 D[r,{xm,3}]-8 lm D[r,xm]-4 D[lm,xm] r);
 corrp=lp (4 v-4 f0-l^2/4 D[r,xp,xm])-l^2/4(3 l^2/4 D[r,{xp,3},xm]-6 lp D[r,xp,xm]-4 D[lp,xp] D[r,xm]);
 corrm=lm (4 v-4 f0-l^2/4 D[r,xp,xm])-l^2/4(3 l^2/4 D[r,xp,{xm,3}]-6 lm D[r,xp,xm]-4 D[lm,xm] D[r,xp]);
 NRH`CheckZero["current SMRintegratedstress: both differentiated b+omega responses solve the constraints",
 {dbpp+dwp+D[corrp,xm]-Fptarget["R"],dbmm+dwm+D[corrm,xp]-Fmtarget["R"]}];
 NRH`CheckZero["current SMNRintegratedstress: differentiating both inverse derivatives gives the stress constraints",
 {l^2/4(lm D[r,{xp,3}]+4 lp D[r,xp,{xm,2}]+3 D[lp,xp] D[r,{xm,2}])+4 lp D[v,xm]
 -2 lp D[a,xp]-D[lp,xp] a-1/4(2 W1[xp,xm] D[r,xp]+D[W1[xp,xm],xp] r)-Fptarget["NR"],
 l^2/4(lp D[r,{xm,3}]+4 lm D[r,{xp,2},xm]+3 D[lm,xm] D[r,{xp,2}])+4 lm D[v,xp]
 -2 lm D[b,xm]-D[lm,xm] b-1/4(2 W1[xp,xm] D[r,xm]+D[W1[xp,xm],xm] r)-Fmtarget["NR"]}];
 rlog=-l^3/4(l^2/8 D[r,{xp,3},{xm,3}]-3/2 lp D[r,xp,{xm,3}]-D[lp,xp] D[r,{xm,3}])
 -l^3/4(l^2/8 D[r,{xm,3},{xp,3}]-3/2 lm D[r,xm,{xp,3}]-D[lm,xm] D[r,{xp,3}])
 -l/2(5 lp lm D[r,xp,xm]+3 D[lm,xm] lp D[r,xp]+3 D[lp,xp] lm D[r,xm]+2 D[lp,xp] D[lm,xm] r);
 nru=l/2(2 lp D[a,xp,xm]+D[lp,xp] D[a,xm]+2 lm D[b,xp,xm]+D[lm,xm] D[b,xp])
 -2 l(lp D[v,{xm,2}]+lm D[v,{xp,2}])+l/2(D[Rp[xp,xm],{xm,2}]+D[Rm[xp,xm],{xp,2}])
 -l^3/4 D[lp D[r,{xm,3}],xp]-l^3/4 D[lm D[r,{xp,3}],xm]-l/8 W1[xp,xm] D[r,xp,xm];
 solutions["R"]={b-l/2 yy D[r,{xp,2}]+z(Rp[xp,xm]+yy(-l/8)(l^2 D[r,{xp,3},xm]-4 lp D[r,xp,xm]-4 D[lp,xp] D[r,xm])),
 c0[xp,xm]+z Hr[xp,xm]+l^2/8 yy^2 D[r,{xp,2},{xm,2}]-l/2 yy(D[b,{xm,2}]+D[a,{xp,2}]-4 D[v,xp,xm]+l^2/4 D[r,{xp,2},{xm,2}])+z yy rlog,
 r+z(4 v-4 f0-l^2/4 D[r,xp,xm]-l/2 yy D[r,xp,xm]),
 a-l/2 yy D[r,{xm,2}]+z(Rm[xp,xm]+yy(-l/8)(l^2 D[r,xp,{xm,3}]-4 lm D[r,xp,xm]-4 D[lm,xm] D[r,xp])),
 v-l/8 yy D[r,xp,xm]+z(l^2/8(D[b,{xm,2}]+D[a,{xp,2}]-4 D[v,xp,xm])-l^2/32(l^2 D[r,{xp,2},{xm,2}]-4 lp D[r,{xm,2}]-4 lm D[r,{xp,2}])-l^3/16 yy D[r,{xp,2},{xm,2}])};
 solutions["NR"]={b+z Rp[xp,xm]-l/2 yy(D[r,{xp,2}]-z D[lp D[r,xm],xp]),
 c0[xp,xm]+z Hn[xp,xm]+l^2/8 yy^2 D[r,{xp,2},{xm,2}]-l/2 yy(D[b,{xm,2}]+D[a,{xp,2}]-4 D[v,xp,xm]+l^2/4 D[r,{xp,2},{xm,2}]-2/l z nru),
 r+z cNR,a+z Rm[xp,xm]-l/2 yy(D[r,{xm,2}]-z D[lm D[r,xp],xm]),
 v+l^2/8 z(lm D[r,{xp,2}]+lp D[r,{xm,2}])-l/8 yy D[r,xp,xm]};
 Do[
 sub=Thread[fieldsH->(Function[{xa,xb,yc},Evaluate[#/.{xp->xa,xm->xb,yy->yc}]]& /@ solutions[branch])];
 ops=Table[If[branch=="R",getE[linR,c],getE[linNR,c]],{c,comps}];
 If[branch=="NR", sub=sub/.NRLpsi];
 stressRules={Derivative[i_,j_][Rp][xp,xm] /; j>=1 :> D[Fptarget[branch],{xp,i},{xm,j-1}],
 Derivative[i_,j_][Rm][xp,xm] /; i>=1 :> D[Fmtarget[branch],{xp,i-1},{xm,j}]};
 res=Map[SeriesZ[#,1]&,applyOp[#,sub]& /@ ops];
 res=res/.stressRules; If[branch=="NR",res=res/.NRLpsi];
 NRH`CheckZero["CURRENT compact "<>branch<>" solution: all ten direct-curvature residuals through z",res];
 Print["  compact ",branch," residuals: ",InputForm[Together /@ res]],
 {branch,{"R","NR"}}];
];

NRH`FileSummary[];
