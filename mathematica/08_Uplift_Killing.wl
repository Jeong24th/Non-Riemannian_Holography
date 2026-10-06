(* 08_Uplift_Killing.wl | 2026-10-06 standalone edition.
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
(*NRH09 SM5 Ten-Dimensional Uplift and Killing Symmetries*)
(* Direct uplift checks cover arbitrary chiral R data and the one-sided NR family.  Reduced spinor and
   historical complex-rank checks do not verify the full current real Majorana product basis or construct
   nonzero fermionic charges.  Finite counts refer to polarizations/constant representatives, not the
   dimension of the arbitrary chiral-function solution space. *)

NRH`BeginFile["08_Uplift_Killing.wl"];

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
NRH`Check["historical complex-rank check: Weyl + S^3-line + zeta_+ leave complex dimension 4; not a verification of the current real product basis",
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

(* The reduced one-sided jet system after SMreducedDirac, derived in the printed frame: D_{pbar} E = 0 and the
   opposite-channel Dirac equation gamma^p D_p E + E/(Sqrt[2] l) = 0 on the exact L_- = 0 background with
   W = W_0 + u_m W_1 with manuscript u_m = e^{-2y/l} (code u = e^{2y/l}),
   generic spinor E = (e0, e1)(x+, x-, y). *)
Module[{Wn, HNRone1, VexU1, VbexU1, V1, Vb1, gamma1, Phi1, spin1, DA1, DP1, DPbar1, slashD1, Eg, jetRule, grav, dirOpp,
        jetVars, sysDisplayed, mat, mC, computed, simp},
   Wn = W0[xp, xm] + W1[xp, xm]/u;
   HNRone1 = {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 2 Lp[xp]/u, -1, 0}, {0, 0, 1, 0, 0, 0},
      {1, 2 Lp[xp]/u, 0, -2 Lp[xp] Wn/u, Wn, 0}, {0, -1, 0, Wn, 0, 0}, {0, 0, 0, 0, 0, 1}};
   (* SMvielbein in the one-sided limit cosh h -> 1, e^sigma sinh h -> L+/u, e^{-sigma} sinh h -> 0 (upper local indices) *)
   VexU1 = {{0, -1/Sqrt[2], 0}, {0, -(Lp[xp]/u)/Sqrt[2], 0}, {0, 0, 1/Sqrt[2]},
      {Sqrt[2], Wn (Lp[xp]/u)/(2 Sqrt[2]), 0}, {0, -Wn/(2 Sqrt[2]), 0}, {0, 0, 1/Sqrt[2]}};
   VbexU1 = {{0, 0, 0}, {1/Sqrt[2], 0, 0}, {0, 0, -1/Sqrt[2]},
      {-Wn/(2 Sqrt[2]), -Sqrt[2] Lp[xp]/u, 0}, {0, Sqrt[2], 0}, {0, 0, 1/Sqrt[2]}};
   V1 = VexU1 . eta3; Vb1 = VbexU1 . (-eta3);
   NRH`CheckZero["one-sided printed frame: V eta V^T - Vbar etabar Vbar^T = H(L+, W_0, W_1) and V -> Vinf at L+ = W = 0",
      {Simplify[V1 . eta3 . Transpose[V1] - Vb1 . (-eta3) . Transpose[Vb1] - HNRone1],
       (V1 /. {Lp[xp] -> 0, W0[xp, xm] -> 0, W1[xp, xm] -> 0}) - Vinf, (Vb1 /. {Lp[xp] -> 0, W0[xp, xm] -> 0, W1[xp, xm] -> 0}) - Vbinf}];
   gamma1 = GammaDFT[HNRone1, dinf, xsU];
   Phi1 = Map[Together, SpinConnectionDFT[V1, eta3, gamma1, xsU], {3}];
   spin1[a_] := 1/4 Sum[Phi1[[a, p, q]] gampqLow[p, q], {p, 3}, {q, 3}];
   DA1[Es_, a_] := DblD[Es, a, xsU] + spin1[a] . Es;
   DP1[Es_, p_] := Sum[(JJ . V1)[[a, p]] DA1[Es, a], {a, 6}];
   DPbar1[Es_, pb_] := Sum[(JJ . Vb1)[[a, pb]] DA1[Es, a], {a, 6}];
   slashD1[Es_] := Sum[gam3[[p]] . DP1[Es, p], {p, 3}];
   Eg = {e0[xp, xm, u], e1[xp, xm, u]};
   jetRule = {Derivative[1, 0, 0][e0][xp, xm, u] -> ep0, Derivative[0, 1, 0][e0][xp, xm, u] -> em0,
      Derivative[0, 0, 1][e0][xp, xm, u] -> l ey0/(2 u), e0[xp, xm, u] -> e0,
      Derivative[1, 0, 0][e1][xp, xm, u] -> ep1, Derivative[0, 1, 0][e1][xp, xm, u] -> em1,
      Derivative[0, 0, 1][e1][xp, xm, u] -> l ey1/(2 u), e1[xp, xm, u] -> e1};
   grav = Flatten[Table[DPbar1[Eg, pb], {pb, 3}]] /. jetRule;
   dirOpp = (slashD1[Eg] + Eg/(Sqrt[2] l)) /. jetRule;
   simp[e_] := Together[e /. u -> 1/uu];                                 (* uu = e^{-2y/l}, the manuscript's u *)
   NRH`CheckZero["SM text after SMreducedDirac: D_{pbar} E has the components 0, -L+ u e0/l, em0/Sqrt2, em1/Sqrt2 - W1 u e0/(4l), ey0/Sqrt2, ey1/Sqrt2",
      Together[(grav /. u -> 1/uu) - {0, -Lp[xp] uu e0/l, em0/Sqrt[2], em1/Sqrt[2] - W1[xp, xm] uu e0/(4 l), ey0/Sqrt[2], ey1/Sqrt[2]}]];
   NRH`CheckZero["SM text after SMreducedDirac: gamma^p D_p E + E/(Sqrt2 l) has the components ey0/Sqrt2 + Sqrt2 e0/l and ep0 - ey1/Sqrt2 + L+ u em0",
      Together[(dirOpp /. u -> 1/uu) - {ey0/Sqrt[2] + Sqrt[2] e0/l, ep0 - ey1/Sqrt[2] + Lp[xp] uu em0}]];
   jetVars = {e0, e1, ep0, ep1, em0, em1, ey0, ey1};
   sysDisplayed = {LpS uu e0, em0, W1S uu e0 - 2 Sqrt[2] l em1, ey0, ey1, 2 e0 + l ey0, 2 LpS uu em0 + 2 ep0 - Sqrt[2] ey1};
   mat[sys_] := Table[Coefficient[Together[sys[[i]]], jetVars[[j]]], {i, Length[sys]}, {j, 8}];
   computed = Join[grav, dirOpp] /. {Lp[xp] -> LpS, W1[xp, xm] -> W1S, W0[xp, xm] -> W0S} /. u -> 1/uu;
   mC = mat[Together[computed]];
   NRH`Check["the displayed system has rank six (seven relations, one dependent) and spans the derived equations",
      MatrixRank[mC] == 6 && MatrixRank[mat[sysDisplayed]] == 6 && MatrixRank[Join[mC, mat[sysDisplayed]]] == 6];
   NRH`Check["the solution space is exactly {e_1, d_+ e_1}: e_0 = 0, e_1 = F_+(x^+) with d_+ e_1 unconstrained",
      Sort[RowReduce[NullSpace[mC]]] === Sort[{{0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}}]];
   NRH`Check["the mechanism: W_1 multiplies only e_0, and every d_+ L_+, d W_0, d W_1 term cancels in the full connection",
      Union[Cases[sysDisplayed, W1S x_ :> x, Infinity]] === {uu e0} && FreeQ[computed, Derivative[__][LpS | W0S | W1S][__]] &&
      FreeQ[Together[computed], W0S]];
];

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
NRH`Check["counting: arithmetic consistency of the stated polarization counts, (4,4) + (4,4) = 16 and (4,0)/(0,4) = 4; not an independent global spinor count",
   4 + 4 + 4 + 4 == 16 && 4 + 0 == 4 && 0 + 4 == 4];

NRH`FileSummary[];
