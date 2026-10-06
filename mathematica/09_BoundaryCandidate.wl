(* 09_BoundaryCandidate.wl | 2026-10-06 standalone edition.
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
(*NRH10 SM6 A Classical Boundary-Theory Candidate*)
(* Scope: the bosonic gauge and extra-transformation checks use non-Abelian covariant derivatives.
   The Grassmann Witt, fermionic and trivial-variation checks use the A=0 realization with partial
   derivatives.  They are not a full non-Abelian fermionic verification or a quantum/holographic test. *)

NRH`BeginFile["09_BoundaryCandidate.wl"];

J4 = ODDJ[2];
H0 = {{0, 0, 1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, -1, 0, 0}};
P0up = J4 . ((J4 + H0)/2) . J4;
Pb0up = J4 . ((J4 - H0)/2) . J4;

DofA = {0, 0, DP1, DP2};
term1 = Sum[P0up[[a, b]] DofA[[a]] phi[b], {a, 4}, {b, 4}];
term2 = Sum[Pb0up[[a, b]] DofA[[a]] phi[b], {a, 4}, {b, 4}];
NRH`Check["the P^(0) term contains only D_+ and the Pbar^(0) term only D_-",
   FreeQ[term1, DP2] && FreeQ[term2, DP1] && ! FreeQ[term1, DP1] && ! FreeQ[term2, DP2]];
NRH`Check["each term selects a single independent component of phi_A",
   Length[Union[Cases[term1, phi[_], Infinity]]] == 1 &&
   Length[Union[Cases[term2, phi[_], Infinity]]] == 1 &&
   Cases[term1, phi[_], Infinity] =!= Cases[term2, phi[_], Infinity]];

g2op = Sqrt[2] {{0, 0}, {1, 0}}; g2om = -Sqrt[2] {{0, 1}, {0, 0}};
NRH`Check["gamma^oplus (gamma^bar-ominus) have rank one: one surviving chiral component each",
   MatrixRank[g2op] == 1 && MatrixRank[g2om] == 1];

mat[f_] := {{f[1][xp, xm], f[2][xp, xm]}, {f[3][xp, xm], f[4][xp, xm]}};
phiP = mat[pp]; phiM = mat[pm]; Ap = mat[aP]; Am = mat[aM]; Lam = mat[lam];
DD[X_, mu_] := D[X, {xp, xm}[[mu]]] - I (({Ap, Am}[[mu]]) . X - X . ({Ap, Am}[[mu]]));
Lbos = Tr[DD[phiP, 1] . DD[phiM, 2]];
deltaGauge[X_] := I (Lam . X - X . Lam);
deltaGaugeA[mu_] := D[Lam, {xp, xm}[[mu]]] - I (({Ap, Am}[[mu]]) . Lam - Lam . ({Ap, Am}[[mu]]));
dLgauge = D[
   Tr[DD2[phiP + t deltaGauge[phiP], 1, Ap + t deltaGaugeA[1], Am + t deltaGaugeA[2]] .
      DD2[phiM + t deltaGauge[phiM], 2, Ap + t deltaGaugeA[1], Am + t deltaGaugeA[2]]] /.
   DD2[X_, mu_, AAp_, AAm_] :> (D[X, {xp, xm}[[mu]]] - I (({AAp, AAm}[[mu]]) . X - X . ({AAp, AAm}[[mu]]))),
   t] /. t -> 0;
NRH`CheckZero["delta_Lambda L = 0 for the non-abelian bosonic sector",
   Together[Expand[dLgauge]]];

sp = {{0, 1}, {0, 0}}; sz = {{1, 0}, {0, -1}}; i2 = IdentityMatrix[2];
th[1] = KroneckerProduct[sp, i2, i2, i2];
th[2] = KroneckerProduct[sz, sp, i2, i2];
th[3] = KroneckerProduct[sz, sz, sp, i2];
th[4] = KroneckerProduct[sz, sz, sz, sp];
NRH`CheckZero["Jordan-Wigner generators anticommute (exact Grassmann algebra)",
   Flatten[Table[th[i] . th[j] + th[j] . th[i], {i, 4}, {j, 4}]]];
id16 = IdentityMatrix[16];

psiP = th[1] pf[xp, xm];
psiM = th[2] mf[xp, xm];
phiPa = id16 bp[xp, xm];
phiMa = id16 bm[xp, xm];
LC = phiHold;
Lcand[phP_, phM_, psP_, psM_] :=
   D[phP, xp] . D[phM, xm] + psP . D[psP, xm] + psM . D[psM, xp];
L0 = Lcand[phiPa, phiMa, psiP, psiM];

vL = {vpf[xp], vmf[xm]};
wD[X_] := vL[[1]] D[X, xp] + vL[[2]] D[X, xm];
dPhiP = wD[phiPa]; dPhiM = wD[phiMa];
dPsiP = wD[psiP] + 1/2 D[vL[[1]], xp] psiP;
dPsiM = wD[psiM] + 1/2 D[vL[[2]], xm] psiM;
dLWitt = D[Lcand[phiPa + t dPhiP, phiMa + t dPhiM, psiP + t dPsiP, psiM + t dPsiM], t] /. t -> 0;
NRH`CheckZero["delta_v L = d_+(v^+ L) + d_-(v^- L) (unit-weight scalar density)",
   Map[Together, Expand[dLWitt - D[vL[[1]] L0, xp] - D[vL[[2]] L0, xm]], {2}]];

xs2 = {xp, xm};
bDblD[e_, m_] := If[m <= 2, 0, D[e, xs2[[m - 2]]]];
GenLie2[xiUp_, HH_] := Module[{dim = 4, xiLow, amat},
   xiLow = J4 . xiUp;
   amat = Table[bDblD[xiUp[[c]], m] - Sum[J4[[c, dd]] bDblD[xiLow[[m]], dd], {dd, 4}], {m, 4}, {c, 4}];
   Sum[xiUp[[c]] bDblD[HH, c], {c, 4}] + amat . HH + HH . Transpose[amat]];
lie2 = GenLie2[{0, 0, vg1[xp, xm], vg2[xp, xm]}, H0];
NRH`Check["Lhat_xi H^(0) = 0 <=> d_- v^+ = 0 = d_+ v^-",
   Module[{eqs = DeleteCases[Union[Flatten[lie2]], 0]},
      Union[Together[eqs /. {Derivative[0, 1][vg1][xp, xm] -> DV1, Derivative[1, 0][vg2][xp, xm] -> DV2}]] ===
      Union[{2 DV1, -2 DV1, 2 DV2, -2 DV2}] ||
      (And @@ (PossibleZeroQ[# /. {Derivative[0, 1][vg1][xp, xm] -> 0, Derivative[1, 0][vg2][xp, xm] -> 0}] & /@ eqs)) &&
      ! FreeQ[eqs, Derivative[0, 1][vg1][xp, xm]] && ! FreeQ[eqs, Derivative[1, 0][vg2][xp, xm]]]];

epsP = th[3] ef[xp];
epsM = th[4] gf[xm];

dphiM1 = psiP . epsP;
dpsiP1 = 1/2 epsP . D[phiPa, xp];
dL1 = D[Lcand[phiPa, phiMa + t dphiM1, psiP + t dpsiP1, psiM], t] /. t -> 0;
NRH`CheckZero["delta_{eps+} L = d_-( psi^+ delta_{eps+} psi^+ )",
   Map[Together, Expand[dL1 - D[psiP . dpsiP1, xm]], {2}]];

dphiP2 = psiM . epsM;
dpsiM2 = 1/2 epsM . D[phiMa, xm];
dL2 = D[Lcand[phiPa + t dphiP2, phiMa, psiP, psiM + t dpsiM2], t] /. t -> 0;
NRH`CheckZero["delta_{eps-} L = d_+( psi^- delta_{eps-} psi^- )",
   Map[Together, Expand[dL2 - D[psiM . dpsiM2, xp]], {2}]];
NRH`Check["the fermionic parameters carry arbitrary chiral profiles (one function per chirality)",
   ! FreeQ[dL1, ef] && ! FreeQ[dL2, gf]];

dphiPz = zb[xm] DD[phiM, 2];
dphiMz = zf[xp] DD[phiP, 1];
dLz = D[Tr[DD2[phiP + t dphiPz, 1] . DD2[phiM + t dphiMz, 2]] /.
   DD2[X_, mu_] :> (D[X, {xp, xm}[[mu]]] - I (({Ap, Am}[[mu]]) . X - X . ({Ap, Am}[[mu]]))), t] /. t -> 0;
NRH`CheckZero["delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]",
   Together[Expand[dLz
      - 1/2 D[zb[xm] Tr[DD[phiM, 2] . DD[phiM, 2]], xp]
      - 1/2 D[zf[xp] Tr[DD[phiP, 1] . DD[phiP, 1]], xm]]]];

psiP2 = th[1] pf1[xp, xm] + th[3] pf3[xp, xm];
psiM2 = th[2] mf2[xp, xm] + th[4] mf4[xp, xm];
a0 = af0[xp, xm]; aPl = afp[xp, xm]; aMi = afm[xp, xm];
dpsiPa = a0 D[psiM2, xp] + aPl D[psiP2, xm];
dpsiMa = a0 D[psiP2, xm] + aMi D[psiM2, xp];
dLa = D[Lcand[phiPa, phiMa, psiP2 + t dpsiPa, psiM2 + t dpsiMa], t] /. t -> 0;
NRH`CheckZero["delta_alpha L = d_-(psi^+ delta_alpha psi^+) + d_+(psi^- delta_alpha psi^-) identically (off shell)",
   Map[Together, Expand[dLa - D[psiP2 . dpsiPa, xm] - D[psiM2 . dpsiMa, xp]], {2}]];
NRH`Check["the transformations are built from the fermion equations of motion d_- psi^+ = 0 = d_+ psi^-, hence vanish on shell",
   Module[{onshell = {Derivative[0, 1][pf1][xp, xm] -> 0, Derivative[0, 1][pf3][xp, xm] -> 0,
       Derivative[1, 0][mf2][xp, xm] -> 0, Derivative[1, 0][mf4][xp, xm] -> 0}},
      Union[Flatten[{dpsiPa, dpsiMa} /. onshell]] === {0}]];
NRH`Check["the two-component fermions used here are genuinely Grassmann: psi^+ d_- psi^+ is nonzero off shell",
   ! (Union[Flatten[Expand[psiP2 . D[psiP2, xm]]]] === {0})];

NRH`FileSummary[];
