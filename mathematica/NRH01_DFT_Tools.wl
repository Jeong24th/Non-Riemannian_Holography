(* ::Package:: *)

(* ::Title:: *)
(*NRH01 — Shared DFT Toolbox*)


(* ::Text:: *)
(*This file is the shared toolbox for the Mathematica verification of*)
(*"Long Strings and Non-Riemannian Hair" (NR_Holography.tex, Letter + Supplemental Material).*)
(**)
(*It implements, in plain matrix language, the semi-covariant Double Field Theory (DFT)*)
(*machinery in the conventions of the paper and of the DFT core lecture notes*)
(*[arXiv reference: Park, "Introduction to Double Field Theory / core conventions"]:*)
(**)
(*  -  the O(D,D) invariant metric J_{MN} in the ordered doubled basis x^M = (dual x~_mu ; physical x^mu),*)
(*  -  the projectors P = (J+H)/2 and Pbar = (J-H)/2,*)
(*  -  the DFT Christoffel connection Gamma_{CAB}   [lecture notes Eq. (2.41); SM display of the paper],*)
(*  -  its field strength R_{CDAB}                  [Eq. (2.90)],*)
(*  -  the semi-covariant Riemann curvature S_{ABCD} [Eq. (2.89)],*)
(*  -  the Ricci curvature S_{AB} = S^{C}{}_{ACB}    [Eq. (2.101)],*)
(*  -  the scalar curvature S_(0), computed BOTH from the closed form [Eq. (2.103)]*)
(*     and from the S_{ABCD} contraction [Eq. (2.102)] as an internal cross-check,*)
(*  -  the Einstein curvature G_{AB} = 4(P S Pbar)_{[AB]} - (1/2) J_{AB} S_(0)  [Eq. (2.115)],*)
(*  -  the generalized Lie derivative of H_{MN} and d,*)
(*  -  the Gamma^2 Lagrangian density               [Eq. (2.122)],*)
(*  -  the double-vielbein spin connection Phi_{Apq} (as in the published Python guards).*)
(**)
(*The section files NRH02-NRH08 Get[] this file and verify the displayed equations of the*)
(*paper in the order in which they appear.  Every check is registered through the small*)
(*NRH` framework below, so that each file ends with an explicit PASS/FAIL summary and*)
(*NRH00_RunAll.wl can print a grand total.*)
(**)
(*The connection routine follows the structure of the SymPy implementation in*)
(*calculations/verify_10d_killing_spinor.py of the public verification archive, with the*)
(*trace vector of the third structure fixed by the dilaton compatibility of the connection*)
(*(see the implementation note above GammaDFT).  Antisymmetrization brackets have unit*)
(*weight, e.g. X_[AB] = (X_AB - X_BA)/2.  Every property that defines the connection -*)
(*nabla P = 0, the dilaton trace, vanishing totally antisymmetric part - and the curvature*)
(*identities (algebraic Bianchi, pair symmetry, projective property, the commutator*)
(*identity that defines R_{CDAB}) are themselves verified in NRH02 on exact backgrounds.*)


(* ::Section:: *)
(*Bookkeeping framework (context NRH`, survives ClearAll["Global`*"])*)


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

(* A check that an expression (scalar, list, or matrix) is identically zero. *)
NRH`CheckZero[label_String, expr_] := Module[{z},
   z = NRH`ZeroQ[expr];
   NRH`Record[label, TrueQ[z]]];

(* A check of a boolean statement. *)
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

(* Robust zero test: cheap normal forms first, Simplify only as a fallback. *)
NRH`ZeroQ[expr_] := Module[{flat, t},
   flat = Flatten[{expr}];
   AllTrue[flat,
      Function[e,
         t = Together[Expand[e]];
         If[t === 0, True,
            t = Together[ExpandAll[TrigToExp[t]]];
            If[t === 0, True, PossibleZeroQ[Simplify[t]]]]]]];


(* ::Section:: *)
(*O(D,D) metric, projectors, and the section rule*)


(* ::Text:: *)
(*The ordered doubled basis is x^M = (x~_mu ; x^mu):  the first Dphys slots are the dual*)
(*(winding) coordinates, the last Dphys slots the physical coordinates.  In this basis the*)
(*O(D,D) metric is off-block-diagonal, and the section condition  \[PartialD]~^mu = 0  says that all*)
(*fields depend only on the physical coordinates: the doubled derivative \[PartialD]_M acts as 0 on*)
(*the first block and as \[PartialD]/\[PartialD]x^mu on the second.*)


ODDJ[nphys_Integer] := ArrayFlatten[{{0, IdentityMatrix[nphys]}, {IdentityMatrix[nphys], 0}}];

(* Doubled derivative of a scalar/matrix expression with respect to doubled index m.
   An entry of xs may be either a plain coordinate symbol (ordinary D), or a pure
   function implementing a weighted/compound derivative.  The latter is used to work
   with rational radial variables, e.g. u = Exp[2y/l] with  d/dy = (2u/l) d/du,
   which keeps every intermediate expression rational and fast. *)
DblD[expr_, m_Integer, xs_List] := Module[{n = Length[xs], op},
   If[m <= n, 0*expr,
      op = xs[[m - n]];
      If[Head[op] === Function, op[expr], D[expr, op]]]];

(* Gradient as a doubled (lower-index) vector. *)
DblGrad[expr_, xs_List] := Table[DblD[expr, m, xs], {m, 1, 2 Length[xs]}];


(* ::Section:: *)
(*DFT Christoffel connection Gamma_{CAB}  (torsionless, Eq. (2.41))*)


(* ::Text:: *)
(*GammaDFT returns the list {Gamma[[C]]}_{C=1..2D} of matrices in the last two (antisymmetric)*)
(*indices, Gamma_{CAB}.  The three structures are*)
(*   2 (P dP_C P Pbar)_[AB]*)
(* + 2 (Pbar_[A^D Pbar_B]^E - P_[A^D P_B]^E) dP_{D,EC}*)
(* - 4/(Dphys-1) (Pbar_C[A Pbar_B]^D + P_C[A P_B]^D) (dd_D + (P d^E P Pbar)_[ED]).*)
(*This is a direct port of semi_covariant_connection() from the tested Python guard.*)


(* ::Text:: *)
(*Implementation note.  The first two structures are assembled literally.  The trace vector*)
(*X_D multiplying the third structure is fixed, unambiguously, by the defining dilaton*)
(*compatibility of the torsionless semi-covariant connection,*)
(*    Gamma^B{}_{BA} = -2 dd_A     (equivalently  nabla_A e^{-2d} = 0),*)
(*because the J-trace of the third structure equals  -2 X_A  identically.  Hence*)
(*    X_A = dd_A + (1/2) J^{BC} (Gamma_structures12)_{CBA} ,*)
(*which is the meaning of the symbolic combination  dd_D + (P dP^E P Pbar)_[ED]  in the*)
(*printed formula.  With this X the connection satisfies ALL of its defining properties*)
(*(nabla P = 0, nabla J = 0, dilaton compatibility, vanishing totally antisymmetric part),*)
(*and the curvature checks of the section files then hold exactly - including for chiral,*)
(*x-dependent backgrounds.  A naive row-projected reading of (P dP^E P Pbar)_[ED] satisfies*)
(*nabla P = 0 but violates the dilaton trace on x-dependent backgrounds and is NOT used.*)


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
         (* second structure, summed over the derivative index in the physical block *)
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


(* ::Section:: *)
(*Curvatures:  R_{CDAB},  Ricci S_{AB},  scalar S_(0),  Einstein G_{AB}*)


(* ::Text:: *)
(*R_{CDAB} = dGamma_{BCD,A} - dGamma_{ACD,B} + Gamma_{AC}{}^E Gamma_{BED} - Gamma_{BC}{}^E Gamma_{AED}.*)
(*We store it as R4[[a,b]] = matrix in (C,D).  The Ricci curvature is the contraction*)
(*S_{AB} = S^{C}{}_{ACB} of the semi-covariant Riemann tensor*)
(*S_{ABCD} = (R_{ABCD} + R_{CDAB} - Gamma^E{}_{AB} Gamma_{ECD})/2.*)
(*Raising one doubled index with J in the ordered basis simply swaps the two blocks:*)
(*partner[m] = m + D  (m <= D)  or  m - D  (m > D).*)


RiemannR4[gamma_List, xs_List] := Module[{n = Length[xs], dim, JJ},
   dim = 2 n; JJ = ODDJ[n];
   Table[
      If[b <= a, ConstantArray[0, {dim, dim}],   (* antisymmetric in the derivative pair; fill below *)
         Map[Together,
            DblD[gamma[[b]], a, xs] - DblD[gamma[[a]], b, xs]
            + gamma[[a]] . JJ . gamma[[b]] - gamma[[b]] . JJ . gamma[[a]], {2}]],
      {a, 1, dim}, {b, 1, dim}]
   // (# - Transpose[#, {2, 1, 3, 4}] &)];

Partner[m_Integer, n_Integer] := If[m <= n, m + n, m - n];

(* Ricci S_{AB}; r4 = RiemannR4 output. *)
RicciS[gamma_List, r4_List, xs_List] := Module[{n = Length[xs], dim, gg},
   dim = 2 n;
   (* GammaGamma_{EACB} = sum_F Gamma^F_{EA} Gamma_{FCB} = sum_F Gamma[[partner[F]]]_{EA} Gamma[[F]]_{CB} *)
   Table[
      Together[Sum[Module[{e = Partner[c, n]},
         (r4[[c, b, e, a]] + r4[[e, a, c, b]]
            - Sum[gamma[[Partner[f, n], e, a]]*gamma[[f, c, b]], {f, 1, dim}])/2],
         {c, 1, dim}]],
      {a, 1, dim}, {b, 1, dim}]];

(* Scalar curvature from the closed form (2.103): uses only H_{MN} and d. *)
ScalarS0[HH_, dd_, xs_List] := Module[
   {n = Length[xs], dim, JJ, Hup, Hmix, gradd, term},
   dim = 2 n; JJ = ODDJ[n];
   Hup = JJ . HH . JJ;    (* H^{AB} *)
   Hmix = HH . JJ;        (* H_A{}^B *)
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

(* Scalar curvature from the S_{ABCD} contraction (2.102): independent cross-check.  *)
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

(* Projected Ricci (P S Pbar)_{AB} and Einstein curvature G_{AB}. *)
ProjectedRicci[HH_, ricci_, xs_List] := Module[{n = Length[xs], JJ, P, Pb},
   JJ = ODDJ[n]; P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   P . JJ . ricci . JJ . Pb];

EinsteinG[HH_, ricci_, s0_, xs_List] := Module[{n = Length[xs], JJ, psp},
   JJ = ODDJ[n];
   psp = ProjectedRicci[HH, ricci, xs];
   4*(psp - Transpose[psp])/2 - 1/2*JJ*s0];


(* ::Section:: *)
(*Generalized Lie derivative*)


(* ::Text:: *)
(*For a doubled vector with upper components xi^M = (lambda_mu ; xi^mu) in the ordered basis:*)
(*  (Lhat_xi H)_{MN} = xi^C dH_{MN,C} + (dxi^C_{,M} - d^C xi_{M,}) H_{CN} + (M <-> N on the second slot),*)
(*  Lhat_xi d = xi^A dd_A - (1/2) dxi^A_A.*)


GenLieH[xiUp_List, HH_, xs_List] := Module[
   {n = Length[xs], dim, JJ, xiLow, dxiUp, dxiLow, amat},
   dim = 2 n; JJ = ODDJ[n];
   xiLow = JJ . xiUp;
   dxiUp = Table[DblD[xiUp[[c]], m, xs], {m, 1, dim}, {c, 1, dim}];   (* dxiUp[[M,C]] = d_M xi^C *)
   dxiLow = Table[DblD[xiLow[[c]], m, xs], {m, 1, dim}, {c, 1, dim}];
   (* amat[[M,C]] = d_M xi^C - d^C xi_M ;  d^C = J^{CD} d_D *)
   amat = Table[dxiUp[[m, c]] - Sum[JJ[[c, dd2]]*dxiLow[[dd2, m]], {dd2, 1, dim}], {m, 1, dim}, {c, 1, dim}];
   Sum[xiUp[[c]]*DblD[HH, c, xs], {c, 1, dim}] + amat . HH + HH . Transpose[amat]];

GenLieD[xiUp_List, dd_, xs_List] := Module[{n = Length[xs], dim},
   dim = 2 n;
   Sum[xiUp[[a]]*DblD[dd, a, xs], {a, 1, dim}] - 1/2*Sum[DblD[xiUp[[a]], a, xs], {a, 1, dim}]];


(* ::Section:: *)
(*Riemannian packaging:  H_{MN} and e^{-2d} from (g, B, phi)*)


(* ::Text:: *)
(*In the Riemannian parametrization the generalized metric and dilaton are*)
(*   H = ( g^{-1} , -g^{-1} B ;  B g^{-1} , g - B g^{-1} B ),    e^{-2d} = Sqrt[-det g] e^{-2 phi}.*)


RiemannianH[g_, B_] := Module[{gi = Inverse[g]},
   ArrayFlatten[{{gi, -gi . B}, {B . gi, g - B . gi . B}}]];

RiemannianDilaton[g_, phi_] := phi - 1/4 Log[-Det[g]];  (* d with e^{-2d} = Sqrt[-g] e^{-2 phi} *)


(* ::Section:: *)
(*Gamma^2 Lagrangian density (2.122) and its boundary vector B^A (2.121)*)


Gamma2Density[HH_, dd_, gamma_List, xs_List] := Module[
   {n = Length[xs], dim, JJ, P, Pb, Pup, Pbup, gup},
   dim = 2 n; JJ = ODDJ[n];
   P = (JJ + HH)/2; Pb = (JJ - HH)/2;
   Pup = JJ . P . JJ; Pbup = JJ . Pb . JJ;
   (* Gamma_{AC}{}^{E} = gamma[[a,c,f]] J^{fE} : contract with partner. *)
   Exp[-2 dd]*Sum[(Pup[[a, c]]*Pup[[b, d]] - Pbup[[a, c]]*Pbup[[b, d]])*
        Sum[gamma[[a, c, Partner[e, n]]]*gamma[[b, d, e]]
            - gamma[[a, b, Partner[e, n]]]*gamma[[d, c, e]]
            + 1/2*gamma[[Partner[e, n], a, b]]*gamma[[e, c, d]], {e, 1, dim}],
      {a, 1, dim}, {b, 1, dim}, {c, 1, dim}, {d, 1, dim}]];

GammaBVector[HH_, dd_, xs_List] := Module[{n = Length[xs], dim, JJ, Hup},
   dim = 2 n; JJ = ODDJ[n]; Hup = JJ . HH . JJ;
   Table[4*Sum[Hup[[a, b]]*DblD[dd, b, xs], {b, 1, dim}]
         - Sum[DblD[Hup[[a, b]], b, xs], {b, 1, dim}], {a, 1, dim}]];


(* ::Section:: *)
(*Double-vielbein spin connection Phi_{Apq}  (port of spin_connection() from the guard)*)


(* ::Text:: *)
(*V is a (2D x k) double vielbein with flat metric eta (k x k);  gamma is the GammaDFT list.*)
(*Phi_{Apq} = V^B{}_p ( dV_{Bq,A} + Gamma_{AB}{}^{C} V_{Cq} ),  antisymmetrized in (p,q),*)
(*returned as a list over A of (k x k) matrices with both flat indices lowered.*)


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


(* ::Section:: *)
(*Convenience: full curvature pipeline*)


(* ::Text:: *)
(*DFTCurvature[H, d, xs] returns an Association with the connection, Ricci, scalar (closed form),*)
(*projected Ricci, and Einstein curvature.  For the 3-dimensional saddles of the paper this*)
(*runs exactly; for the 10-dimensional uplift the section files use the same functions*)
(*with numerical high-precision spot checks where full symbolics would be too slow.*)


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
