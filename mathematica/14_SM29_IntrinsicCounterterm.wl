(* 14_SM29_IntrinsicCounterterm.wl | 2026-10-08 standalone edition.
   Reproduce SMboundarycurvaturecandidate (SM29), not an assumed general-state
   nonlinear completion. No Get, Needs, packages, or external input files.
   Run in a fresh Wolfram kernel; all Global` definitions are cleared.

   CONVENTIONS: intrinsic D=2, doubled order (tilde x+,tilde x-;x+,x-),
   section partial_tilde=0, eta={ {0,-1},{-1,0} }, etabar=-eta.
   All densities below omit the common action factor -1/(32 Pi G).
   eps is source amplitude; u=Exp[-2 Y/l]. First expand in eps to degree 2,
   then substitute the full linear source-to-cutoff map, then retain u^-1
   and u^0 (including every power of Y). Never take Y->infinity first.
   Equalities of integrated densities are modulo tangential total derivatives,
   with compact support or periodic boundary conditions. Lp,Lm,W1 are constants.
   Arbitrary boundary functions are left symbolic, with exact coefficients.

   Curvature scalar: core lecture (2.103). Mixed Ricci: intrinsic D=2
   DFT Christoffel connection, independently checked by the constant-background
   variation formula and the Euler derivative of the scalar-curvature action.
   The affine dilaton and constrained quadratic coset completion are essential.
*)
ClearAll["Global`*", "NRH`*"];
NRH`$AllResults = {}; NRH`$FileResults = {};
NRH`$CurrentFile = "14_SM29_IntrinsicCounterterm.wl";
NRH`Record[label_String, ok_] := Module[{b = TrueQ[ok]},
 AppendTo[NRH`$FileResults, {NRH`$CurrentFile, label, b}];
 AppendTo[NRH`$AllResults, {NRH`$CurrentFile, label, b}];
 Print[If[b, "  [PASS] ", "  [FAIL] "], label]; b];
NRH`ZeroQ[e_] := AllTrue[Flatten[{e}], Together[Expand[#]] === 0 &];
NRH`CheckZero[label_String, e_] := NRH`Record[label, NRH`ZeroQ[e]];
NRH`Check[label_String, e_] := NRH`Record[label, e];
NRH`FileSummary[] := Module[{bad = Select[NRH`$FileResults, !TrueQ[#[[3]]] &]},
 Print["----------------------------------------------------------------"];
 Print["  ", NRH`$CurrentFile, ": ", Length[NRH`$FileResults] - Length[bad],
   "/", Length[NRH`$FileResults], " checks passed."];
 If[bad =!= {}, Print["  FAILED: ", bad[[All, 2]]];
   If[$FrontEnd === Null, Exit[1], Abort[]]]; True];
Print["SMboundarycurvaturecandidate: intrinsic D=2 quadratic derivation"];
Print["u=Exp[-2 Y/l]; all printed densities omit -1/(32 Pi G)."];

(* Exact polynomial operations at the orders needed by the nondecaying action. *)
Trunc[e_, z_, n_Integer] := Total[Table[Coefficient[Expand[e], z, k] z^k, {k, 0, n}]] // Expand;
TU[e_] := Trunc[e, u, 1];
TE[e_] := Trunc[e, eps, 2];
TUmat[e_] := Map[TU, e, {2}];
TEmat[e_] := Map[TE, e, {2}];
Nondecay[e_] := Expand[Coefficient[Expand[u e], u, 0]/u + Coefficient[Expand[u e], u, 1]];
(* Explicit version avoids treating logarithms Y as powers of u. *)
DivPart[e_] := Expand[Coefficient[Expand[u e], u, 0]/u +
 Coefficient[Expand[u e], u, 1] - (Coefficient[Expand[u e], u, 1] /. Y -> 0)];
Der[e_, m_Integer : 0, n_Integer : 0] := D[e, {xp, m}, {xm, n}];
DB[e_, i_Integer] := Switch[i, 1 | 2, 0 e, 3, D[e, xp], 4, D[e, xm]];
jmat = ArrayFlatten[{{0, IdentityMatrix[2]}, {IdentityMatrix[2], 0}}];
eta = {{0, -1}, {-1, 0}}; etab = -eta;
fields = {a0, b0, r0, c0, v0, hs, rp, rm};
freeFields = {fPP, fMM, fMP, fPM, fd};

(* Bilinear integration by parts: transfer all derivatives to the field later
   in heads. Unlike replacing a density by zero, this retains its functional. *)
JetInfo[z_, heads_List] := Module[{t, p},
 t = Which[
   MatchQ[z, Derivative[_Integer, _Integer][_Symbol][xp, xm]],
     z /. Derivative[m_, n_][f_][xp, xm] :> {f, m, n},
   MatchQ[z, _Symbol[xp, xm]], {Head[z], 0, 0}, True, Return[{}]];
 p = FirstPosition[heads, t[[1]]];
 If[MissingQ[p], {}, {First[p], t[[2]], t[[3]]}]];
IBPTerm[t_, heads_List] := Module[{fac, jets = {}, cf = 1, info, base, pow, q, i, m, n},
 fac = If[Head[t] === Times, List @@ t, {t}];
 Do[base = If[Head[z] === Power, z[[1]], z];
   pow = If[Head[z] === Power, z[[2]], 1]; info = JetInfo[base, heads];
   If[info === {}, cf *= z,
     If[!IntegerQ[pow] || pow < 0 || pow > 2, Print["Invalid bilinear jet: ", z]; Abort[]];
     jets = Join[jets, ConstantArray[info, pow]]], {z, fac}];
 q = Length[jets];
 Which[q == 0, cf,
 q == 1, {i, m, n} = First[jets]; If[m + n == 0, cf heads[[i]][xp, xm], 0],
 q == 2, jets = Sort[jets]; m = Total[jets[[All, 2]]]; n = Total[jets[[All, 3]]];
   If[jets[[1, 1]] == jets[[2, 1]] && OddQ[m + n], 0,
     (-1)^Total[jets[[1, {2, 3}]]] cf heads[[jets[[1, 1]]]][xp, xm]
       Der[heads[[jets[[2, 1]]]][xp, xm], m, n]],
 True, Print["Expected a polynomial at most quadratic in sources: ", t]; Abort[]]];
IBP[e_, heads_List] := Module[{ex = Expand[e]},
 Expand[Total[IBPTerm[#, heads] & /@ If[Head[ex] === Plus, List @@ ex, {ex}]]]];
VariationalD[lag_, field_] := Module[{fn = Head[field], jets, inds, vars, poly},
 jets = DeleteDuplicates[Join[{field}, Cases[lag,
   z : Derivative[__][f_][__] /; f === fn :> z, Infinity]]];
 inds = (If[# === field, {0, 0}, List @@ Head[Head[#]]] & /@ jets);
 vars = Table[Unique["jet"], {Length[jets]}]; poly = lag /. Thread[jets -> vars];
 Expand[Sum[(-1)^Total[inds[[k]]] Der[D[poly, vars[[k]]] /. Thread[vars -> jets],
   inds[[k, 1]], inds[[k, 2]]], {k, Length[jets]}]]];

(* Core lecture (2.103), used at finite cutoff in D=2, not bulk D=3. *)
ScalarDFT[hh_, dd_] := Module[{hu = jmat . hh . jmat, hm = hh . jmat, gd},
 gd = Table[DB[dd, i], {i, 4}];
 Sum[hu[[i, j]] (1/8 Sum[DB[hu[[k, q]], i] DB[hh[[k, q]], j], {k, 4}, {q, 4}]
   + 1/2 Sum[DB[hm[[i, q]], k] DB[hm[[j, k]], q], {k, 4}, {q, 4}]
   - 4 gd[[i]] gd[[j]] + 4 DB[gd[[j]], i]), {i, 4}, {j, 4}]
 - Sum[DB[DB[hu[[i, j]], i], j], {i, 4}, {j, 4}]
 + 4 Sum[DB[hu[[i, j]], i] gd[[j]], {i, 4}, {j, 4}]];

(* Linear D=2 DFT Christoffel symbols about an intrinsic constant background.
   Lowered Gamma_CAB; unit-weight antisymmetrization. Trace is nabla_A d=0.
   The factor -4/(D-1) is therefore -4, not the bulk D=3 value -2. *)
GammaLinear[hbg_, dh_, dd_] := Module[
 {p = (jmat + hbg)/2, pb = (jmat - hbg)/2, pm, pbm, dp,
  gam12, tr, xx, pmx, pbmx, t1, t2, cqa, cqba},
 pm = p . jmat; pbm = pb . jmat; dp = Table[DB[dh/2, i], {i, 4}];
 gam12 = Table[t1 = pm . dp[[c]] . jmat . pb; t1 = t1 - Transpose[t1];
   t2 = Sum[cqba = (pbm . dp[[q]])[[All, c]]; cqa = (pm . dp[[q]])[[All, c]];
     Outer[Times, pbm[[All, q]], cqba] - Outer[Times, cqba, pbm[[All, q]]]
     - Outer[Times, pm[[All, q]], cqa] + Outer[Times, cqa, pm[[All, q]]], {q, 3, 4}];
   TUmat[t1 + t2], {c, 4}];
 tr = Table[TU[Sum[jmat[[b, q]] gam12[[q, b, a]], {b, 4}, {q, 4}]], {a, 4}];
 xx = Table[DB[dd, a], {a, 4}] + tr/2; pmx = pm . xx; pbmx = pbm . xx;
 Table[TUmat[gam12[[c]] - 2 (Outer[Times, pb[[c]], pbmx] - Outer[Times, pbmx, pb[[c]]]
     + Outer[Times, p[[c]], pmx] - Outer[Times, pmx, p[[c]]])], {c, 4}]];
Partner[k_Integer] := If[k <= 2, k + 2, k - 2];
LinearCurvature[gam_, hbg_] := Module[{r4, ric, p, pb, pup, pbup, sc},
 r4 = Table[DB[gam[[b]], a] - DB[gam[[a]], b], {a, 4}, {b, 4}];
 ric = Table[Sum[(r4[[c, b, Partner[c], a]] + r4[[Partner[c], a, c, b]])/2,
   {c, 4}], {a, 4}, {b, 4}];
 p = (jmat + hbg)/2; pb = (jmat - hbg)/2; pup = jmat . p . jmat; pbup = jmat . pb . jmat;
 sc = TU[Sum[(pup[[a, c]] pup[[b, d]] - pbup[[a, c]] pbup[[b, d]])
   (r4[[c, d, a, b]] + r4[[a, b, c, d]])/2, {a, 4}, {b, 4}, {c, 4}, {d, 4}]];
 {sc, TUmat[p . jmat . ric . jmat . pb]}];

(* Constant-state saddle and aligned frame to O(u). R data are derived from
   exact g,B,e,ebar (SMexactRdata). NR data expand NRHcompact and SMvielbein:
   chi=2 Sqrt[Lp Lm] u+O(u^3), W=W1 u+O(u^3). The polynomial form also
   covers Lp=0 or Lm=0 without dividing by either constant. *)
IntrinsicBackground[branch_String] := Module[{g, bb, e, eb, gi, hh, vv, vb},
 If[branch === "R",
   g = {{2 Lp, -1/u - Lp Lm u}, {-1/u - Lp Lm u, 2 Lm}};
   bb = {{0, -1/u - Lp Lm u}, {1/u + Lp Lm u, 0}};
   e = {{1, -Lp}, {-u Lm, 1/u}}; eb = {{1/u, -u Lp}, {-Lm, 1}};
   gi = Inverse[g]; hh = ArrayFlatten[{{gi, -gi . bb}, {bb . gi, g - bb . gi . bb}}];
   vv = Join[Transpose[Inverse[e]], e . eta + bb . Transpose[Inverse[e]]]/Sqrt[2];
   vb = Join[Transpose[Inverse[eb]], eb . etab + bb . Transpose[Inverse[eb]]]/Sqrt[2];
   {hh, vv, vb} = Map[Normal[Series[#, {u, 0, 1}]] &, {hh, vv, vb}, {3}],
   hh = {{0, 0, 1, -2 Lm u}, {0, 0, 2 Lp u, -1},
     {1, 2 Lp u, 0, W1 u}, {-2 Lm u, -1, W1 u, 0}};
   vv = {{1/Sqrt[2], 0}, {Lp u/Sqrt[2], 0}, {0, -Sqrt[2]},
     {W1 u/(2 Sqrt[2]), Sqrt[2] Lm u}};
   vb = {{0, Lm u/Sqrt[2]}, {0, 1/Sqrt[2]}, {-Sqrt[2] Lp u, -W1 u/(2 Sqrt[2])},
     {Sqrt[2], 0}}];
 {hh, vv, vb}];

(* SMRstresspluscoefficients, SMRstressminuscoefficients, SMRtypecoefficients,
   SMRhaircoefficients, SMRdilatoncoefficients, SMRintegratedstress;
   SMNRcompactsolution, SMNRcompactU, SMNRintegratedstress, constant Lpm.
   These supply all O(u) terms, including logarithmic radial responses.
   Omitting them generally changes the finite match and the obstruction.
   The R integration constant cs here is -4 f0 in the manuscript. *)
SourceMap[branch_String] := Module[
 {a = a0[xp, xm], b = b0[xp, xm], r = r0[xp, xm], c = c0[xp, xm], v = v0[xp, xm],
  h0, h1, h1b, h2, h2l, d1, d2, d2l, fp, fm, cuts, radial, dcut, drad},
 h0 = {{b, c}, {r, a}};
 h1 = {{-l/2 Der[r, 2, 0], l/2 (4 Der[v, 1, 1] - Der[a, 2, 0] - Der[b, 0, 2] - l^2/4 Der[r, 2, 2])},
   {0, -l/2 Der[r, 0, 2]}};
 h1b = {{0, l^2/8 Der[r, 2, 2]}, {0, 0}}; d1 = -l/8 Der[r, 1, 1];
 If[branch === "R",
   h2 = {{rp[xp, xm], hs[xp, xm]}, {cs + 4 v - l^2/4 Der[r, 1, 1], rm[xp, xm]}};
   h2l = {{l/2 Lp Der[r, 1, 1] - l^3/8 Der[r, 3, 1],
     l/2 Der[rm[xp, xm], 2, 0] + l/2 Der[rp[xp, xm], 0, 2] + l^3 Der[v, 2, 2]
     - l^3/4 Der[a, 3, 1] + l^5/16 Der[r, 3, 3] - l^3/4 Der[b, 1, 3]
     + l Lp Der[a, 1, 1] + l Lm Der[b, 1, 1] - 2 l Lp Der[v, 0, 2] - 2 l Lm Der[v, 2, 0]
     - l^3/4 Lm Der[r, 3, 1] - l/2 Lp Lm Der[r, 1, 1] - l^3/4 Lp Der[r, 1, 3]},
     {-l/2 Der[r, 1, 1], l/2 Lm Der[r, 1, 1] - l^3/8 Der[r, 1, 3]}};
   d2 = l^2/8 Der[a, 2, 0] - l^2/2 Der[v, 1, 1] - l^4/32 Der[r, 2, 2]
     + l^2/8 Lp Der[r, 0, 2] + l^2/8 Lm Der[r, 2, 0] + l^2/8 Der[b, 0, 2];
   d2l = -l^3/16 Der[r, 2, 2];
   fp = -2 Lp Der[a, 1, 0] - l^2 Der[v, 2, 1] + 4 Lp Der[v, 0, 1]
     + l^2/4 Der[a, 3, 0] - l^4/8 Der[r, 3, 2] + l^2 Lp Der[r, 1, 2]
     + l^2/4 Der[b, 1, 2] - 2 Lp Lm Der[r, 1, 0] + l^2/4 Lm Der[r, 3, 0];
   fm = -2 Lm Der[b, 0, 1] - l^2 Der[v, 1, 2] + 4 Lm Der[v, 1, 0]
     + l^2/4 Der[b, 0, 3] - l^4/8 Der[r, 2, 3] + l^2/4 Der[a, 2, 1]
     + l^2 Lm Der[r, 2, 1] - 2 Lp Lm Der[r, 0, 1] + l^2/4 Lp Der[r, 0, 3],
   h2 = {{rp[xp, xm], hs[xp, xm]}, {cs, rm[xp, xm]}};
   h2l = {{l/2 Lp Der[r, 1, 1],
     l/2 Der[rm[xp, xm], 2, 0] + l/2 Der[rp[xp, xm], 0, 2]
     + l Lp Der[a, 1, 1] + l Lm Der[b, 1, 1] - 2 l Lp Der[v, 0, 2] - 2 l Lm Der[v, 2, 0]
     - l^3/4 Lm Der[r, 3, 1] - l^3/4 Lp Der[r, 1, 3] - l/8 W1 Der[r, 1, 1]},
     {0, l/2 Lm Der[r, 1, 1]}};
   d2 = l^2/8 (Lp Der[r, 0, 2] + Lm Der[r, 2, 0]); d2l = 0;
   fp = -2 Lp Der[a, 1, 0] + 4 Lp Der[v, 0, 1] - W1/2 Der[r, 1, 0]
     + l^2 Lp Der[r, 1, 2] + l^2/4 Lm Der[r, 3, 0];
   fm = -2 Lm Der[b, 0, 1] + 4 Lm Der[v, 1, 0] - W1/2 Der[r, 0, 1]
     + l^2 Lm Der[r, 2, 1] + l^2/4 Lp Der[r, 0, 3]];
 cuts = h0 + Y h1 + Y^2 h1b + u (h2 + Y h2l);
 radial = h1 + 2 Y h1b + u (-2/l (h2 + Y h2l) + h2l);
 dcut = v + Y d1 + u (d2 + Y d2l); drad = d1 + u (-2/l (d2 + Y d2l) + d2l);
 <|"Fields" -> {cuts[[1, 1]], cuts[[2, 2]], cuts[[2, 1]], cuts[[1, 2]], dcut},
   "Radial" -> {radial[[1, 1]], radial[[2, 2]], radial[[2, 1]], radial[[1, 2]], drad},
   "Fp" -> fp, "Fm" -> fm|>];
SubSources[e_, sm_] := Module[{heads, values},
 heads = Join[freeFields, {gPP, gMM, gMP, gPM, gd}];
 values = Join[sm["Fields"], sm["Radial"]];
 Expand[e /. {
   Derivative[m_Integer, n_Integer][ff_][xp, xm] /; MemberQ[heads, ff] :>
     Der[values[[First[FirstPosition[heads, ff]]]], m, n],
   ff_[xp, xm] /; MemberQ[heads, ff] :> values[[First[FirstPosition[heads, ff]]]]}]];
UseConstraints[e_, sm_] := e /. {
 Derivative[m_Integer, n_Integer][rp][xp, xm] /; n >= 1 :> Der[sm["Fp"], m, n - 1],
 Derivative[m_Integer, n_Integer][rm][xp, xm] /; m >= 1 :> Der[sm["Fm"], m - 1, n]};
ReduceSources[e_, sm_] := IBP[UseConstraints[IBP[Nondecay[SubSources[e, sm]], fields], sm], fields] // Expand;

(* Derive each branch once at arbitrary constant state, then specialize to its vacuum. *)
RunBranch[branch_String] := Module[
 {data, hh0, vv, vb, hf, xx, dh, dh2, he, sc, s1, s2, gam, curv, qconn, qdirect,
  p, pb, hu0, om, pdhpb, boxmat, qflat, qup, qsq, sboxs, linLag,
  measure, cand, ct, unsub, sm, diff, rawren, vac, div, expected, finiteExpected,
  vacuumRules = {Lp -> 0, Lm -> 0, W1 -> 0}, a, b, r, v, zero = ConstantArray[0, {4, 4}]},
 Print["\n=== ", branch, " branch: intrinsic free-cutoff calculation ==="];
 data = IntrinsicBackground[branch];
 If[!NRH`Check[branch <> " background/frame matrix dimensions",
   ListQ[data] && Length[data] == 3 && And @@ (MatrixQ /@ data) &&
   (Dimensions /@ data) === {{4, 4}, {4, 2}, {4, 2}}], NRH`FileSummary[]];
 {hh0, vv, vb} = data;
 NRH`CheckZero[branch <> " background H J H = J through O(u)", TUmat[hh0 . jmat . hh0 - jmat]];
 NRH`CheckZero[branch <> " frame projector completeness", TUmat[vv . eta . Transpose[vv] + vb . etab . Transpose[vb] - jmat]];
 NRH`CheckZero[branch <> " frame reconstruction", TUmat[vv . eta . Transpose[vv] - vb . etab . Transpose[vb] - hh0]];
 hf = {{fPP[xp, xm], fPM[xp, xm]}, {fMP[xp, xm], fMM[xp, xm]}};
 xx = vv . eta . hf . etab . Transpose[vb]; dh = TUmat[xx + Transpose[xx]];
 dh2 = TUmat[-1/2 hh0 . jmat . dh . jmat . dh]; he = hh0 + eps dh + eps^2 dh2;
 NRH`CheckZero[branch <> " coset quadratic completion", TUmat[TEmat[he . jmat . he - jmat]]];
 NRH`CheckZero[branch <> " fluctuation/frame round trip", TUmat[Transpose[jmat . vv] . dh . jmat . vb - hf]];
 p = (jmat + hh0)/2; pb = (jmat - hh0)/2; hu0 = jmat . hh0 . jmat;
 NRH`CheckZero[branch <> " quadratic completion has no mixed projection", TUmat[p . jmat . dh2 . jmat . pb]];
 sc = TU[TE[ScalarDFT[he, -Y/l + eps fd[xp, xm]]]];
 NRH`CheckZero[branch <> " intrinsic background scalar vanishes", sc /. eps -> 0];
 s1 = Coefficient[sc, eps, 1]; s2 = Coefficient[sc, eps, 2];
 Print["S1 = ", InputForm[s1]]; Print["S2 = ", InputForm[s2]];
 gam = GammaLinear[hh0, dh, fd[xp, xm]]; curv = LinearCurvature[gam, hh0];
 NRH`CheckZero[branch <> " scalar from Christoffel curvature equals explicit DFT scalar at linear order", s1 - curv[[1]]];
 qconn = curv[[2]];
 om = Table[Sum[jmat[[n, k]] DB[dh[[n, s]], k], {n, 4}, {k, 4}]
   - 2 Sum[hh0[[s, k]] jmat[[k, n]] DB[fd[xp, xm], n], {k, 4}, {n, 4}], {s, 4}];
 pdhpb = TUmat[p . jmat . dh . jmat . pb];
 boxmat = Sum[hu0[[i, j]] DB[DB[pdhpb, i], j], {i, 4}, {j, 4}];
 qdirect = TUmat[Table[-boxmat[[m, n]]/4 + 1/2 Sum[
   (p . jmat)[[m, i]] (pb . jmat)[[n, s]] DB[om[[s]], i]
   - (p . jmat)[[m, s]] (pb . jmat)[[n, i]] DB[om[[s]], i], {i, 4}, {s, 4}], {m, 4}, {n, 4}]];
 NRH`CheckZero[branch <> " mixed Ricci: Christoffel vs independent variation formula", qconn - qdirect];
 qflat = TUmat[Transpose[jmat . vv] . qconn . jmat . vb];
 Print["S_(a,bbar) at order eps = ", InputForm[qflat]];
 linLag = Expand[(s2 - 2 fd[xp, xm] s1)/u];
 qup = Table[VariationalD[linLag, hf[[i, j]]] u/2, {i, 2}, {j, 2}];
 NRH`CheckZero[branch <> " mixed Ricci: independent Euler derivative of scalar Hessian", TUmat[eta . qup . etab] - qflat];
 qsq = TU[Sum[eta[[i, k]] etab[[j, q]] qflat[[i, j]] qflat[[k, q]], {i, 2}, {j, 2}, {k, 2}, {q, 2}]];
 NRH`CheckZero[branch <> " Q_AB Q^AB agrees with flat contraction",
   qsq - TU[Sum[qconn[[i, j]] (jmat . qconn . jmat)[[i, j]], {i, 4}, {j, 4}]]];
 sboxs = TU[s1 Sum[hu0[[i, j]] DB[DB[s1, i], j], {i, 4}, {j, 4}]];
 Print["(S1)^2 + 8 (Q1)^2 = ", InputForm[TU[s1^2 + 8 qsq]]];
 Print["S1 Box0 S1 = ", InputForm[sboxs]];
 (* S0=Q0=0. The quadratic S Box S term needs only Box0, while l S
    needs both S2 and the linear measure correction -2 fd S1. *)
 measure = (1 - 2 eps fd[xp, xm] + 2 eps^2 fd[xp, xm]^2)/u;
 NRH`CheckZero[branch <> " SMvolumecounterterm quadratic coefficient", Coefficient[Expand[(8/l) measure], eps, 2] - 16 fd[xp, xm]^2/(l u)];
 cand = Nondecay[Coefficient[Expand[measure (l (eps s1 + eps^2 s2)
   + l^3/16 eps^2 TU[s1^2] - l^2 Y/8 eps^2 TU[s1^2 + 8 qsq]
   + l^3 Y^2/32 eps^2 sboxs)], eps, 2]];
 a = fMM[xp, xm]; b = fPP[xp, xm]; r = fMP[xp, xm]; v = fd[xp, xm];
 ct = (l/2 a Der[r, 2, 0] + l/2 b Der[r, 0, 2] - 2 l v Der[r, 1, 1]
   + l^3/16 r Der[r, 2, 2])/u + 2 Y (Lp a + Lm b) Der[r, 1, 1]
   + l Y^2/2 r (Lp Der[r, 1, 3] + Lm Der[r, 3, 1]);
 If[branch === "R", ct += l^2 Y v Der[r, 2, 2] + l^3 Y^2/16 r Der[r, 3, 3]];
 Print["SM29 minus common volume, free-cutoff quadratic density (mod parts) = ", InputForm[IBP[cand, freeFields]]];
 Print["SM29 - (volume + branch source counterterm), free-cutoff density = ", InputForm[IBP[cand - ct, freeFields]]];
 unsub = -((-(a gPP[xp, xm] + b gMM[xp, xm] + r gPM[xp, xm] + fPM[xp, xm] gMP[xp, xm])/2
   + 8 v gd[xp, xm])/u - 4/l v (Lp a + Lm b + If[branch === "R", Lp Lm r + fPM[xp, xm], W1 r/4]));
 sm = SourceMap[branch];
 NRH`CheckZero[branch <> " source substitution preserves undifferentiated cutoff field",
   SubSources[fMM[xp, xm], sm] - sm["Fields"][[2]]];
 NRH`CheckZero[branch <> " source substitution commutes with tangential derivatives",
   SubSources[Der[fMP[xp, xm], 2, 1], sm] - Der[sm["Fields"][[3]], 2, 1]];
 diff = ReduceSources[cand - ct, sm]; rawren = ReduceSources[unsub + ct, sm];
 NRH`CheckZero[branch <> " branch quadratic counterterm cancels raw UV divergence (constant state)", DivPart[rawren]];
 NRH`CheckZero[branch <> " branch quadratic counterterm cancels raw UV divergence (vacuum)", DivPart[rawren /. vacuumRules]];
 vac = Expand[diff /. vacuumRules];
 Print[branch, " vacuum: candidate - (volume + source counterterm) = ", InputForm[vac]];
 NRH`CheckZero[branch <> " vacuum: SM29 divergent match", DivPart[vac]];
 If[branch === "NR", NRH`CheckZero["NR vacuum: SM29 finite and divergent equality modulo parts", vac],
   a = a0[xp, xm]; b = b0[xp, xm]; r = r0[xp, xm]; v = v0[xp, xm];
   finiteExpected = l/8 (l^2 a Der[r, 3, 1] + l^2 b Der[r, 1, 3] - 8 l^2 r Der[v, 2, 2]
     - 16 a Der[v, 2, 0] - 16 b Der[v, 0, 2] + 64 v Der[v, 1, 1]);
   NRH`CheckZero["R vacuum: explicitly computed finite local contact mismatch", vac - finiteExpected];
   NRH`Check["R vacuum: finite local mismatch is not zero", !NRH`ZeroQ[vac]]];
 a = a0[xp, xm]; b = b0[xp, xm]; r = r0[xp, xm];
 expected = -2 Y (Lp a + Lm b) Der[r, 1, 1]
   + (l Y^2/2 - l^2 Y/4) r (Lp Der[r, 1, 3] + Lm Der[r, 3, 1]);
 div = DivPart[diff];
 Print[branch, " constant-state divergent obstruction Drem = ", InputForm[div]];
 NRH`CheckZero[branch <> " constant-state obstruction equals displayed Drem", div - expected];
 NRH`Check[branch <> " generic nonzero Lpm obstructs SM29 cancellation", !NRH`ZeroQ[div]];
 NRH`Check[branch <> " divergent obstruction independent of W1", FreeQ[div, W1]];
 NRH`CheckZero[branch <> " no uncanceled power divergence in candidate difference", Coefficient[Expand[u diff], u, 0]];
 <|"Difference" -> diff, "VacuumDifference" -> vac, "Divergence" -> div,
   "ScalarLinear" -> s1, "ScalarQuadratic" -> s2, "MixedRicciLinear" -> qflat,
   "CandidateFree" -> cand, "SourceCountertermFree" -> ct|>];

nrResult = RunBranch["NR"];
rResult = RunBranch["R"];
NRH`CheckZero["same generic constant-state obstruction on R and NR", nrResult["Divergence"] - rResult["Divergence"]];
(* Independent compact-support/torus witness. Choose a0=A Cos[theta], r0=R Cos[theta],
   b0=0. Its independent A R term cannot cancel a term proportional to R^2.
   Taking xp,xm of period 2 Pi and kplus=kminus=1 gives a nonzero integral. *)
witnessDensity = nrResult["Divergence"] /. {
 a0 -> Function[{xp, xm}, aa Cos[xp + xm]],
 b0 -> Function[{xp, xm}, 0], r0 -> Function[{xp, xm}, rr Cos[xp + xm]]};
mixedWitness = aa rr Coefficient[Coefficient[Expand[witnessDensity], aa, 1], rr, 1];
NRH`CheckZero["cosine witness obtained from computed obstruction", mixedWitness - 2 Y Lp aa rr Cos[xp + xm]^2];
witnessIntegral = Integrate[mixedWitness, {xp, 0, 2 Pi}, {xm, 0, 2 Pi}];
Print["Periodic A*R obstruction witness integral = ", InputForm[witnessIntegral]];
NRH`CheckZero["obstruction is not a tangential total derivative: periodic witness", witnessIntegral - 4 Pi^2 Y Lp aa rr];
NRH`Check["periodic witness is generically nonzero", !NRH`ZeroQ[witnessIntegral]];
Print["Conclusion: vacuum quadratic matching only. This does not establish a general-state or finite-source nonlinear counterterm."];
Print["R contact adjustment: subtract the printed finite mismatch density from the SM29 quadratic density to match the branch source-counterterm convention."];
NRH`FileSummary[];
