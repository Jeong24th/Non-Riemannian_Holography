(* ::Package:: *)

(* ::Title:: *)
(*NRH08 — SM 7, "Classical Boundary-Theory Candidate"*)


(* ::Text:: *)
(*This file verifies the classical statements of SM 7:*)
(**)
(*  SM (120)-(121) [SMcandidateD, SMcandidateaction]  the covariant form of the candidate*)
(*      action reduces, through the type-(1,1) boundary projectors, to*)
(*      Tr(D+ phi^+ D- phi^- + psi^+ D- psi^+ + psi^- D+ psi^-):  the P-projector term*)
(*      contains only D_+ derivatives and the Pbar term only D_-, each selecting one*)
(*      independent component of phi_A;*)
(*  SM (122) [SMcandidategauge]  invariance under the internal gauge transformations*)
(*      (verified with a genuinely non-abelian 2 x 2 internal algebra);*)
(*  SM (123)-(124) [SMcandidateWitt(L)]  the two commuting classical Witt symmetries:*)
(*      delta_v L = d_lambda (v^lambda L),  including the weight-1/2 fermions;*)
(*  the equivalence  Lhat_xi H^(0) = 0  <=>  chirality of v^pm  for pure boundary*)
(*      diffeomorphisms (quoted below SM (124));*)
(*  SM (125)-(126) [SMcandidatefermionic(L)]  the chiral fermionic transformations with*)
(*      arbitrary Grassmann-odd eps^pm(x^pm):  delta L is the displayed total derivative.*)
(*      Grassmann arithmetic is handled EXACTLY by realizing the four odd generators as*)
(*      Jordan-Wigner matrices, so anticommutation is ordinary matrix algebra;*)
(*  SM (127)-(128) [SMcandidateextra(L)]  the extra Grassmann-even chiral transformation*)
(*      and its total-derivative variation (non-abelian);*)
(*  SM (129)-(130) [SMcandidatetrivial(L)]  the equation-of-motion redundancies delta_alpha and*)
(*      their off-shell total-derivative variation, with two-component Grassmann fermions.*)


ClearAll["Global`*"];
Get[FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH08_SM_BoundaryCandidate.wl"];


(* ::Section:: *)
(*SM (120)-(121): projector selection on the type-(1,1) boundary*)


J4 = ODDJ[2];
H0 = {{0, 0, 1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, -1, 0, 0}};
P0up = J4 . ((J4 + H0)/2) . J4;
Pb0up = J4 . ((J4 - H0)/2) . J4;
(* on the section, D_A = (0, 0; D+, D-) *)
DofA = {0, 0, DP1, DP2};  (* formal derivative symbols acting on phi_B *)
term1 = Sum[P0up[[a, b]] DofA[[a]] phi[b], {a, 4}, {b, 4}];
term2 = Sum[Pb0up[[a, b]] DofA[[a]] phi[b], {a, 4}, {b, 4}];
NRH`Check["SM(121): the P^(0) term contains only D_+ and the Pbar^(0) term only D_-",
   FreeQ[term1, DP2] && FreeQ[term2, DP1] && ! FreeQ[term1, DP1] && ! FreeQ[term2, DP2]];
NRH`Check["SM(121): each term selects a single independent component of phi_A",
   Length[Union[Cases[term1, phi[_], Infinity]]] == 1 &&
   Length[Union[Cases[term2, phi[_], Infinity]]] == 1 &&
   Cases[term1, phi[_], Infinity] =!= Cases[term2, phi[_], Infinity]];
(* the 2d boundary gamma matrices single out one chiral fermion component each *)
g2op = Sqrt[2] {{0, 0}, {1, 0}}; g2om = -Sqrt[2] {{0, 1}, {0, 0}};
NRH`Check["SM(121): gamma^oplus (gamma^bar-ominus) have rank one: one surviving chiral component each",
   MatrixRank[g2op] == 1 && MatrixRank[g2om] == 1];


(* ::Section:: *)
(*SM (122): non-abelian internal gauge invariance (bosonic sector)*)


(* internal algebra: 2 x 2 matrix-valued fields; D_mu X = d_mu X - i[A_mu, X] *)
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
NRH`CheckZero["SM(122): delta_Lambda L = 0 for the non-abelian bosonic sector",
   Together[Expand[dLgauge]]];


(* ::Section:: *)
(*SM (123)-(124): the classical Witt symmetries (with weight-1/2 fermions)*)


(* ::Text:: *)
(*Grassmann-odd fields are realized exactly through Jordan-Wigner matrices: four odd*)
(*generators theta_i with {theta_i, theta_j} = 0; an arbitrary odd field is theta_i times*)
(*an ordinary function (sufficient by multilinearity of every identity checked below).*)
(*The internal algebra is taken abelian in the fermionic checks, as the transformations*)
(*SM (123), (125), (127) do not involve the internal commutator.*)


sp = {{0, 1}, {0, 0}}; sz = {{1, 0}, {0, -1}}; i2 = IdentityMatrix[2];
th[1] = KroneckerProduct[sp, i2, i2, i2];
th[2] = KroneckerProduct[sz, sp, i2, i2];
th[3] = KroneckerProduct[sz, sz, sp, i2];
th[4] = KroneckerProduct[sz, sz, sz, sp];
NRH`CheckZero["Jordan-Wigner generators anticommute (exact Grassmann algebra)",
   Flatten[Table[th[i] . th[j] + th[j] . th[i], {i, 4}, {j, 4}]]];
id16 = IdentityMatrix[16];

psiP = th[1] pf[xp, xm];    (* psi^+ *)
psiM = th[2] mf[xp, xm];    (* psi^- *)
phiPa = id16 bp[xp, xm];    (* abelian phi^+, phi^- as even elements *)
phiMa = id16 bm[xp, xm];
LC = phiHold;   (* assembled below *)
Lcand[phP_, phM_, psP_, psM_] :=
   D[phP, xp] . D[phM, xm] + psP . D[psP, xm] + psM . D[psM, xp];
L0 = Lcand[phiPa, phiMa, psiP, psiM];

(* Witt transformation with v^+(x^+), v^-(x^-) *)
vL = {vpf[xp], vmf[xm]};
wD[X_] := vL[[1]] D[X, xp] + vL[[2]] D[X, xm];
dPhiP = wD[phiPa]; dPhiM = wD[phiMa];
dPsiP = wD[psiP] + 1/2 D[vL[[1]], xp] psiP;
dPsiM = wD[psiM] + 1/2 D[vL[[2]], xm] psiM;
dLWitt = D[Lcand[phiPa + t dPhiP, phiMa + t dPhiM, psiP + t dPsiP, psiM + t dPsiM], t] /. t -> 0;
NRH`CheckZero["SM(124): delta_v L = d_+(v^+ L) + d_-(v^- L)  (unit-weight scalar density)",
   Map[Together, Expand[dLWitt - D[vL[[1]] L0, xp] - D[vL[[2]] L0, xm]], {2}]];

(* boundary isometry statement: for xi^A = (0; v^mu),  Lhat_xi H^(0) = 0  iff  v^pm chiral *)
xs2 = {xp, xm};
bDblD[e_, m_] := If[m <= 2, 0, D[e, xs2[[m - 2]]]];
GenLie2[xiUp_, HH_] := Module[{dim = 4, xiLow, amat},
   xiLow = J4 . xiUp;
   amat = Table[bDblD[xiUp[[c]], m] - Sum[J4[[c, dd]] bDblD[xiLow[[m]], dd], {dd, 4}], {m, 4}, {c, 4}];
   Sum[xiUp[[c]] bDblD[HH, c], {c, 4}] + amat . HH + HH . Transpose[amat]];
lie2 = GenLie2[{0, 0, vg1[xp, xm], vg2[xp, xm]}, H0];
NRH`Check["Lhat_xi H^(0) = 0  <=>  d_- v^+ = 0 = d_+ v^-",
   Module[{eqs = DeleteCases[Union[Flatten[lie2]], 0]},
      Union[Together[eqs /. {Derivative[0, 1][vg1][xp, xm] -> DV1, Derivative[1, 0][vg2][xp, xm] -> DV2}]] ===
      Union[{2 DV1, -2 DV1, 2 DV2, -2 DV2}] ||
      (And @@ (PossibleZeroQ[# /. {Derivative[0, 1][vg1][xp, xm] -> 0, Derivative[1, 0][vg2][xp, xm] -> 0}] & /@ eqs)) &&
      ! FreeQ[eqs, Derivative[0, 1][vg1][xp, xm]] && ! FreeQ[eqs, Derivative[1, 0][vg2][xp, xm]]]];


(* ::Section:: *)
(*SM (125)-(126): the chiral fermionic transformations*)


epsP = th[3] ef[xp];   (* Grassmann-odd eps^+(x^+) *)
epsM = th[4] gf[xm];   (* Grassmann-odd eps^-(x^-) *)

(* delta_{eps+}: delta phi^- = psi^+ eps^+, delta psi^+ = (1/2) eps^+ D_+ phi^+, rest zero *)
dphiM1 = psiP . epsP;
dpsiP1 = 1/2 epsP . D[phiPa, xp];
dL1 = D[Lcand[phiPa, phiMa + t dphiM1, psiP + t dpsiP1, psiM], t] /. t -> 0;
NRH`CheckZero["SM(126): delta_{eps+} L = d_-( psi^+ delta_{eps+} psi^+ )",
   Map[Together, Expand[dL1 - D[psiP . dpsiP1, xm]], {2}]];

(* delta_{eps-}: delta phi^+ = psi^- eps^-, delta psi^- = (1/2) eps^- D_- phi^-, rest zero *)
dphiP2 = psiM . epsM;
dpsiM2 = 1/2 epsM . D[phiMa, xm];
dL2 = D[Lcand[phiPa + t dphiP2, phiMa, psiP, psiM + t dpsiM2], t] /. t -> 0;
NRH`CheckZero["SM(126): delta_{eps-} L = d_+( psi^- delta_{eps-} psi^- )",
   Map[Together, Expand[dL2 - D[psiM . dpsiM2, xp]], {2}]];
NRH`Check["the fermionic parameters carry arbitrary chiral profiles (one function per chirality)",
   ! FreeQ[dL1, ef] && ! FreeQ[dL2, gf]];


(* ::Section:: *)
(*SM (127)-(128): the extra Grassmann-even chiral transformation (non-abelian)*)


(* delta_zeta phi^+ = zetabar(x^-) D_- phi^-,  delta_zeta phi^- = zeta(x^+) D_+ phi^+ *)
dphiPz = zb[xm] DD[phiM, 2];
dphiMz = zf[xp] DD[phiP, 1];
dLz = D[Tr[DD2[phiP + t dphiPz, 1] . DD2[phiM + t dphiMz, 2]] /.
   DD2[X_, mu_] :> (D[X, {xp, xm}[[mu]]] - I (({Ap, Am}[[mu]]) . X - X . ({Ap, Am}[[mu]]))), t] /. t -> 0;
NRH`CheckZero["SM(128): delta_zeta L = (1/2) d_+ Tr[zetabar (D_- phi^-)^2] + (1/2) d_- Tr[zeta (D_+ phi^+)^2]",
   Together[Expand[dLz
      - 1/2 D[zb[xm] Tr[DD[phiM, 2] . DD[phiM, 2]], xp]
      - 1/2 D[zf[xp] Tr[DD[phiP, 1] . DD[phiP, 1]], xm]]]];


(* ::Section:: *)
(*SM (129)-(130) [SMcandidatetrivial(L)]: the equation-of-motion redundancies*)


(* ::Text:: *)
(*For arbitrary Grassmann-even local functions alpha_0, alpha_pm, the transformations*)
(*   delta psi^+ = alpha_0 D_+ psi^- + alpha_+ D_- psi^+,   delta psi^- = alpha_0 D_- psi^+ + alpha_- D_+ psi^-*)
(*change the Lagrangian by the displayed total derivative, off shell and identically.  Here each*)
(*fermion carries two Jordan-Wigner components, so that the Grassmann cancellations*)
(*Tr(d_- psi^+ d_- psi^+) = 0 and Tr(d_+ psi^- d_- psi^+ + d_- psi^+ d_+ psi^-) = 0 are genuinely*)
(*tested rather than trivialized by the nilpotency of a single generator.*)


psiP2 = th[1] pf1[xp, xm] + th[3] pf3[xp, xm];
psiM2 = th[2] mf2[xp, xm] + th[4] mf4[xp, xm];
a0 = af0[xp, xm]; aPl = afp[xp, xm]; aMi = afm[xp, xm];
dpsiPa = a0 D[psiM2, xp] + aPl D[psiP2, xm];
dpsiMa = a0 D[psiP2, xm] + aMi D[psiM2, xp];
dLa = D[Lcand[phiPa, phiMa, psiP2 + t dpsiPa, psiM2 + t dpsiMa], t] /. t -> 0;
NRH`CheckZero["SM(130): delta_alpha L = d_-(psi^+ delta_alpha psi^+) + d_+(psi^- delta_alpha psi^-) identically (off shell)",
   Map[Together, Expand[dLa - D[psiP2 . dpsiPa, xm] - D[psiM2 . dpsiMa, xp]], {2}]];
NRH`Check["SM(129): the transformations are built from the fermion equations of motion d_- psi^+ = 0 = d_+ psi^-, hence vanish on shell",
   Module[{onshell = {Derivative[0, 1][pf1][xp, xm] -> 0, Derivative[0, 1][pf3][xp, xm] -> 0,
       Derivative[1, 0][mf2][xp, xm] -> 0, Derivative[1, 0][mf4][xp, xm] -> 0}},
      Union[Flatten[{dpsiPa, dpsiMa} /. onshell]] === {0}]];
NRH`Check["the two-component fermions used here are genuinely Grassmann: psi^+ d_- psi^+ is nonzero off shell",
   ! (Union[Flatten[Expand[psiP2 . D[psiP2, xm]]]] === {0})];

NRH`FileSummary[];
