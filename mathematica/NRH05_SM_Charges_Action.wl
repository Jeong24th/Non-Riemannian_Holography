(* ::Package:: *)

(* ::Title:: *)
(*NRH05 — SM 2 "Covariant Charges and Asymptotic Algebras" and SM 3 "Renormalized On-Shell Action"*)


(* ::Text:: *)
(*This file verifies, in paper order:*)
(**)
(*  SM 2 (charges):*)
(*   - the Riemannian surface-potential components SM (50) [SMRpotentialcomponents]:*)
(*        lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+'',  and its minus mirror;*)
(*   - the Brown-Henneaux central-charge identity behind Eq. (9)/SM (50): the charge*)
(*        cocycle reduces to a total derivative plus -(l^2/4) eps1 d^3 eps2;*)
(*   - the non-Riemannian charge one-form SM (52) [SMCPSresult]:*)
(*        k^{-y}[eps+] = (4/l) eps+ dL+  (and mirror), with all Thetahat components*)
(*        vanishing at the boundary;*)
(*   - the state-dependent falloffs SM (53) [SMNRchargefalloffs];*)
(*   - the componentwise W_1 cancellation SM (54) [SMNRchargecancellation]:  the charge*)
(*        one-form contains neither W_1 nor delta W_1 nor the opposite-chirality delta L;*)
(*   - the centerless algebra SM (56) [SMCPSalgebra]: same-chirality C-brackets close on*)
(*        the vector representative up to closed B-gauge parameters of vanishing potential,*)
(*        and the NR cocycle density is a total derivative (c_charge = 0).*)
(**)
(*  SM 3 (action):*)
(*   - the Gamma^2 identity SM (57) [SMgamma2]:  e^{-2d} S_(0) = L_{Gamma^2} + d_M(e^{-2d} B^M),*)
(*        verified exactly on BOTH saddles (arbitrary chiral L_pm; arbitrary W_0, W_1 hair*)
(*        on the non-Riemannian side);*)
(*   - the flux SM (58) [SMgamma2flux]:  B^y = 4 d_y d,  e^{-2d} B^y = -2 d_y e^{-2d}*)
(*        = -(4/l)(e^{2y/l} + mu e^{-2y/l})  with  mu = L+L- (R) and mu = L+L-/2 (NR),*)
(*        independently of the hair;*)
(*   - the cutoff computation SM (60) [SMgamma2cutoff] and the renormalized value SM (61);*)
(*   - the endpoint loci: Killing horizon u^2 = L+L- (R) vs q = 1, u^2 = L+L-/2 (NR).*)


ClearAll["Global`*"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH05_SM_Charges_Action.wl"];

JJ = ODDJ[3];


(* ::Section:: *)
(*Shared machinery: Noether potential K^{AB}[X] and boundary vector B^A*)


(* ::Text:: *)
(*K^{AB}[X] is the DFT Noether surface potential (Park-Rey-Rim-Sakatani, arXiv:1507.07545*)
(*Eq. (A.4)); B^A = 4 H^{AB} d_B d - d_B H^{AB} is the Gamma^2 boundary vector, and*)
(*Khat^{AB} = K^{AB} + 2 X^{[A} B^{B]}.  The implementation below is a direct port of the*)
(*machine-verified routine in the published Python archive (dft_covariant_phase_space.py).*)


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


(* ::Section:: *)
(*SM (50): Riemannian surface-potential components and the central charge*)


xsU = {xp, xm, Function[e, (2 u/l) D[e, u]]};
fB = u + Lp[xp] Lm[xm]/u;
gB = {{2 Lp[xp], -fB, 0}, {-fB, 2 Lm[xm], 0}, {0, 0, 1}};
BB = fB {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
HR = Map[Together, RiemannianH[gB, BB], {2}];
dR = -1/2 Log[u (1 - Lp[xp] Lm[xm]/u^2)];

(* full Eq. (8) generators, including the radial tails of xi^pm (Riemannian phase space) *)
xiPlus = {0, l^2/(2 u) Lm[xm] D[ep[xp], {xp, 2}], -l/2 D[ep[xp], xp],
   ep[xp], l^2/(4 u) D[ep[xp], {xp, 2}], -l/2 D[ep[xp], xp]};
xiMinus = {-l^2/(2 u) Lp[xp] D[em[xm], {xm, 2}], 0, +l/2 D[em[xm], xm],
   l^2/(4 u) D[em[xm], {xm, 2}], em[xm], -l/2 D[em[xm], xm]};

KfullP = Together[Exp[-2 dR] KhatComp[HR, dR, xiPlus, 5, 6, xsU]];
NRH`CheckZero["SM(50): lim e^{-2d} Khat^{-y}[eps+] = (4/l) eps+ L+ - 2 l eps+''",
   Together[Limit[KfullP, u -> Infinity]
      - (4/l ep[xp] Lp[xp] - 2 l D[ep[xp], {xp, 2}])]];
(* the display quotes Khat^{y+}; we compute the component pair (y, x+) accordingly *)
KfullM = Together[Exp[-2 dR] KhatComp[HR, dR, xiMinus, 6, 4, xsU]];
NRH`CheckZero["SM(50): lim e^{-2d} Khat^{y+}[eps-] = -(4/l) eps- L- + 2 l eps-''  (mirror)",
   Together[Limit[KfullM, u -> Infinity]
      - (-(4/l) em[xm] Lm[xm] + 2 l D[em[xm], {xm, 2}])]];

(* Brown-Henneaux cocycle: the charge-bracket density is
     e1 delta_{e2} L - [e1, e2] L  =  d/dx(e1 e2 L) - (l^2/4) e1 e2''' ,
   so after adding back (l^2/4) e1 e2''' the density must be a total x^+ derivative.
   A one-variable density is a total derivative iff all its Euler-Lagrange derivatives
   vanish; we test them for L and for both parameters.  The leftover -(l^2/4) e1 e2'''
   is the Brown-Henneaux center, normalized in NRH02 to c = 3l/2G. *)
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


(* ::Section:: *)
(*SM (52)-(54): the non-Riemannian charge one-form*)


(* ::Text:: *)
(*Following the published guard, the exact family is expanded through z^2 (z = e^{-2y/l}),*)
(*which retains every finite boundary term.  The charge one-form for the field-dependent*)
(*parameter X[eps] is*)
(*   k_X = delta(e^{-2d} Khat_X) - e^{-2d} Khat_{delta X} + 2 e^{-2d} X^{[A} Thetahat^{B]},*)
(*and the SM statements are: all e^{-2d} Thetahat^A vanish at the boundary; the finite*)
(*one-form is (4/l) eps^pm delta L_pm; W_1, delta W_1, and the opposite delta L drop out*)
(*componentwise.*)


xsZ = {xp, xm, Function[e, -(2 z/l) D[e, z]]};
HNRz = Module[{h = ConstantArray[0, {6, 6}]},
   h[[1, 4]] = h[[4, 1]] = 1 + 2 Lp[xp] Lm[xm] z^2;
   h[[2, 5]] = h[[5, 2]] = -1 - 2 Lp[xp] Lm[xm] z^2;
   h[[3, 3]] = h[[6, 6]] = 1;
   h[[1, 5]] = h[[5, 1]] = -2 Lm[xm] z;
   h[[2, 4]] = h[[4, 2]] = 2 Lp[xp] z;
   h[[4, 4]] = -2 Lp[xp] W1[xp, xm] z^2;
   h[[4, 5]] = h[[5, 4]] = W1[xp, xm] z;
   h[[5, 5]] = -2 Lm[xm] W1[xp, xm] z^2;
   h];
dNRz = Log[z]/2 + Lp[xp] Lm[xm] z^2/4;
eDenz = (1 - Lp[xp] Lm[xm] z^2/2)/z;

NRH`Check["truncation obeys H J H = J through z^2",
   Module[{c = Expand[HNRz . JJ . HNRz - JJ]},
      AllTrue[Flatten[c], PossibleZeroQ[Coefficient[#, z, 0]] && PossibleZeroQ[Coefficient[#, z, 1]] && PossibleZeroQ[Coefficient[#, z, 2]] &]]];
NRH`CheckZero["SM(53): state-dependent falloffs delta H^-_+ = 2 z dL+, delta H^+_- = -2 z dL-, delta H_{+-} = z dW1",
   {D[HNRz[[2, 4]], Lp[xp]] - 2 z, D[HNRz[[1, 5]], Lm[xm]] + 2 z, D[HNRz[[4, 5]], W1[xp, xm]] - z}];

(* variation along the state directions *)
varyRules = {Lp -> Function[x, Lp[x] + tt dLpF[x]], Lm -> Function[x, Lm[x] + tt dLmF[x]],
   W1 -> Function[{x, y2}, W1[x, y2] + tt dW1F[x, y2]]};
HNRzT = HNRz /. varyRules; dNRzT = dNRz /. varyRules; eDenzT = eDenz /. varyRules;
dH = D[HNRzT, tt] /. tt -> 0; dd0 = D[dNRzT, tt] /. tt -> 0;

(* Thetahat density:  e^{-2d} Theta^A - delta(e^{-2d} B^A),
   Theta^A = 4 H^{AB} d_B delta d - nabla_B delta H^{AB} *)
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
NRH`CheckZero["SM(52): lim e^{-2d} Thetahat^{+,-,y} = 0 at the boundary",
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
NRH`CheckZero["SM(52): k^{-y}[eps+] = (4/l) eps+ dL+",
   Together[kPlus - 4/l ep[xp] dLpF[xp]]];
NRH`CheckZero["SM(52): k^{+y}[eps-] = (4/l) eps- dL-",
   Together[kMinus - 4/l em[xm] dLmF[xm]]];
NRH`Check["SM(54): W_1, delta W_1, and the opposite-chirality delta L all drop out componentwise",
   FreeQ[{kPlus, kMinus}, W1] && FreeQ[{kPlus, kMinus}, dW1F] &&
   FreeQ[kPlus, dLmF] && FreeQ[kMinus, dLpF]];


(* ::Section:: *)
(*SM (55)-(56): bracket closure and the centerless algebra*)


CBracket[x_, y_, xs_] := Module[{xd = JJ . x, yd = JJ . y, dim = 6},
   Table[
      Sum[x[[b]] DblD[y[[a]], b, xs] - y[[b]] DblD[x[[a]], b, xs], {b, dim}]
      + 1/2 Sum[yd[[b]] Sum[JJ[[a, c]] DblD[x[[b]], c, xs], {c, dim}]
              - xd[[b]] Sum[JJ[[a, c]] DblD[y[[b]], c, xs], {c, dim}], {b, dim}],
      {a, dim}]];

xiPe = Function[{e}, {0, l^2 z Lm[xm] D[e, {xp, 2}]/2, -l D[e, xp]/2, e, 0, -l D[e, xp]/2}];
alphaP = e1[xp] D[e2[xp], xp] - e2[xp] D[e1[xp], xp];
bracketDiff = Together[CBracket[xiPe[e1[xp]], xiPe[e2[xp]], xsZ] - xiPe[alphaP]];
NRH`Check["SM(56): the same-chirality C-bracket closes up to a closed B-gauge parameter (slot x~+ only)",
   Together[bracketDiff[[2 ;; 6]]] === {0, 0, 0, 0, 0} && ! PossibleZeroQ[bracketDiff[[1]]]];
NRH`CheckZero["SM(56): the leftover reducibility parameter is chiral and closed:  d_- and d_y of it vanish",
   {D[bracketDiff[[1]], xm], D[bracketDiff[[1]], z]}];
NRH`CheckZero["SM(56): the closed B-gauge parameter carries no surface potential",
   Limit[Together[eDenz KhatComp[HNRz, dNRz, {zp[xp], 0, 0, 0, 0, 0}, 5, 6, xsZ]], z -> 0]];
(* centerless cocycle: without the anomalous term the density is a total derivative *)
NRH`CheckZero["SM(56): NR cocycle  e1 (e2 L' + 2 L e2') - alpha L  =  d/dx (e1 e2 L)   =>  c_charge = 0",
   Together[e1[xp] (e2[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e2[xp], xp])
      - alphaP Lp[xp] - D[e1[xp] e2[xp] Lp[xp], xp]]];


(* ::Section:: *)
(*SM (57)-(61): the Gamma^2 identity, the flux, and the renormalized action*)


(* Riemannian saddle, exact in u *)
gammaR = GammaDFT[HR, dR, xsU];
NRH`CheckZero["SM(57) on R:  e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)",
   Together[Exp[-2 dR] ScalarS0[HR, dR, xsU]
      - Gamma2Density[HR, dR, gammaR, xsU]
      - Sum[DblD[Exp[-2 dR] GammaBVector[HR, dR, xsU][[m]], m, xsU], {m, 6}]]];
BvecR = GammaBVector[HR, dR, xsU];
NRH`CheckZero["SM(58) on R:  B^y = 4 d_y d  and  e^{-2d}B^y = -(4/l)(u + L+L-/u)",
   {Together[BvecR[[6]] - 4 (2 u/l) D[dR, u]],
    Together[Exp[-2 dR] BvecR[[6]] + 4/l (u + Lp[xp] Lm[xm]/u)]}];

(* Non-Riemannian saddle, exact in chi with arbitrary W(x+, x-, chi) *)
chy = -(2 Sqrt[2]/l) Sinh[ch/Sqrt[2]];
chp = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psip][xp]/psip[xp];
chm = -Sqrt[2] Sinh[ch/Sqrt[2]] Derivative[1][psim][xm]/psim[xm];
xsNR = {Function[e, D[e, xp] + chp D[e, ch]], Function[e, D[e, xm] + chm D[e, ch]],
   Function[e, chy D[e, ch] + D[e, Ysym]]};
esig = psim[xm]/psip[xp];
HNRchi = {{0, 0, 0, Cosh[ch], -Sinh[ch]/esig, 0}, {0, 0, 0, esig Sinh[ch], -Cosh[ch], 0},
   {0, 0, 1, 0, 0, 0},
   {Cosh[ch], esig Sinh[ch], 0, -W[xp, xm, ch] esig Sinh[ch], W[xp, xm, ch] Cosh[ch], 0},
   {-Sinh[ch]/esig, -Cosh[ch], 0, W[xp, xm, ch] Cosh[ch], -W[xp, xm, ch] Sinh[ch]/esig, 0},
   {0, 0, 0, 0, 0, 1}};
dNRchi = -Ysym/l + Log[Cosh[ch/(2 Sqrt[2])]];
NRHZeroNR[label_, e_] := NRH`CheckZero[label,
   Together[ExpandAll[TrigToExp[e /. ch -> 2 Sqrt[2] Log[T]]] /. Log[T] -> LT]];

gammaNRc = GammaDFT[HNRchi, dNRchi, xsNR];
NRHZeroNR["SM(57) on NR (arbitrary W): e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)",
   Exp[-2 dNRchi] ScalarS0[HNRchi, dNRchi, xsNR]
      - Gamma2Density[HNRchi, dNRchi, gammaNRc, xsNR]
      - Sum[DblD[Exp[-2 dNRchi] GammaBVector[HNRchi, dNRchi, xsNR][[m]], m, xsNR], {m, 6}]];
BvecNR = GammaBVector[HNRchi, dNRchi, xsNR];
NRH`Check["SM(58)-(59) on NR: B^pm and B^y contain no W  (the hair never enters the flux)",
   FreeQ[Together[BvecNR], W]];
NRHZeroNR["SM(58) on NR: B^y = 4 d_y d",
   BvecNR[[6]] - 4 (chy D[dNRchi, ch] + D[dNRchi, Ysym])];
(* e^{-2d} = u - (L+L-/2)/u exactly, so e^{-2d}B^y = -(4/l)(u + mu/u), mu = L+L-/2 *)
eNRu = u - Lp[xp] Lm[xm]/(2 u);
NRH`CheckZero["SM(58) on NR: -2 d_y e^{-2d} = -(4/l)(u + (L+L-/2)/u)   [mu-dichotomy]",
   Together[-2 (2 u/l) D[eNRu, u] + 4/l (u + Lp[xp] Lm[xm]/(2 u))]];

(* SM (60)-(61): cutoff algebra and the renormalized value *)
SrenY = 1/(16 Pi G) (4/l (Exp[2 Y/l] + mu Exp[-2 Y/l]) - 8/l Sqrt[mu] - 4/l (Exp[2 Y/l] - mu Exp[-2 Y/l]));
NRH`CheckZero["SM(60): the regulated combination equals (8 mu/l) e^{-2Y/l} - (8/l) Sqrt[mu]",
   Together[SrenY - 1/(16 Pi G) (8 mu/l Exp[-2 Y/l] - 8/l Sqrt[mu])]];
NRH`CheckZero["SM(61): Y -> Infinity limit gives  S_ren = -(8 Sqrt[mu])/(16 pi G l) Int d^2x",
   Limit[SrenY, Y -> Infinity, Assumptions -> l > 0 && mu > 0] + 8 Sqrt[mu]/(16 Pi G l)];
NRH`CheckZero["endpoints: e^{-2d} = 0 at u^2 = L+L- (R horizon) and u^2 = L+L-/2 (NR, q = 1)",
   {Together[Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]],
    Together[eNRu /. u -> Sqrt[Lp[xp] Lm[xm]/2]]}];

NRH`FileSummary[];
