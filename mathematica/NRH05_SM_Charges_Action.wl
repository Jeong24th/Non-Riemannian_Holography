(* ::Title:: *)
(*NRH05 SM Charges Action*)

ClearAll["Global`*"];
Get[FileNameJoin[{If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]], "NRH01_DFT_Tools.wl"}]];
NRH`BeginFile["NRH05_SM_Charges_Action.wl"];

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

gammaR = GammaDFT[HR, dR, xsU];
NRH`CheckZero["on R: e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)",
   Together[Exp[-2 dR] ScalarS0[HR, dR, xsU]
      - Gamma2Density[HR, dR, gammaR, xsU]
      - Sum[DblD[Exp[-2 dR] GammaBVector[HR, dR, xsU][[m]], m, xsU], {m, 6}]]];
BvecR = GammaBVector[HR, dR, xsU];
NRH`CheckZero["on R: B^y = 4 d_y d and e^{-2d}B^y = -(4/l)(u + L+L-/u)",
   {Together[BvecR[[6]] - 4 (2 u/l) D[dR, u]],
    Together[Exp[-2 dR] BvecR[[6]] + 4/l (u + Lp[xp] Lm[xm]/u)]}];

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
NRHZeroNR["on NR (arbitrary W): e^{-2d} S_(0) = L_Gamma2 + d_M(e^{-2d} B^M)",
   Exp[-2 dNRchi] ScalarS0[HNRchi, dNRchi, xsNR]
      - Gamma2Density[HNRchi, dNRchi, gammaNRc, xsNR]
      - Sum[DblD[Exp[-2 dNRchi] GammaBVector[HNRchi, dNRchi, xsNR][[m]], m, xsNR], {m, 6}]];
BvecNR = GammaBVector[HNRchi, dNRchi, xsNR];
NRH`Check["on NR: B^pm and B^y contain no W (the hair never enters the flux)",
   FreeQ[Together[BvecNR], W]];
NRHZeroNR["on NR: B^y = 4 d_y d",
   BvecNR[[6]] - 4 (chy D[dNRchi, ch] + D[dNRchi, Ysym])];

eNRu = u - Lp[xp] Lm[xm]/(2 u);
NRH`CheckZero["on NR: -2 d_y e^{-2d} = -(4/l)(u + (L+L-/2)/u) [mu-dichotomy]",
   Together[-2 (2 u/l) D[eNRu, u] + 4/l (u + Lp[xp] Lm[xm]/(2 u))]];

SrenY = 1/(16 Pi G) (4/l (Exp[2 Y/l] + mu Exp[-2 Y/l]) - 8/l Sqrt[mu] - 4/l (Exp[2 Y/l] - mu Exp[-2 Y/l]));
NRH`CheckZero["the regulated combination equals (8 mu/l) e^{-2Y/l} - (8/l) Sqrt[mu]",
   Together[SrenY - 1/(16 Pi G) (8 mu/l Exp[-2 Y/l] - 8/l Sqrt[mu])]];
NRH`CheckZero["Y -> Infinity limit gives S_ren = -(8 Sqrt[mu])/(16 pi G l) Int d^2x",
   Limit[SrenY, Y -> Infinity, Assumptions -> l > 0 && mu > 0] + 8 Sqrt[mu]/(16 Pi G l)];
NRH`CheckZero["endpoints: e^{-2d} = 0 at u^2 = L+L- (R horizon) and u^2 = L+L-/2 (NR, q = 1)",
   {Together[Exp[-2 dR] /. u -> Sqrt[Lp[xp] Lm[xm]]],
    Together[eNRu /. u -> Sqrt[Lp[xp] Lm[xm]/2]]}];

NRH`CheckZero["footnote: with delta T = eps T' + 2 T eps' - l^2 eps''' the combination L + T/4 obeys the law with -(l^2/4) eps'''",
   Module[{dL = e1[xp] D[Lp[xp], xp] + 2 Lp[xp] D[e1[xp], xp],
      dT = e1[xp] D[TT[xp], xp] + 2 TT[xp] D[e1[xp], xp] - l^2 D[e1[xp], {xp, 3}], comb},
      comb = Lp[xp] + TT[xp]/4;
      Together[dL + dT/4 - (e1[xp] D[comb, xp] + 2 comb D[e1[xp], xp] - l^2/4 D[e1[xp], {xp, 3}])]]];

NRH`FileSummary[];
