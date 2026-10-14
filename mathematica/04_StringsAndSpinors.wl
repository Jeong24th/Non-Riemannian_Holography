(* Run this file in a fresh kernel: wolframscript -file 04_StringsAndSpinors.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* SM (58): eliminate beta and betaBar by their own field equations.
   dxPlus=partial x^+, barDxPlus=bar-partial x^+ (similarly -,y).
   F=exp(2y/l)+L_+ L_- exp(-2y/l) is nonzero. Overall 1/(2 Pi alpha') omitted. *)
firstOrder=(dy barDy+2 Lplus dxPlus barDxPlus+2 Lminus dxMinus barDxMinus
           +beta barDxPlus+betaBar dxMinus+beta betaBar/(2 F));
auxiliary=First[Solve[{D[firstOrder,beta]==0,D[firstOrder,betaBar]==0},{beta,betaBar}]];
check["worldsheet auxiliary elimination",(firstOrder/.auxiliary)
      -(dy barDy+2 Lplus dxPlus barDxPlus+2 Lminus dxMinus barDxMinus-2 F dxMinus barDxPlus)];
(* SM (60)-(62): the clock covectors reconstruct the W-dependent block. *)
X={Cosh[chi/2],-Exp[-sigma] Sinh[chi/2]};
barX={-Exp[sigma] Sinh[chi/2],Cosh[chi/2]};
check["X wedge barX = dx+ wedge dx-",Det[{X,barX}]-1];
check["hair metric block",calW (Outer[Times,X,barX]+Outer[Times,barX,X])
      -calW {{-Exp[sigma] Sinh[chi],Cosh[chi]},{Cosh[chi],-Exp[-sigma] Sinh[chi]}}];

(* SM82 uses lambda matrices, distinct from the SM105 gamma representation. *)
lambda={Sqrt[2] {{0,0},{1,0}},-Sqrt[2] {{0,1},{0,0}},DiagonalMatrix[{1,-1}]};
eta={{0,-1,0},{-1,0,0},{0,0,1}};
check["three-dimensional Clifford algebra",Table[
 lambda[[a]].lambda[[b]]+lambda[[b]].lambda[[a]]-2 eta[[a,b]] IdentityMatrix[2],{a,3},{b,3}]];
(* SM85 -> SM87: k_+=s_{+,i}s_{+,j}; the two Hill equations are inputs.
   This is not a derivation of the full 10d Killing-spinor system or SM89. *)
kPlus=sPlusI[xPlus] sPlusJ[xPlus];
hillRules=Flatten[Table[{Derivative[2][s][xPlus]->2 Lplus[xPlus] s[xPlus]/l^2,
 Derivative[3][s][xPlus]->D[2 Lplus[xPlus] s[xPlus]/l^2,xPlus]},{s,{sPlusI,sPlusJ}}]];
check["Hill bilinear is a local isometry",(kPlus Lplus'[xPlus]
 +2 Lplus[xPlus] D[kPlus,xPlus]-l^2 D[kPlus,{xPlus,3}]/4)/.hillRules];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
