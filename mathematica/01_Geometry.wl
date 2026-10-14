(* Run this file in a fresh kernel: wolframscript -file 01_Geometry.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (1)-(2), (4a), (6): generalized metrics and shared boundary.
   This short check does not rederive the full tensor EDFE. *)
$Assumptions=l>0 && Element[{y,Lplus,Lminus},Reals];
calJ=ArrayFlatten[{{0,IdentityMatrix[3]},{IdentityMatrix[3],0}}];
(* Lplus,Lminus denote L_+(x^+),L_-(x^-); phi_0=0. *)
F=Exp[2 y/l]+Lplus Lminus Exp[-2 y/l];
g={{2 Lplus,-F,0},{-F,2 Lminus,0},{0,0,1}};
B=F {{0,-1,0},{1,0,0},{0,0,0}};
calHR=Map[Together,ArrayFlatten[{{Inverse[g],-Inverse[g].B},
                              {B.Inverse[g],g-B.Inverse[g].B}}],{2}];
check["R: H J^(-1) H = J",calHR.Inverse[calJ].calHR-calJ];
check["R: dilaton measure squared = -det(g)",
      (Exp[2 y/l]-Lplus Lminus Exp[-2 y/l])^2+Det[g]];

(* chi, sigma and calW are unrestricted in this algebraic identity. *)
calHNR={{0,0,0,Cosh[chi],-Exp[-sigma] Sinh[chi],0},
 {0,0,0,Exp[sigma] Sinh[chi],-Cosh[chi],0},
 {0,0,1,0,0,0},
 {Cosh[chi],Exp[sigma] Sinh[chi],0,-calW Exp[sigma] Sinh[chi],calW Cosh[chi],0},
 {-Exp[-sigma] Sinh[chi],-Cosh[chi],0,calW Cosh[chi],-calW Exp[-sigma] Sinh[chi],0},
 {0,0,0,0,0,1}};
check["NR: H J^(-1) H = J",calHNR.Inverse[calJ].calHNR-calJ];
check["NR: inverse-metric block has rank one",calHNR[[1;;3,1;;3]]-DiagonalMatrix[{0,0,1}]];
calHInfinity={{0,0,0,1,0,0},{0,0,0,0,-1,0},{0,0,1,0,0,0},
           {1,0,0,0,0,0},{0,-1,0,0,0,0},{0,0,0,0,0,1}};
check["R boundary",Map[Limit[#,y->Infinity]&,calHR,{2}]-calHInfinity];
check["NR boundary at W0=0",(calHNR/.{chi->0,calW->0})-calHInfinity];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
