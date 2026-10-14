(* Run this file in a fresh kernel: wolframscript -file 01_Geometry.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (2), (4a), (6): generalized metrics and their shared boundary.
   This short check does not rederive the full tensor EDFE. *)
j=ArrayFlatten[{{0,IdentityMatrix[3]},{IdentityMatrix[3],0}}];
(* u=exp(2y/l); Lp,Lm are held fixed in this radial calculation. *)
g={{2 Lp,-u-Lp Lm/u,0},{-u-Lp Lm/u,2 Lm,0},{0,0,1}};
b=(u+Lp Lm/u) {{0,-1,0},{1,0,0},{0,0,0}};
hR=Map[Together,ArrayFlatten[{{Inverse[g],-Inverse[g].b},
                              {b.Inverse[g],g-b.Inverse[g].b}}],{2}];
check["R: H J H = J",hR.j.hR-j];
check["R: dilaton measure squared = -det(g)",(u-Lp Lm/u)^2+Det[g]];

(* chi, sigma and W are unrestricted in this algebraic identity. *)
hNR={{0,0,0,Cosh[chi],-Exp[-sigma] Sinh[chi],0},
 {0,0,0,Exp[sigma] Sinh[chi],-Cosh[chi],0},
 {0,0,1,0,0,0},
 {Cosh[chi],Exp[sigma] Sinh[chi],0,-W Exp[sigma] Sinh[chi],W Cosh[chi],0},
 {-Exp[-sigma] Sinh[chi],-Cosh[chi],0,W Cosh[chi],-W Exp[-sigma] Sinh[chi],0},
 {0,0,0,0,0,1}};
check["NR: H J H = J",hNR.j.hNR-j];
check["NR: inverse-metric block has rank one",hNR[[1;;3,1;;3]]-DiagonalMatrix[{0,0,1}]];
hInfinity={{0,0,0,1,0,0},{0,0,0,0,-1,0},{0,0,1,0,0,0},
           {1,0,0,0,0,0},{0,-1,0,0,0,0},{0,0,0,0,0,1}};
check["R boundary",Map[Limit[#,u->Infinity]&,hR,{2}]-hInfinity];
check["NR boundary at W0=0",(hNR/.{chi->0,W->0})-hInfinity];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
