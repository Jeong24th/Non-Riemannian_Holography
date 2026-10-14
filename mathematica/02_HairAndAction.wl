(* Run this file in a fresh kernel: wolframscript -file 02_HairAndAction.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (8)-(9): integrate the stated radial hair equation.
   A=u_+ u_+''+u_- u_-'', P=u_+'+u_-', M=u_+'-u_-'.
   These source combinations are constant with respect to chi. *)
integralPrime=4 Sqrt[2] (Sinh[chi]/Sinh[chi/Sqrt[2]]-Sqrt[2]);
hair=W0+W1 chi/(2 Sqrt[product])+l^2 Sqrt[product]/16 (
  A radialIntegral[chi]-2 P^2 (Exp[chi]-1-chi)
                          +2 M^2 (Exp[-chi]-1+chi));
radialSource=l^2 Sqrt[product]/16 (
  A D[integralPrime,chi]-2 P^2 Exp[chi]+2 M^2 Exp[-chi]);
check["hair solves the reduced radial equation",
 (D[hair,{chi,2}]/.radialIntegral''[chi]->D[integralPrime,chi])-radialSource];
check["W0 boundary value",((hair/.chi->0)/.radialIntegral[0]->0)-W0];
check["W1 radial slope",((D[hair,chi]/.chi->0)/.radialIntegral'[0]->0)-W1/(2 Sqrt[product])];
check["integral starts at order chi^3",Limit[integralPrime/chi^2,chi->0]-2/3];

(* SM (2)-(5): action per boundary coordinate volume.
   mu=Lp Lm on R; mu=Lp Lm/2 on NR. u=exp(2Y/l).
   Input: the on-shell radial flux and no extra interior boundary action. *)
measure=u-mu/u;
flux=-2 (2 u/l) D[measure,u];
check["radial boundary flux",flux+4/l (u+mu/u)];
bare=(-flux-8 Sqrt[mu]/l)/(16 Pi GN);
counterterm=-4 measure/(16 Pi GN l);
renormalized=Together[bare+counterterm];
check["renormalized saddle action",
 Limit[renormalized,u->Infinity]+8 Sqrt[mu]/(16 Pi GN l)];
Print["S_ren per boundary volume = ",-8 Sqrt[mu]/(16 Pi GN l)];
(* W does not enter the flux/measure inputs; this is not a new EDFE derivation. *)

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
