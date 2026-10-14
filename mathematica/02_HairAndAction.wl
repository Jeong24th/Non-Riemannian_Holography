(* Run this file in a fresh kernel: wolframscript -file 02_HairAndAction.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (5), (8)-(9): calPi is the product Pi=L_+ L_-;
   calA,calI,calW denote calligraphic A,I,W (Pi and I are reserved in WL).
   xPlus/xMinus mean x^+/x^-; hold them fixed when differentiating chi. *)
calPi=Lplus[xPlus] Lminus[xMinus];
u[1]=Lplus[xPlus]^(-1/2); u[-1]=Lminus[xMinus]^(-1/2);
calA[1]=D[u[1],{xPlus,2}]/u[1];
calA[-1]=D[u[-1],{xMinus,2}]/u[-1];
calIPrime=4 Sqrt[2] (Sinh[chi]/Sinh[chi/Sqrt[2]]-Sqrt[2]);
calW=(W0[xPlus,xMinus]+W1[xPlus,xMinus] chi/(2 Sqrt[calPi])
 +l^2 Sqrt[calPi]/16 Sum[calA[s] u[s]^2 calI[chi]
 -2 s (D[u[1],xPlus]+s D[u[-1],xMinus])^2 (Exp[s chi]-1-s chi),{s,{-1,1}}]);
(* Input: reduced radial EDFE; the full tensor reduction is outside this file. *)
radialSource=l^2 Sqrt[calPi]/16 Sum[calA[s] u[s]^2 D[calIPrime,chi]
 -2 s (D[u[1],xPlus]+s D[u[-1],xMinus])^2 Exp[s chi],{s,{-1,1}}];
check["hair solves the reduced radial equation",
 (D[calW,{chi,2}]/.calI''[chi]->D[calIPrime,chi])-radialSource];
check["W0 boundary value",((calW/.chi->0)/.calI[0]->0)-W0[xPlus,xMinus]];
check["W1 radial slope",((D[calW,chi]/.chi->0)/.calI'[0]->0)-W1[xPlus,xMinus]/(2 Sqrt[calPi])];
check["integral starts at order chi^3",Limit[calIPrime/chi^2,chi->0]-2/3];

(* SM (2)-(5): action per boundary coordinate volume.
   mu=L_+ L_- on R; mu=L_+ L_-/2 on NR; constant L_pm>0, l>0.
   Input: the on-shell radial flux and no extra interior boundary action. *)
$Assumptions=l>0 && mu>0;
measure=Exp[2 y/l]-mu Exp[-2 y/l];
flux=-2 D[measure,y];
check["radial boundary flux",flux+4/l (Exp[2 y/l]+mu Exp[-2 y/l])];
Sren=((-flux/.y->Y)-8 Sqrt[mu]/l-4 (measure/.y->Y)/l)/(16 Pi G);
check["renormalized saddle action",
 Limit[Sren,Y->Infinity]+8 Sqrt[mu]/(16 Pi G l)];
Print["S_ren per boundary volume = ",-8 Sqrt[mu]/(16 Pi G l)];
(* W does not enter the flux/measure inputs; this is not a new EDFE derivation. *)

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
