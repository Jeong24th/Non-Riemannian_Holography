(* Run this file in a fresh kernel: wolframscript -file 03_TwoPointKernels.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (20), SM (27), (43), (47)-(49), (53)-(55).
   Work in the + sector; interchange + and - for the other sector.
   Inputs: the quadratic-action response operators and SM47 prescription;
   this checks the specified particular kernels, not an interior completion. *)
deltaPlus=xPlus-xPrimePlus; deltaMinus=xMinus-xPrimeMinus;
inverseDerivative=1/(2 Pi I deltaPlus); (* partial_-^(-1) delta^(2) *)
calKplusNR[f_]:=-2 Lplus[xPlus] D[f,xPlus]-Lplus'[xPlus] f;
calKplusR[f_]:=calKplusNR[f]+l^2 D[f,{xPlus,3}]/4;
(* SM27: Hessian factor 1/(4i), applied to the 1/(16 Pi G l) response. *)
AplusNR=Together[calKplusNR[inverseDerivative]/(64 Pi I G l)];
AplusR=Together[calKplusR[inverseDerivative]/(64 Pi I G l)];
MplusNR=Together[(-2 W1[xPlus,xMinus] D[inverseDerivative,xPlus]
 -D[W1[xPlus,xMinus],xPlus] inverseDerivative
 +l^2 Lminus[xMinus] D[inverseDerivative,{xPlus,3}])/(256 Pi I G l)];
check["NR diagonal kernel",AplusNR-(Lplus'[xPlus]/deltaPlus-2 Lplus[xPlus]/deltaPlus^2)/(128 Pi^2 G l)];
check["R central term",AplusR-AplusNR-3 l/(256 Pi^2 G deltaPlus^4)];
check["NR mixed hair kernel",MplusNR-(D[W1[xPlus,xMinus],xPlus]/deltaPlus
 -2 W1[xPlus,xMinus]/deltaPlus^2+6 l^2 Lminus[xMinus]/deltaPlus^4)/(512 Pi^2 G l)];

(* Specified regular R vacuum: differentiate before taking separated points. *)
logKernel=I/(2 Pi deltaPlus deltaMinus);
hairSelf=Together[(l^6/32) D[logKernel,{xPlus,3},{xMinus,3}]/(64 Pi I G l)];
check["R-vacuum hair self-correlator",hairSelf-9 l^5/(1024 Pi^2 G deltaPlus^4 deltaMinus^4)];
Print["A_+^NR = ",AplusNR,"; A_+^R = ",AplusR];
Print["M_+^NR = ",MplusNR,"; R-vacuum hair = ",hairSelf];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
