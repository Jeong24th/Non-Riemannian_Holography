(* Run this file in a fresh kernel: wolframscript -file 03_TwoPointKernels.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* Letter (20), SM (43), (47)-(49).
   Inputs: source-normalized response operators and the stated Lorentzian
   inverse-derivative prescription. No homogeneous/interior completion. *)
inverseDerivative=1/(2 Pi I delta);
normalization=1/(16 Pi GN l);
wardNR[f_]:=2 LL D[f,delta]+Lprime f;
wardR[f_]:=wardNR[f]-l^2 D[f,{delta,3}]/4;
hairResponse[f_]:=2 WW D[f,delta]+Wprime f-l^2 Lopposite D[f,{delta,3}];
diagonalNR=Together[-normalization wardNR[inverseDerivative]/(4 I)];
diagonalR=Together[-normalization wardR[inverseDerivative]/(4 I)];
mixedNR=Together[-normalization hairResponse[inverseDerivative]/(16 I)];
check["NR diagonal kernel",diagonalNR-(Lprime/delta-2 LL/delta^2)/(128 Pi^2 GN l)];
check["R central term",diagonalR-diagonalNR-3 l/(256 Pi^2 GN delta^4)];
check["NR mixed hair kernel",mixedNR-(Wprime/delta-2 WW/delta^2+6 l^2 Lopposite/delta^4)/(512 Pi^2 GN l)];

(* Specified regular R vacuum: differentiate before taking separated points. *)
logKernel=I/(2 Pi deltaPlus deltaMinus);
hairSelf=Together[(l^6/32) D[logKernel,{deltaPlus,3},{deltaMinus,3}]/(64 Pi I GN l)];
check["R-vacuum hair self-correlator",hairSelf-9 l^5/(1024 Pi^2 GN deltaPlus^4 deltaMinus^4)];
Print["NR diagonal = ",diagonalNR,"; R diagonal = ",diagonalR];
Print["NR mixed = ",mixedNR,"; R-vacuum hair = ",hairSelf];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
