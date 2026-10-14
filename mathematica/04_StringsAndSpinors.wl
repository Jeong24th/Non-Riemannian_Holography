(* Run this file in a fresh kernel: wolframscript -file 04_StringsAndSpinors.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* SM (58): eliminate beta and betaBar by their own field equations.
   plus/minus = partial x^plus/minus; barred = bar-partial x^plus/minus. *)
firstOrder=radial+2 Lp plus barredPlus+2 Lm minus barredMinus
           -beta barredMinus-betaBar plus+beta betaBar/(2 F);
auxiliary=First[Solve[{D[firstOrder,beta]==0,D[firstOrder,betaBar]==0},{beta,betaBar}]];
check["worldsheet auxiliary elimination",(firstOrder/.auxiliary)
      -(radial+2 Lp plus barredPlus+2 Lm minus barredMinus-2 F plus barredMinus)];
(* SM (60)-(62): the clock covectors reconstruct the W-dependent block. *)
clock={Cosh[chi/2],-Exp[-sigma] Sinh[chi/2]};
barClock={-Exp[sigma] Sinh[chi/2],Cosh[chi/2]};
check["clock wedge is dx+ wedge dx-",Det[{clock,barClock}]-1];
check["hair metric block",W (Outer[Times,clock,barClock]+Outer[Times,barClock,clock])
      -W {{-Exp[sigma] Sinh[chi],Cosh[chi]},{Cosh[chi],-Exp[-sigma] Sinh[chi]}}];

(* SM (89): unprimed reduced NR-vacuum spinor channels.
   The projected operator is an input; this is not a full 10d SDFT derivation. *)
dirac[e_]:=DiagonalMatrix[{1,-1}].e/(Sqrt[2] l)+{{0,0},{1,0}}.D[e,x];
positive={Sqrt[2] f[x],l f'[x]}; negative={0,g[x]};
check["arbitrary chiral f channel",dirac[positive]-positive/(Sqrt[2] l)];
check["arbitrary chiral g channel",dirac[negative]+negative/(Sqrt[2] l)];

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
