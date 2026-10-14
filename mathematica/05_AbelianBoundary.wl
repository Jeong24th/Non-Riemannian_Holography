(* Run this file in a fresh kernel: wolframscript -file 05_AbelianBoundary.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* SM (94)-(108): unrestricted calB; no field chirality/self-duality constraint.
   Constant DFT vielbeins of Letter (4b), gamma matrices of SM105. *)
gammaOplus=Sqrt[2] {{0,1},{0,0}}; gammaOminus=Sqrt[2] {{0,0},{-1,0}};
barGammaOplus=gammaOplus; barGammaOminus=-gammaOminus;
gamma3=DiagonalMatrix[{1,-1}]; barGamma3=-gamma3;
chargeC={{0,1},{1,0}}; barChargeC=chargeC; (* C,bar C; C is reserved in WL. *)
calB={{c[xPlus,xMinus],bPlus[xPlus,xMinus]},{bMinus[xPlus,xMinus],phi[xPlus,xMinus]}};
Dplus[z_]:=(gammaOplus.D[z,xPlus]+gamma3.D[z,xMinus].barGammaOminus)/Sqrt[2];
Dminus[z_]:=(gammaOplus.D[z,xPlus]-gamma3.D[z,xMinus].barGammaOminus)/Sqrt[2];
check["D_+ nilpotency and RR gauge invariance",Dplus[Dplus[calB]]];
check["D_- nilpotency",Dminus[Dminus[calB]]];
calF=Dplus[calB]; barCalF=Inverse[barChargeC].Transpose[calF].chargeC;
check["bosonic action reduces to one scalar",
      -Tr[calF.barCalF]/2-D[phi[xPlus,xMinus],xPlus] D[phi[xPlus,xMinus],xMinus]];
check["unbarred fermion kinetic term selects psi-minus",
      chargeC.gammaOplus/Sqrt[2]-DiagonalMatrix[{0,1}]];
check["barred fermion kinetic term selects psi-plus",
      barChargeC.barGammaOminus/Sqrt[2]-DiagonalMatrix[{1,0}]];

(* SM97,103: insertion=iota_v; kosmann=frak L_v at weight 1/2;
   vPlus(xPlus),vMinus(xMinus) are the arbitrary chiral vector components. *)
insertion[z_]:=(-vPlus[xPlus] gammaOminus.z+vMinus[xMinus] gamma3.z.barGammaOplus)/Sqrt[2];
kosmann[z_]:=(vPlus[xPlus] D[z,xPlus]+vMinus[xMinus] D[z,xMinus]
 +(vPlus'[xPlus]+vMinus'[xMinus]) z/2
 +vPlus'[xPlus] gamma3.z/2-vMinus'[xMinus] z.barGamma3/2);
check["Cartan anticommutator equals Kosmann",Dplus[insertion[calB]]+insertion[Dplus[calB]]-kosmann[calB]];
Print["calF = ",MatrixForm[calF]];
(* Null components have no propagator fixed by this action. Full Grassmann
   SUSY/closure checks remain available in the linked detailed release. *)

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
