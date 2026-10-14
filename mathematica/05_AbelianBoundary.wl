(* Run this file in a fresh kernel: wolframscript -file 05_AbelianBoundary.wl
   Zero residual means the stated identity holds. *)
ClearAll["Global`*"];
passed=0; total=0;
check[name_, residual_] := (total++; If[
  And @@ (TrueQ[Simplify[#]==0]& /@ Flatten[{residual}]),
  passed++; Print["[PASS] ",name], Print["[FAIL] ",name]]);

(* SM (94)-(108): unrestricted B; no field chirality/self-duality constraint.
   Constant DFT vielbeins of Letter (4b), gamma matrices of SM105. *)
gammaPlus=Sqrt[2] {{0,1},{0,0}};
gammaMinus=Sqrt[2] {{0,0},{-1,0}};
barGammaMinus=-gammaMinus;
chirality=DiagonalMatrix[{1,-1}]; chargeConjugation={{0,1},{1,0}};
bispinor={{c[x,y],bplus[x,y]},{bminus[x,y],phi[x,y]}};
rr[z_,sign_:1]:=(gammaPlus.D[z,x]+sign chirality.D[z,y].barGammaMinus)/Sqrt[2];
check["D_plus is nilpotent (also RR gauge invariance)",rr[rr[bispinor]]];
check["D_minus is nilpotent",rr[rr[bispinor,-1],-1]];
fieldStrength=rr[bispinor];
barStrength=chargeConjugation.Transpose[fieldStrength].chargeConjugation;
check["bosonic action reduces to one scalar",
      -Tr[fieldStrength.barStrength]/2-D[phi[x,y],x] D[phi[x,y],y]];
check["unbarred fermion kinetic term selects psi-minus",
      chargeConjugation.gammaPlus/Sqrt[2]-DiagonalMatrix[{0,1}]];
check["barred fermion kinetic term selects psi-plus",
      chargeConjugation.barGammaMinus/Sqrt[2]-DiagonalMatrix[{1,0}]];

(* SM103: Cartan identity for arbitrary chiral v^plus(x), v^minus(y). *)
insertion[z_]:=-vp[x] gammaMinus.z/Sqrt[2]+vm[y] chirality.z.gammaPlus/Sqrt[2];
kosmann[z_]:=vp[x] D[z,x]+vm[y] D[z,y]+(vp'[x]+vm'[y]) z/2
             +vp'[x] chirality.z/2+vm'[y] z.chirality/2;
check["Cartan anticommutator equals Kosmann",rr[insertion[bispinor]]+insertion[rr[bispinor]]-kosmann[bispinor]];
Print["F = ",MatrixForm[fieldStrength]];
(* Null components have no propagator fixed by this action. Full Grassmann
   SUSY/closure checks remain available in the linked detailed release. *)

Print[passed,"/",total," checks passed."];
If[$FrontEnd===Null, Exit[If[passed==total,0,1]]];
