# September 2026 SM revision

Current manuscript: `57960AAD42F2A9150C0958D1C253F08B0939B76C0E646EE250A7D824D6319B68` (2026-09-30).
Previous documentation pin: `EFB2D75953E9306C92165B532E855C590CA1FF8929FC61BF8EC8531BC522CD94` (2026-09-19).

The current SM3 replaces the longer derivation with six compact subsections. The public
calculations retain useful expanded checks, including equations no longer displayed in
the manuscript. Their algebra is unchanged; comments and assertion descriptions are
updated, and affected notebooks contain the same inputs as their `.wl` counterparts.
There has been no new Wolfram execution. Current coverage is audited in
[EQUATION_LEDGER.md](EQUATION_LEDGER.md); the historical run remains separately dated.

## Current SM3 organization

| Current subsection | Existing files and scope |
|---|---|
| 3.1 Constrained Field Variations | NRH05 frame/connection checks; full universal Box identities are cited, not directly tested. |
| 3.2 Renormalized On-Shell Action | NRH06 radial momenta and ordered asymptotic bilinear; the new curved-index Hessian is not directly compared. |
| 3.3 Linearized EDFE with an Interior Boundary Condition | NRH05 direct component linearization; K3 regularity and the interior determination of H_R are new coverage gaps. |
| 3.4 Riemannian Branch | NRH05/06 earlier differential hierarchy and coefficient ingredients; compact-expression bridges are not public checks. |
| 3.5 Non-Riemannian Branch | NRH05/06 earlier hierarchy and particular kernels, with the interior response left free. |
| 3.6 Covariant Charges and Asymptotic Algebras | NRH07 boundary charge and algebra checks. |

## Notation and mathematical scope

- The radial profile is now calligraphic I(chi); LaTeX labels `MainGprofile` and
  `NRGprofile` are unchanged. Existing GG/Gp Wolfram symbols are retained (see README);
  bare `I` is the protected imaginary unit and is not a radial function.
- Box is the full DFT tensor operator with connection and curvature terms. Setting
  the projected Ricci tensor to zero does not discard the curvature terms inside Box.
- The explicit Einstein-tensor decomposition is now included before the linearized
  equations. The new off-shell Codazzi equation is not proved by exact-saddle checks.
- Code r is the third, type-changing source h_(ominus baroplus)^(0); c is the fourth
  W0-source channel h_(oplus barominus)^(0). Code cs on R equals -4 f0. Rp/Rm are the
  full order-u stress responses and H is the unshifted free hair response.
- The checked NR mixed-kernel coefficient remains 1/(512 pi^2 G l). No factor-two
  scalar renormalization or change of the quarter-Hessian convention is made.
- Assertion wording is narrowed where a predicate checks only a sufficient log-free
  sector, an undressed longitudinal contraction, a historical complex rank, or count arithmetic.
  These wording corrections do not enlarge the mathematical predicates.
- Supersymmetry counts distinguish real polarizations, arbitrary chiral function
  families and nonzero integrable supercharges. The latter are not constructed here.

## Label migration from the previous documentation

Numbers in the first column belong to the 2026-09-19 source. “No current display”
means the label was removed or its content absorbed/replaced; it does not delete an
existing regression calculation or certify a new replacement. Unlabeled displays
are omitted from this migration table; the current map includes them.

| Previous no. | Stable label | Current no. |
|---:|---|---:|
| 1 | `Rfields` | 1 |
| 2 | `RDFTfields` | 2 |
| 3 | `Rboundary` | 3 |
| 4a | `boundaryH` | 4a |
| 4b | `Rboundaryframe` | 4b |
| 5 | `NRvariables` | 5 |
| 6 | `NRHcompact` | 6 |
| 7 | `NRdilaton` | 7 |
| 8 | `NRWgeneral` | 8 |
| 9 | `MainGprofile` | 9 |
| 10 | `RGKPW` | 10 |
| 11 | `Rrenvariation` | 11 |
| 12 | `Mainmomenta` | 12 |
| 13a | `RKdef` | 14a |
| 13b | `RDFTconservation` | 14b |
| 14 | `Rcontinuity` | No current display |
| 15 | `Rkilling` | 15 |
| 16 | `RVirasoro` | 16 |
| 17 | `NRasympt` | 17 |
| 18 | `Mainonepoints` | 18 |
| 19 | `Mainframevariation` | No current display |
| 20 | `Mainfivebyfive` | 19 |
| 21 | `Rcorrelators` | 20 |
| 22 | `Mainworldsheet` | 21 |
| SM 1 | `SMgamma2` | SM 1 |
| SM 2 | `SMgamma2flux` | SM 2 |
| SM 3 | `SMmudefinition` | SM 3 |
| SM 4 | `SMgamma2cutoff` | SM 4 |
| SM 5 | `SMgamma2value` | SM 5 |
| SM 7 | `NRhill` | SM 7 |
| SM 8 | `NRradialchange` | SM 8 |
| SM 9 | `NRradialoperator` | SM 9 |
| SM 10 | `NRchiODE` | SM 10 |
| SM 11 | `NRsource` | SM 11 |
| SM 12 | `NRg` | SM 12 |
| SM 13 | `NRGprofile` | SM 13 |
| SM 16 | `SMBtransform` | SM 16 |
| SM 17 | `SMWshift` | SM 17 |
| SM 18 | `SMprojectors` | No current display |
| SM 19 | `SMmixedfluctuation` | No current display |
| SM 20 | `SMcosetreconstruction` | SM 18 |
| SM 21 | `SMframevariation` | SM 19 |
| SM 22 | `SMresponsedef` | No current display |
| SM 23 | `SMbackgroundconnection` | SM 21 |
| SM 24 | `SMconnectionvariation` | No current display |
| SM 25 | `SMsixprojectors` | No current display |
| SM 26 | `SMflatmetrics` | No current display |
| SM 27 | `SMinfinityvielbein` | SM 20 |
| SM 28 | `SMFG` | No current display |
| SM 29 | `SMsources` | No current display |
| SM 30 | `SMW0completion` | No current display |
| SM 31 | `SMW0variation` | No current display |
| SM 32 | `SMbackgroundexpansion` | No current display |
| SM 33 | `SMexactboxEDFE` | SM 39 |
| SM 34 | `SMboxdefinition` | SM 25 |
| SM 35 | `SMmixedboxdefinition` | SM 26 |
| SM 37 | `SMboxcurvature` | SM 28 |
| SM 38 | `SMexactRdata` | SM 40 |
| SM 39 | `SMexactNRdata` | No current display |
| SM 40 | `SMexactcomponentrecipe` | No current display |
| SM 42 | `SMexactfieldorder` | No current display |
| SM 43 | `SMexactnormalform` | No current display |
| SM 44 | `SMexactcoefficientpolynomial` | No current display |
| SM 45 | `SMexactTaylor` | No current display |
| SM 46 | `SMexactTaylorSolution` | No current display |
| SM 47 | `SMexactCauchyconstraints` | No current display |
| SM 48 | `SMexactconstraintpropagation` | No current display |
| SM 49 | `SMoperatorseries` | No current display |
| SM 50 | `SMradialansatz` | SM 49 |
| SM 51 | `SMNRradialintegration` | No current display |
| SM 52 | `SMtotalorderhierarchy` | No current display |
| SM 53 | `SMcoupledorders` | No current display |
| SM 54 | `SMmixingexample` | No current display |
| SM 55 | `SMEzerocomponents` | No current display |
| SM 56 | `SMEoneNRcomponents` | No current display |
| SM 57 | `SMRoperatorrule` | No current display |
| SM 58 | `SMDeltaEonecomponents` | No current display |
| SM 59 | `SMleadingdilaton` | No current display |
| SM 60 | `SMleadingstress` | No current display |
| SM 61 | `SMNRsol` | No current display |
| SM 62 | `SMlogexample` | No current display |
| SM 63 | `SMRh2L_omopb` | No current display |
| SM 64 | `SMRdd2L` | No current display |
| SM 65 | `SMRh2L_opopb` | No current display |
| SM 66 | `SMRh2L_omomb` | No current display |
| SM 67 | `SMRdd2` | No current display |
| SM 68 | `SMRh2L_opomb` | No current display |
| SM 69 | `SMNRh2L_omopb` | No current display |
| SM 70 | `SMNRdd2L` | No current display |
| SM 71 | `SMNRdd2` | No current display |
| SM 72 | `SMNRh2L_opopb` | No current display |
| SM 73 | `SMNRh2L_omomb` | No current display |
| SM 74 | `SMNRh2L_opomb` | No current display |
| SM 75 | `SMtypeconstraint` | No current display |
| SM 76 | `SMRconstraint0` | No current display |
| SM 77 | `SMRconstraint1` | No current display |
| SM 78 | `SMNRconstraint0` | No current display |
| SM 79 | `SMNRconstraint1` | No current display |
| SM 80 | `SMgeneralradialsolution` | No current display |
| SM 81 | `SMlogfreeconditions` | No current display |
| SM 82 | `SMfullscalarvariation` | SM 22 |
| SM 83 | `SMfullcurvaturevariation` | SM 23 |
| SM 84 | `SMfullOmega` | SM 24 |
| SM 85 | `GammaDFT` | SM 29 |
| SM 86 | `variation` | SM 30 |
| SM 87a | `SMBdefinition` | SM 31a |
| SM 87b | `SMAdefinition` | SM 31b |
| SM 88 | `SMrenvariation` | SM 33 |
| SM 89 | `SMvolumecounterterm` | SM 32 |
| SM 90a | `SMmetricradialdictionary` | No current display |
| SM 90b | `SMscalarradialdictionary` | No current display |
| SM 91 | `SMsecondvariation` | SM 35 |
| SM 92 | `SMmomentumvariation` | No current display |
| SM 93 | `SMBmomentumvariation` | No current display |
| SM 94 | `SMcurvedmomentumdefinition` | No current display |
| SM 95 | `SMprojectedconnectionmomentum` | No current display |
| SM 96 | `SMbackgroundmomenta` | No current display |
| SM 97 | `SMsamechiralitymomenta` | No current display |
| SM 98 | `SMonept` | No current display |
| SM 99 | `SMNRonept` | No current display |
| SM 100 | `SMstatefluctuations` | SM 59 |
| SM 101 | `SMderivativectNR` | SM 36 |
| SM 102 | `SMderivativectR` | SM 37 |
| SM 103 | `SMctlocality` | No current display |
| SM 104 | `SMctadjoint` | No current display |
| SM 105 | `SMboundarycurvaturedefinitions` | No current display |
| SM 106 | `SMboundarycurvaturecandidate` | No current display |
| SM 107 | `SMboundarycurvatureresidual` | No current display |
| SM 108 | `SMdirectsources` | No current display |
| SM 109 | `SMdirectbackground` | No current display |
| SM 110 | `SMframemomentumcheck` | No current display |
| SM 111 | `SMconnectionmomentumcheck` | No current display |
| SM 112 | `SMlinearizedmomenta` | No current display |
| SM 113 | `SMdirectcutoffbilinear` | No current display |
| SM 114 | `SMfinitepartexplicit` | No current display |
| SM 115 | `SMfullfinitehessian` | No current display |
| SM 116 | `SMcontactconstraintreduction` | No current display |
| SM 117 | `SMdirectstressrows` | No current display |
| SM 118 | `SMdirectremainingrows` | No current display |
| SM 119 | `SMdirectkernel` | No current display |
| SM 120 | `SMRscalarcontactrow` | No current display |
| SM 121 | `SMNRscalarcontactrow` | No current display |
| SM 122 | `SMLorentziandeltaconvention` | SM 46 |
| SM 123 | `SMfullhessianclosure` | No current display |
| SM 124 | `SMlinearWeylanomaly` | No current display |
| SM 125 | `SMstresssourceexample` | No current display |
| SM 126 | `SMgeneralWardoperators` | No current display |
| SM 127 | `SMstressparticular` | No current display |
| SM 128 | `SMindependentsourcegenerator` | No current display |
| SM 129 | `SMhairparticularsolution` | No current display |
| SM 130 | `SMgeneralparticularkernels` | No current display |
| SM 131 | `SMPBHkernel` | No current display |
| SM 132 | `SMgeneralpositionkernels` | SM 60 |
| SM 133 | `SMRtwopt` | SM 47 |
| SM 134 | `SMNRhairhessian` | SM 61 |
| SM 135 | `SMRhairuniformization` | No current display |
| SM 136 | `SMRhairconditional` | No current display |
| SM 137 | `SMRhairconstantL` | No current display |
| SM 138 | `SMonepointequivalence` | No current display |
| SM 139 | `SMCPSform` | SM 62 |
| SM 140 | `SMRpotentialcomponents` | SM 63 |
| SM 141 | `SMCPSresult` | SM 64 |
| SM 142 | `SMNRchargefalloffs` | No current display |
| SM 143 | `SMNRchargecancellation` | No current display |
| SM 144 | `SMNRcharge` | SM 65 |
| SM 145 | `SMCPSalgebra` | No current display |
| SM 146 | `SMdyg` | SM 66 |
| SM 147 | `SMphysicalsectionA` | SM 67 |
| SM 148 | `SMRfirstorder` | SM 68 |
| SM 149 | `SMceff` | SM 69 |
| SM 151 | `SMlongstringE` | SM 71 |
| SM 152 | `SNCtau` | SM 72 |
| SM 153 | `SMdygconstraints` | SM 73 |
| SM 154 | `SMdygGO` | SM 74 |
| SM 155 | `SMGO` | SM 75 |
| SM 156 | `SMvertex` | SM 76 |
| SM 157 | `SMWisB` | SM 77 |
| SM 158 | `SMdeltaL` | SM 78 |
| SM 159 | `SMFTsector` | SM 79 |
| SM 160 | `SMWZWlevel` | SM 80 |
| SM 161 | `SMFTweight` | SM 81 |
| SM 162 | `SMBRSTradial` | SM 82 |
| SM 163 | `SMBRSTmomentum` | SM 83 |
| SM 164 | `SMBRSTvertex` | SM 84 |
| SM 165 | `SMBRSTcentral` | SM 85 |
| SM 166 | `SMWgaugeobstruction` | SM 86 |
| SM 167 | `SMBRSTfusion` | SM 87 |
| SM 168 | `SMupliftblocks` | SM 88 |
| SM 169 | `SMDFTKilling` | SM 89 |
| SM 170 | `SMtypeIIKS` | SM 90 |
| SM 171 | `SMsusyclosure` | SM 91 |
| SM 172 | `SMsemicov` | SM 92 |
| SM 173 | `SMuplift` | SM 93 |
| SM 174 | `SMRDFTKilling` | SM 94 |
| SM 175 | `SMRlocaliso` | SM 95 |
| SM 178 | `SMRcomponentHill` | SM 98 |
| SM 179 | `SMRlocalKS` | SM 99 |
| SM 181 | `SMvielbein` | SM 101 |
| SM 182 | `SMvielbeincheck` | SM 102 |
| SM 183 | `SMNRlocalstabilizer` | SM 103 |
| SM 184 | `SMexactiso` | SM 104 |
| SM 185 | `SMweighteddilaton` | SM 105 |
| SM 186 | `SMkillingspinor` | SM 106 |
| SM 187 | `SMcomplexblocks` | SM 107 |
| SM 189 | `SMgammaten` | SM 109 |
| SM 190 | `SMgammaMajorana` | SM 110 |
| SM 191 | `SMgammabarred` | SM 111 |
| SM 192 | `SMreducedDirac` | SM 112 |
| SM 193 | `SMinternalprojectors` | SM 113 |
| SM 195 | `SMspinorcountchain` | SM 115 |
| SM 196 | `SMhairyKS` | SM 116 |
| SM 199 | `SMcandidateD` | SM 119 |
| SM 200 | `SMcandidateaction` | SM 120 |
| SM 201 | `SMcandidategauge` | SM 121 |
| SM 202 | `SMcandidateWitt` | SM 122 |
| SM 203 | `SMcandidateWittL` | SM 123 |
| SM 204 | `SMcandidatefermionic` | SM 124 |
| SM 205 | `SMcandidatefermionicL` | SM 125 |
| SM 206 | `SMcandidateextra` | SM 126 |
| SM 207 | `SMcandidateextraL` | SM 127 |
| SM 208 | `SMcandidatetrivial` | SM 128 |
| SM 209 | `SMcandidatetrivialL` | SM 129 |

## Newly labeled displays

| Current no. | Label |
|---:|---|
| 13 | `Codazzi` |
| SM 34 | `SMradialmetricvariation` |
| SM 38 | `SMEinsteintensor` |
| SM 41 | `SMRvacuumdilatoncoefficient` |
| SM 42 | `SMRvacuumsolution` |
| SM 43 | `SMRvacuumregularscalar` |
| SM 44 | `SMRvacuumresponse` |
| SM 45 | `SMRvacuumhessian` |
| SM 48 | `SMRvacuumhaircorrelator` |
| SM 50 | `SMRstresspluscoefficients` |
| SM 51 | `SMRstressminuscoefficients` |
| SM 52 | `SMRtypecoefficients` |
| SM 53 | `SMRhaircoefficients` |
| SM 54 | `SMRdilatoncoefficients` |
| SM 55 | `SMRintegratedstress` |
| SM 56 | `SMNRcompactsolution` |
| SM 57 | `SMNRintegratedstress` |
| SM 58 | `SMNRcompactU` |

These labels identify the current statements; their independent verification status
is recorded in the ledger. New labels do not imply new machine checks.
