(* ::Title:: *)
(*NRH00 RunAll*)

NRH`$Dir = If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]];
NRH`$Files = {
   "NRH02_Letter.wl",
   "NRH03_SM1_Action.wl",
   "NRH04_SM2_RadialBranches.wl",
   "NRH05_SM3_Linearized.wl",
   "NRH06_SM3_Renormalization.wl",
   "NRH07_SM3_Charges.wl",
   "NRH08_SM4_Worldsheet.wl",
   "NRH09_SM5_UpliftKilling.wl",
   "NRH10_SM6_BoundaryCandidate.wl"};
NRH`$AllResults = {};
NRH`$DeferExit = True;

Scan[Get[FileNameJoin[{NRH`$Dir, #}]] &, NRH`$Files];

NRH`$DeferExit = False;
NRH`GrandSummary[];
