(* ::Title:: *)
(*NRH00 RunAll*)

NRH`$Dir = If[FileExistsQ[FileNameJoin[{DirectoryName[$InputFileName], "NRH01_DFT_Tools.wl"}]], DirectoryName[$InputFileName], NotebookDirectory[]];
NRH`$Files = {
   "NRH02_Letter_Riemannian.wl",
   "NRH03_Letter_NonRiemannian.wl",
   "NRH04_SM_LinearResponse.wl",
   "NRH05_SM_Charges_Action.wl",
   "NRH06_SM_Worldsheet.wl",
   "NRH07_SM_KillingSpinors.wl",
   "NRH08_SM_BoundaryCandidate.wl"};
NRH`$AllResults = {};
NRH`$DeferExit = True;

Scan[Get[FileNameJoin[{NRH`$Dir, #}]] &, NRH`$Files];

NRH`$DeferExit = False;
NRH`GrandSummary[];
