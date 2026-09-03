(* ::Package:: *)

(* ::Title:: *)
(*NRH00 — Run every verification file*)


(* ::Text:: *)
(*Runs the complete Mathematica verification suite for "Long Strings and Non-Riemannian*)
(*Hair" in paper order and prints a grand PASS/FAIL summary.  Each section file clears*)
(*the Global` context and reloads the shared toolbox, so they can equally be run alone.*)
(**)
(*Usage (command line):    wolframscript -file NRH00_RunAll.wl*)
(*Usage (notebook):        open this file in Mathematica and evaluate it.*)
(*The whole suite completes in a few minutes on a laptop (243 checks; see REFERENCE_RUN.md).*)


NRH`$Dir = If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]];
NRH`$Files = {
   "NRH02_Letter_Riemannian.wl",
   "NRH03_Letter_NonRiemannian.wl",
   "NRH04_SM_LinearResponse.wl",
   "NRH05_SM_Charges_Action.wl",
   "NRH06_SM_Worldsheet.wl",
   "NRH07_SM_KillingSpinors.wl",
   "NRH08_SM_BoundaryCandidate.wl"};
NRH`$AllResults = {};
NRH`$DeferExit = True;   (* report everything; exit code is decided by the grand summary *)

Scan[Get[FileNameJoin[{NRH`$Dir, #}]] &, NRH`$Files];

NRH`$DeferExit = False;
NRH`GrandSummary[];
