# Mathematica Reference Execution

This record describes the latest complete reference execution applicable to the
symbolic predicates in this directory.

- Environment: Mathematica 13.2 (Windows), `wolframscript`
- Entry point: `NRH00_RunAll.wl`
- Command: `wolframscript -file NRH00_RunAll.wl`
- Result: 243/243 exact checks passed
  (NRH02 30, NRH03 32, NRH04 49, NRH05 38, NRH06 28, NRH07 52, NRH08 14)
- Exit behavior: a failed check produces a nonzero command-line exit status
- Typical runtime: a few minutes (190 s wall time for the recorded execution); the longest single file is
  `NRH04_SM_LinearResponse.wl` (the second-order expansion of the Gamma^2 density with
  four fluctuation channels and the exact variation identity on the BTZ background)

The `.nb` notebooks are generated cell-for-cell from the `.wl` sources by a kernel-side
converter; evaluating a notebook runs exactly the same checks.

This repository does not bundle a Wolfram runtime or license and does not contain an
automated hosted execution.  Independent reproduction therefore requires Mathematica
or Wolfram Engine and a fresh run of the command above.
