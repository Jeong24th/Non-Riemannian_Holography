# Mathematica Reference Execution

This record describes the latest complete reference execution applicable to the
symbolic predicates in this directory.

**Pending re-execution (2026-09-04).**  The two ten-dimensional curvature checks in
`NRH07_SM_KillingSpinors.wl` were changed from exact rational-parameter probes to fully
symbolic identities (arbitrary chiral `L_±(x^±)`; one-sided hairy family with arbitrary
`W_0`, `W_1`), and the stale comment in `NRH01_DFT_Tools.wl` about "numerical spot
checks" was removed.  The number of checks is unchanged (five in that section), but the
execution below predates this change and must be repeated, and the `.nb` notebooks must
be regenerated from the `.wl` sources, before the record is current.

- Environment: Mathematica 13.2 (Windows), `wolframscript`
- Entry point: `NRH00_RunAll.wl`
- Command: `wolframscript -file NRH00_RunAll.wl`
- Result: 247/247 exact checks passed
  (NRH02 30, NRH03 32, NRH04 49, NRH05 38, NRH06 32, NRH07 52, NRH08 14)
- Exit behavior: a failed check produces a nonzero command-line exit status
- Typical runtime: a few minutes (169 s wall time for the recorded execution); the longest single file is
  `NRH04_SM_LinearResponse.wl` (the second-order expansion of the Gamma^2 density with
  four fluctuation channels and the exact variation identity on the BTZ background)

The `.nb` notebooks are generated cell-for-cell from the `.wl` sources by a kernel-side
converter; evaluating a notebook runs exactly the same checks.

This repository does not bundle a Wolfram runtime or license and does not contain an
automated hosted execution.  Independent reproduction therefore requires Mathematica
or Wolfram Engine and a fresh run of the command above.
