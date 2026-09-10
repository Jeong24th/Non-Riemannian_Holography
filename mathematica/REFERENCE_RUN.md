# Mathematica reference execution

Verified 2026-09-10 against manuscript SHA-256 `FCBE00570616A89820C3F969E8E97BA74454B8555DDA5EBC25ED56C268C5BE6F`.

- Environment: 13.2.1 for Microsoft Windows (64-bit) (January 27, 2023); `Windows-x86-64`.
- Command: `wolframscript -file mathematica/NRH00_RunAll.wl`.
- Command-line result: **270/270**, exit 0; 253.71 seconds.
- Regenerated notebook: `NotebookEvaluate[NRH00_RunAll.nb]` returned `True`,
  **270/270**; 252.97 seconds in the Mathematica front end.
- All nine notebooks have the same held input expressions as their `.wl` sources.
- Section counts: NRH02 30, NRH03 32, NRH04 72, NRH05 38, NRH06 32,
  NRH07 52, NRH08 14.
- The four existing manuscript-contract scripts also passed with the current
  source supplied privately: NR linearization, R falloff, Gamma-squared action
  and doubled-yet-gauged reduction. No Python predicates were changed.

The original 247 checks are retained. The 23 additions cover moving-frame momentum
variations, independent local fluctuation arguments, finite-state response through
exp(-2y/l), logarithmic finite coefficients and Ward-derived A/M kernels. Negative
controls distinguish frozen frames and the opposite source sign.

The exact symbolic data remain arbitrary functions in each check's stated domain.
The 5x5 matrix's undetermined remainder, local contacts and Green-function completion
are not evaluated. See EQUATION_LEDGER.md for per-equation restrictions, including
the scope of the ten-dimensional and Killing-spinor checks.
