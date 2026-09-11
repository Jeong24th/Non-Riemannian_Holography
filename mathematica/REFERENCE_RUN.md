# Verification execution record

## Recorded Mathematica run before the SM1 rewrite

Executed 2026-09-10 against manuscript SHA-256 `FCBE00570616A89820C3F969E8E97BA74454B8555DDA5EBC25ED56C268C5BE6F`.

- Environment: 13.2.1 for Microsoft Windows (64-bit) (January 27, 2023); `Windows-x86-64`.
- Command: `wolframscript -file mathematica/NRH00_RunAll.wl`.
- Command-line result: **270/270**, exit 0; 253.71 seconds.
- Regenerated notebook: `NotebookEvaluate[NRH00_RunAll.nb]` returned `True`,
  **270/270**; 252.97 seconds in the Mathematica front end.
- All nine notebooks have the same held input expressions as their `.wl` sources.
- Section counts: NRH02 30, NRH03 32, NRH04 72, NRH05 38, NRH06 32,
  NRH07 52, NRH08 14.
- The four existing manuscript-contract scripts also passed with that
  source snapshot supplied privately: NR linearization, R falloff, Gamma-squared action
  and doubled-yet-gauged reduction. No Python predicates were changed.

This is the recorded run on the named source snapshot, not a fresh execution on
the current manuscript. The Wolfram runtime was unavailable on the proofreading
host; all `.wl` and `.nb` files remain byte-identical.

The original 247 checks are retained. The 23 additions cover moving-frame momentum
variations, independent local fluctuation arguments, finite-state response through
exp(-2y/l), logarithmic finite coefficients and Ward-derived A/M kernels. Negative
controls distinguish frozen frames and the opposite source sign.

The exact symbolic data remain arbitrary functions in each check's stated domain.
The 5x5 matrix's undetermined remainder, local contacts and Green-function completion
are not evaluated. See EQUATION_LEDGER.md for per-equation restrictions, including
the scope of the ten-dimensional and Killing-spinor checks.

## Current manuscript proofreading and Python execution

Current source SHA-256: `C733463D16EA773C6DDA09143E44C1823F7C3CD46F7D771F14E53C8ED722C2DE` (2026-09-11).
Letter (1)-(22), SM (1)-(202); all 226 numbered displays retain their formulas
through the September 11 proofreading pass. All 62 bibliography entries follow
first citation order. The map was checked against the matching LaTeX build.

Algebra runs repeated on 2026-09-11 without the manuscript source:

| Check | Result |
|---|---|
| NR common-vacuum linearization | PASS |
| R fixed-frame exact projection | PASS |
| Gamma-squared action and normalization | PASS |
| Doubled-yet-gauged worldsheet reduction, strict dependency pin | PASS |
| One-sided reduced Killing system (`most-general`) | PASS |

The worldsheet LaTeX comparison passes against the current source. The other
three legacy comparisons fail: they look for the removed fixed-frame solution,
old response labels and earlier scope sentences. Their failures are not counted
as passes and do not validate the rewritten general-source derivation. No Python
verification code was changed or added in this release.

The older exact-projection scripts and the current saddle-frame tangent h use
different objects. Their coefficients must not be identified by name alone.
Likewise the public reduced spinor script has its own frame convention; its pass
does not certify the updated full ten-dimensional real basis.
