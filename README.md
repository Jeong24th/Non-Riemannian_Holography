# Non-Riemannian Holography — core checks

Five short Mathematica programs verify selected results of *Non-Riemannian Holography: Long Strings and Soft Hair*. Each is **31–39 lines**, with its assumptions, calculation and expected zero residual together. There is no tensor toolbox or repeated intermediate derivation to read.

Version **`submission-2026-10-14-core`**. Manuscript SHA-256: `31AE91B7E6FD951C532B84C77967A7F77F4E3D2BED463B3B414FD9C8AC928579`.

| File | What it checks | Manuscript |
|---|---|---|
| [01 Geometry](mathematica/01_Geometry.wl) | Generalized-metric constraints, rank and shared boundary | (2), (4a), (6) |
| [02 Hair and action](mathematica/02_HairAndAction.wl) | Reduced radial hair equation and finite saddle action | (8)–(9), SM2–SM5 |
| [03 Two-point kernels](mathematica/03_TwoPointKernels.wl) | Particular diagonal/mixed kernels and regular R-vacuum hair correlator | (20), SM43, SM47–SM49 |
| [04 Strings and spinors](mathematica/04_StringsAndSpinors.wl) | Auxiliary-field elimination, clock covectors and reduced chiral spinor channels | SM58, SM60–SM62, SM89 |
| [05 Abelian boundary](mathematica/05_AbelianBoundary.wl) | RR nilpotency, kinetic reduction and Cartan–Kosmann identity | SM94–SM108 |

## Mathematica

Run any file independently in a fresh kernel:

```sh
wolframscript -file mathematica/01_Geometry.wl
```

Or open its matching `.nb`: the short source is displayed as readable code, followed by one input cell that runs the adjacent `.wl`. Keep the pair in the same directory. The notebooks do not hide the calculation in a `ToExpression` string. A failed command-line check exits with status 1.

## Python counterpart and execution record

[verify_core.py](checks/verify_core.py) independently performs the same **27 algebraic checks** using SymPy. On Windows:

```powershell
python -m venv .venv-verify
.venv-verify\Scripts\python -m pip install -r requirements-verification.txt
.venv-verify\Scripts\python run_verification.py
```

On macOS/Linux substitute `.venv-verify/bin/python`. With an activated environment, use `python run_verification.py`; add `--mathematica` to run all five WL files in separate fresh kernels. Logs go to `verification-output/`.

**Recorded result: 27/27 Python checks passed** (Python 3.12.10, SymPy 1.14.0). The five rewritten Mathematica programs have **not been executed in Mathematica**, which is unavailable on the preparation PC. Notebook source matching is checked, but notebook front-end execution is not claimed. See [receipt](VERIFICATION.json) and [Python log](results/python.txt).

## Scope

These are short checks of the displayed results under their stated inputs. They do not rederive the full tensor EDFE, all counterterms or charges, a completed general-state Green function, or the full ten-dimensional Killing-spinor system. File 05 does not define correlators for null components or establish a holographic duality.

The larger suite and its detailed equation ledger remain available in the [previous detailed release](https://github.com/Jeong24th/Non-Riemannian_Holography/tree/submission-2026-10-14). That release includes 46 fresh SM5 Python checks and clearly dated historical Mathematica results. Its long programs and intermediate calculations are intentionally absent from this core release. No manuscript or private reference files are included.
