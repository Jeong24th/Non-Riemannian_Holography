# Non-Riemannian Holography — core checks

Five short Mathematica programs verify selected results of *Non-Riemannian Holography: Long Strings and Soft Hair*. Each follows the manuscript's notation and calculation order, with its inputs, calculation and expected zero residual together. There is no tensor toolbox to load.

Version **`submission-2026-10-14-core-v2`**. Manuscript SHA-256: `7BFD34C2D514F78EA169FD56C258A0052F022952AC39F36809F9D0E0634ACF0E`.

| File | What it checks | Manuscript |
|---|---|---|
| [01 Geometry](mathematica/01_Geometry.wl) | Generalized-metric constraints, rank and shared boundary | (1)–(2), (4a), (6) |
| [02 Hair and action](mathematica/02_HairAndAction.wl) | Reduced radial hair equation and finite saddle action | (8)–(9), SM2–SM5 |
| [03 Two-point kernels](mathematica/03_TwoPointKernels.wl) | SM27 Hessian normalization applied to the specified response kernels | (20), SM27, SM43, SM47–SM49, SM53–SM55 |
| [04 Strings and spinors](mathematica/04_StringsAndSpinors.wl) | Auxiliary-field elimination, clock forms, Clifford algebra and Hill bilinears | SM58, SM60–SM62, SM82, SM85, SM87 |
| [05 Abelian boundary](mathematica/05_AbelianBoundary.wl) | RR nilpotency, kinetic reduction and Cartan–Kosmann identity | SM94–SM108 |

Code uses `xPlus/xMinus` for the coordinates and `Oplus/Ominus` for flat indices; `calW`, `calB`, `calF`, `gamma3` and `barGamma3` denote the corresponding manuscript symbols. `calPi` is the product $L_+L_-$, while Mathematica's `Pi` is the numerical constant. File 02 retains $u_\pm=L_\pm^{-1/2}$ and the sum over $s=\pm1$. File 04 uses the SM4 $\lambda$ matrices; file 05 uses the separate SM5 $\gamma$ representation. No chirality or self-duality condition is imposed on the bispinor field.

This version corrects the first core release's SM58 multiplier orientation: $+\beta\bar\partial x^+ +\bar\beta\partial x^-$ gives $-2F\partial x^-\bar\partial x^+$. It also replaces the assumed reduced NR-spinor operator by the explicit SM85-to-SM87 Hill calculation. No independent verification of the full SM89 solutions is claimed here.

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

These are short checks of the displayed results under their stated inputs. File 02 takes the reduced radial EDFE and on-shell flux as inputs, then checks the integration and finite action. File 03 starts from the specified quadratic-action response operators and Lorentzian prescription. File 04 takes the Hill equations as inputs and derives their bilinear isometry relation. They do not rederive the full tensor EDFE, all counterterms or charges, a completed general-state Green function, or the full ten-dimensional Killing-spinor system. File 05 checks nilpotency and kinetic/Kosmann identities; it does not independently verify the Grassmann SUSY algebra, define correlators for null components, or establish a holographic duality.

The larger suite and its detailed equation ledger remain available in the [previous detailed release](https://github.com/Jeong24th/Non-Riemannian_Holography/tree/submission-2026-10-14). That release includes 46 fresh SM5 Python checks and clearly dated historical Mathematica results. Its long programs and intermediate calculations are intentionally absent from this core release. No manuscript or private reference files are included.
