# Non-Riemannian Holography: verification

Executable calculations for *Non-Riemannian Holography: Long Strings and Soft Hair* by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

Submission version: **`submission-2026-10-14`**. Manuscript SHA-256: `1DADD2A6F4A65E2D142A1086CE907176B8C4677951D877F3076A6077D03C9B55`. This repository contains verification code; the manuscript and private reference material are not included.

## Run

Python 3.12 and SymPy 1.14:

```powershell
python -m venv .venv-verify
.venv-verify\Scripts\python -m pip install -r requirements-verification.txt
.venv-verify\Scripts\python run_verification.py
```

On macOS/Linux use `.venv-verify/bin/python` in place of `.venv-verify\Scripts\python`. With an already activated environment, the command is simply `python run_verification.py`.

The default run validates file hashes and all 14 embedded notebook payloads, then executes the current SM5 checks in an empty temporary directory. It writes logs and `results.json` to `verification-output/` and exits nonzero on failure.

With Mathematica/WolframScript installed, run the full suite in separate fresh kernels:

```sh
python run_verification.py --mathematica
# Or run any standalone calculation directly:
wolframscript -file mathematica/14_SM29_IntrinsicCounterterm.wl
```

Use `--wolframscript PATH` if needed. Each `.wl` contains all its definitions; each `.nb` embeds the matching program. No manuscript input, private package or download is required by a calculation.

## What to run for each section

| Manuscript | Programs | Scope |
|---|---|---|
| Letter 2–3; SM1 | 01, 02, 03 | Exact R/NR saddles, EDFE, radial hair solution and renormalized saddle action |
| Letter 4; SM2 | 04, 05, 10–14 | Linearized solutions, counterterms, asymptotic transformations and specified two-point kernels |
| Letter 4 charges | 06 | Covariant surface charges and central terms |
| SM3 | 07 | Worldsheet reductions, clock covectors, boundary and radial vertex checks |
| SM4 | 08 | Uplift equations, Clifford algebra and specified reduced Killing-spinor sectors |
| SM5 | [Python](checks/verify_sm5_abelian_rr.py), 09 | Abelian RR complex, dilaton cancellation, Kosmann/Witt and physical supersymmetry |

The [file guide](mathematica/README.md), [equation map](mathematica/MANUSCRIPT_MAP.md) and [coverage ledger](mathematica/EQUATION_LEDGER.md) give precise locations and limitations.

## Execution status

- **Current SM5 Python: 46/46 exact checks passed**, Python 3.12.10 / SymPy 1.14.0. Odd fields and their derivatives are independent exterior-algebra generators. Wrong-sign negative controls are included.
- **13 unchanged Mathematica programs: 437 historical passing checks** from October 8, with matching LF-normalized hashes. They were not rerun on October 14.
- **Rewritten Mathematica 09: not executed.** Mathematica is unavailable on the preparation PC. Notebook payload equality is a source-integrity check, not kernel execution.

See the [receipt](mathematica/RUN_RESULTS.json) and [logs/status](mathematica/REFERENCE_RUN.md). The older 451/451 total includes the retired non-Abelian model and is not the current suite's result.

General-state interior conditions, homogeneous completion, contacts and spatial quotient data are not fixed by the particular two-point checks. SM30 is checked at quadratic order on specified vacua and fails as a universal general-state counterterm. SM5 verifies a classical propagating CFT; it does not define correlators of fields absent from the fixed-background action or prove a holographic duality. Historical files under `checks/` and `evidence/`, other than the current SM5 entry point, remain secondary provenance and can target earlier manuscript versions.
