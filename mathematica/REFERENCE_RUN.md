# Execution record — submission-2026-10-14

Manuscript SHA-256: `1DADD2A6F4A65E2D142A1086CE907176B8C4677951D877F3076A6077D03C9B55`.

| Evidence | Result | Scope |
|---|---|---|
| Current SM5 Python, Python 3.12.10 / SymPy 1.14.0 | **46/46 passed** | Fresh exact checks; [log](runs/2026-10-14/sm5-python.txt) |
| Mathematica programs 01–08 and 10–14 | **437 historical passes** | Their LF-normalized hashes match the October 8 executed files; not rerun October 14 |
| Rewritten Mathematica 09 | **Not executed** | No Mathematica installation on the preparation PC |
| 14 notebook/WL pairs | Payloads match | Source comparison only, not Mathematica or front-end execution |

The [current machine-readable receipt](RUN_RESULTS.json) records the exact code hashes and individual historical records. The [October 8 full receipt](runs/2026-10-08/RUN_RESULTS.json) has 451/451 passes, including 14 checks of the retired program 09. **451/451 is not a result for this submission version.** Its old source and notebooks are available at Git commit `ebf7f5b76046f379d84ecb10eed4d1d1a8e2e4ed`.

Reproduce the current Python checks and source-integrity checks with `python run_verification.py`. Add `--mathematica` to run every WL file in an isolated directory and fresh kernel. Current SM5 checks use arbitrary function jets, not numerical samples; the exterior algebra keeps field derivatives independent. The test includes both fermion kinetic reductions and wrong-sign controls. The coverage ledger distinguishes checked fixed-background components from cited general identities and uncomputed quantum/background-response questions.
