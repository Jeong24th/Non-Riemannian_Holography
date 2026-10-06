# Non-Riemannian Holography verification

Reproducibility calculations supporting *Non-Riemannian Holography: Long Strings and Soft Hair* by Shaun D. Hampton, Hyun-Cheol Kim, Jae-Hyuk Oh, and Jeong-Hyuck Park.

The October 6 update supplies **13 standalone Mathematica programs and matching standalone notebooks**, aligned with the revised five-section Supplemental Material. Every program embeds its own tensor tools and can run alone in an empty directory.

Source SHA-256: `F8AC68BB80B6C3D572160FD83705FE7506FCF277B00FB95021D39C2AF1E8E239` (2026-10-06). 139 numbered displays, including subequations.

The suite includes exact R/NR EDFE, asymptotic symmetry transformations, the Noether surface charge Q[epsilon], exact L_pm=0 linearized solutions, general-background asymptotic LEDFE checks, counterterms and finite renormalized action, and two-point functions. It preserves radial and charge derivations as supporting calculations after their sections were removed from the manuscript.

```sh
wolframscript -file mathematica/11_ZeroL_ExactLEDFE.wl
```

Use [the file guide](mathematica/README.md) to choose a calculation. Each matching `.nb` embeds the full source; it does not need the `.wl` or any other file. Mathematica 13.2.1 was used for the recorded execution. See [fresh execution results](mathematica/REFERENCE_RUN.md), [equation coverage](mathematica/EQUATION_LEDGER.md), and [current numbering](mathematica/MANUSCRIPT_MAP.md).

The checks retain their stated restrictions: general L_pm solutions are verified through exp(-2y/l), whereas the specified L_pm=0 solutions have exact full-radius residual checks. General homogeneous responses, contact terms, spatial quotient and NR interior completion are not uniquely determined. The regular R-vacuum K3 response and its separated-point hair correlator are explicitly calculated. A passing suite is not a claim that every equation or interpretation has been proved.

Historical Python checks under `checks/` and `evidence/` are retained for provenance; they are not the current release entry point and some manuscript-string guards target older sources. No new Python verification code is published in this update. The manuscript, protected coauthor notes and private working files are not included. `MANIFEST.sha256` hashes tracked payload bytes after Git normalization.

Fresh execution: **377/377 checks passed**, 13 programs, 13 standalone notebook payloads validated.
