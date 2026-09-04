# Mathematica Reference Execution

This record describes the latest complete reference execution applicable to the
symbolic predicates in this directory.

**Verified 2026-09-04.** Both the canonical `.wl` runner and the regenerated
`NRH00_RunAll.nb`, evaluated through the Mathematica front end, passed 247/247.
The ten-dimensional identities retain arbitrary chiral `L_±(x^±)` on the
Riemannian branch and arbitrary `L_+(x^+)`, `W_0`, `W_1` on the one-sided
non-Riemannian branch. No constant-profile fallback was needed.

- Environment: Mathematica 13.2.1 for Microsoft Windows (64-bit),
  January 27, 2023; `$SystemID = Windows-x86-64`; `wolframscript`
- Entry point: `NRH00_RunAll.wl`
- Command: `wolframscript -file NRH00_RunAll.wl`
- Result: 247/247 exact checks passed
  (NRH02 30, NRH03 32, NRH04 49, NRH05 38, NRH06 32, NRH07 52, NRH08 14)
- Command-line result: exit status 0; 191.30 s wall time
- Notebook result: `NotebookEvaluate` returned `True`; 247/247; 190.97 s
  (including front-end evaluation and evaluated-notebook export)
- Exit behavior: a failed check produces a nonzero command-line exit status

The `.nb` notebooks are generated cell-for-cell from the `.wl` sources by a kernel-side
converter. All nine were regenerated for this execution. Their concatenated Input
cells match the held `.wl` expressions exactly after removing only top-level Nulls
from comment-only/blank lines. The runner notebook loads the canonical section
files, just as the command-line runner does.

The first execution of the symbolic-only revision found 246/247. Its sole failure
exposed an assembly error in `NRH07`: the two flat `R^4` identity blocks were placed
off-diagonally, giving `H_R4 = J_R4` instead of the intended `H_R4 = I_8`.
Moving those two blocks onto the diagonal restores the intended geometry and
makes the arbitrary-chiral projected Ricci identity vanish exactly. No curvature
routine, zero predicate, profile, or check count was weakened. The results above
apply to the corrected assembly; they do not validate the erroneous earlier one.

This repository does not bundle a Wolfram runtime or license and does not contain an
automated hosted execution.  Independent reproduction therefore requires Mathematica
or Wolfram Engine and a fresh run of the command above.
