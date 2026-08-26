# Mathematica verification suite

Symbolic Mathematica verification of the displayed equations of

> S. D. Hampton, H.-C. Kim, J.-H. Oh and J.-H. Park,
> *Long Strings and Non-Riemannian Hair* (Letter + Supplemental Material),

organized **in the order of the paper**: the two Letter sections first, then the seven
Supplemental Material sections.  It complements the Python/SymPy scripts in `checks/`
with an independent implementation in a second computer-algebra system, in the spirit
of the co-author Mathematica notebooks used within the collaboration (not part of this
archive).

Everything here is **exact symbolic computation** — no floating point.  Where a fully
general symbolic run would be prohibitively slow (only the ten-dimensional curvature
probes), the same exact pipeline is evaluated in exact **rational arithmetic** at
explicitly quoted rational parameter values, matching the probe style quoted in the SM
text.

## Files

| File | Verifies | Contents |
|---|---|---|
| `NRH00_RunAll.wl` | — | runs the whole suite, grand PASS/FAIL summary |
| `NRH01_DFT_Tools.wl` | shared toolbox | O(D,D) metric, projectors, torsionless semi-covariant connection, curvatures S_ABCD / Ricci / S_(0) (closed form **and** contraction form), Einstein curvature G_AB, generalized Lie derivative, Gamma^2 density, boundary vector B^A, double-vielbein spin connection, PASS/FAIL framework |
| `NRH02_Letter_Riemannian.wl` | Letter Eqs. (1)–(10) | Banados family in DFT variables; **G_MN = 2 l^-2 J_MN exactly for arbitrary chiral L±(x±)**; boundary data (4); Ward identities (7); asymptotic symmetry (8)–(9) with the anomalous term; c = 3l/2G; two-point normalization; horizon; Gomis–Ooguri channel structure |
| `NRH03_Letter_NonRiemannian.wl` | Letter Eqs. (11)–(16), SM 4 | the exact non-Riemannian matrix (12)–(13); **the full EDFE with a generic hair profile W(x+,x-,chi) collapses to the single radial ODE d²W/dχ² = F of SM (66)–(67), and Eq. (15) solves it** — for arbitrary chiral L±(x±) and arbitrary W₀, W₁(x⁺,x⁻); the SNC clock forms (14); the B-shift removal of W₀; the asymptotic laws (16): δL± anomaly-free, the third derivative in δW₁ |
| `NRH04_SM_LinearResponse.wl` | SM 1 | aligned frame SM (10)–(13); exact fixed-frame projections SM (29)–(30), (44)–(45) incl. h⁽²⁾_{⊕⊖̄} = W₁/2; **the linearized system SM (17) re-derived from the nonlinear G_MN**; its general solution SM (19); log-free conditions SM (20); the crossed Γ²-quadratic action SM (22)–(23) with the coset completion; response normalization SM (31) and the one-point functions; PBH data SM (34) and the two-point normalization (8π)⁻²(c/2) |
| `NRH05_SM_Charges_Action.wl` | SM 2, SM 3 | Noether surface potential K^{AB} (PRRS Eq. (A.4)); Riemannian potential components SM (50); the Brown–Henneaux cocycle; **the NR charge one-form SM (52): k = (4/l) ε δL± with every Θ̂-component vanishing and W₁ dropping out componentwise (SM (53)–(54))**; centerless bracket closure SM (56); the Γ² identity SM (57); the flux SM (58) with the μ-dichotomy; the renormalized action SM (60)–(61); the two interior endpoints |
| `NRH06_SM_Worldsheet.wl` | SM 5 | first-order form SM (75) ↔ E = g − B; c_eff SM (76); the clock kernel and the exact reduced NR Lagrangian routes SM (78)–(79), (94); the GO vertex coefficient SM (81); the winding-probe energy SM (82) with √(−det g₂) = l e^{−2d}; the flux level SM (95); the Fradkin–Tseytlin radial weight SM (96)–(98) with marginal roots {1, e^{−2y/l}}; the contraction bookkeeping behind exact marginality of V_W[W₀] |
| `NRH07_SM_KillingSpinors.wl` | SM 6 | uplift checks SM (101)–(102) factor by factor; **exact rational ten-dimensional probes: S₀⁽¹⁰⁾ = 0 and (P S P̄)⁽¹⁰⁾ = 0 for both uplifts, including one-sided hair**; exact vacuum isometries SM (103)–(104); the stabilizer system SM (105); the complex 3+3+4 Clifford representation, Majorana intertwiner and barred algebra SM (106)–(110) with the count chain; **the vacuum Killing spinor with arbitrary chiral profile SM (111)/(116), both internal channels**; the rank-six hairy jet system SM (117); the Hill system SM (112)–(114) and the global counts by monodromy; the complementary-halves bookkeeping |
| `NRH08_SM_BoundaryCandidate.wl` | SM 7 | projector selection in the candidate action SM (118)–(119); non-abelian gauge invariance SM (120); the Witt symmetries SM (121)–(122); boundary isometry ⇔ chirality; **the chiral fermionic transformations SM (123)–(124) with exact Grassmann arithmetic (Jordan–Wigner realization)**; the extra chiral transformation SM (125)–(126) |

Equation numbers refer to the 2026-08-26 build of the manuscript (Letter Eqs. (1)–(16),
SM (1)–(125)); the LaTeX labels quoted inside the files are stable across rebuilds.

The suite currently comprises **169 exact checks**; the reference run passes 169/169.

## Running

Requires Mathematica (tested with 13.2) or the free Wolfram Engine.

```
wolframscript -file NRH00_RunAll.wl
```

runs everything (about ten minutes) and exits nonzero if any check fails.  Each
`NRH0k_*.wl` file is self-contained (it loads `NRH01_DFT_Tools.wl` from its own
directory) and can be run alone the same way, or opened and evaluated in the Mathematica
front end, where the cell markers render titles and explanations as a notebook.

## Method notes

* **Rational radial variables.**  The Riemannian saddle is handled with u = e^{2y/l}, so
  every identity is a rational-function statement that Mathematica decides exactly.  The
  non-Riemannian saddle is handled in the radial variable chi itself, with the exact
  chain rules dchi/dy = −(2√2/l) sinh(chi/√2) and dchi/dx^± = −√2 sinh(chi/√2) ψ±'/ψ±
  (ψ± := L±^{−1/2}); hyperbolic functions are rationalized by chi → 2√2 log T, under
  which all identities become Laurent-polynomial statements in T and T^{2√2}.
* **The connection.**  The torsionless semi-covariant connection is implemented with its
  trace vector fixed by the defining dilaton compatibility Γ^B_{BA} = −2∂_A d (see the
  implementation note in `NRH01_DFT_Tools.wl`).  All defining properties — ∇P = 0, the
  dilaton trace, the vanishing totally antisymmetric part — and the curvature identities
  (pair symmetry, algebraic Bianchi, projective property, the commutator identity) were
  used to validate the implementation; the EDFE checks then close exactly, including the
  scalar curvature computed two independent ways (closed form vs. S_ABCD contraction).
* **Grassmann arithmetic** in `NRH08` is exact: the four odd parameters/fields are
  realized as Jordan–Wigner matrices, so anticommutation is ordinary matrix algebra and
  no sign bookkeeping is done by hand.

## Provenance

These files were drafted with AI assistance (Anthropic Claude, Fable 5) in the
AI-assisted cross-verification workflow described in the paper's Acknowledgments, and
machine-executed end to end (the 169/169 reference run above).  The physics content
they check is that of the manuscript; the Python guard scripts in `checks/` remain the
primary reproducibility archive, and the two implementations are independent of each
other.
