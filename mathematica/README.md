# Mathematica verification suite

Symbolic Mathematica verification of the displayed equations of

> S. D. Hampton, H.-C. Kim, J.-H. Oh and J.-H. Park,
> *Non-Riemannian Hair in Long-String Holography* (Letter + Supplemental Material),

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

If a term, an abbreviation, or a code symbol is unfamiliar, see the sections
**How to read the files**, **Notation and conventions**, **Glossary of abbreviations**,
and **Code symbol dictionary** below the file table.

## Files

| File | Verifies | Contents |
|---|---|---|
| `NRH00_RunAll.wl` | — | runs the whole suite, grand PASS/FAIL summary |
| `NRH01_DFT_Tools.wl` | shared toolbox | O(D,D) metric, projectors, torsionless semi-covariant connection, curvatures S_ABCD / Ricci / S_(0) (closed form **and** contraction form), Einstein curvature G_AB, generalized Lie derivative, Gamma^2 density, boundary vector B^A, the generalized-metric momentum A^K_MN of the Gamma^2 boundary term, double-vielbein spin connection, PASS/FAIL framework |
| `NRH02_Letter_Riemannian.wl` | Letter labels `Rfields`–`Rcorrelators` (Eqs. (1)–(11)); SM `SMinvariantstress` | Banados family in DFT variables; **G_MN = 2 l^-2 J_MN exactly for arbitrary chiral L±(x±)**; `boundaryH` and `Rboundaryframe`; `RDFTconservation` and `Rcontinuity`; `Rkilling`–`RVirasoro`; the coefficient matching behind `RKdef`; c = 3l/2G; **Eq. (9a) applied nonlinearly to the exact family reproduces the full one-point matrix of Eq. (11)/SM (36) and the vanishing scalar response**; horizon and Gomis–Ooguri channel structure |
| `NRH03_Letter_NonRiemannian.wl` | Letter labels `NRvariables`–`NRasympt` (Eqs. (12)–(16)); SM labels `NRhill`–`NRGprofile`, `SMBtransform`, `SMWshift`, `SNCtau`, `SMWisB`, `SMNRlocalstabilizer` | exact non-Riemannian matrix and dilaton; **the full EDFE with generic W(x+,x-,chi) collapses to `NRchiODE`, and `NRWgeneral` solves it** for arbitrary chiral L± and arbitrary W₀, W₁; SNC clocks; B-shift removal of W₀; anomaly-free δL± and the third derivative in δW₁; the weight-one law of W₀; the sinh(χ/2) expansion and the Hill identity |
| `NRH04_SM_LinearResponse.wl` | SM 1 labels `SMcosetreconstruction`–`SMNRhairhessian` (SM (4)–(45)) | vielbein variations; **the exact Γ²-variation identity SM (11)–(12) verified on BTZ with generic fluctuations**, and its on-shell reduction to SM (13); limiting frame and flat metrics; fixed-frame projections `SMexactprojectionR`/`SMexactprojectionNR`; **`SMNRlin` re-derived from the nonlinear G_MN**; `SMNRsol`, `SMlogfreeconditions`, `SMrwquadratic`/`SMrwvariation`; **the complete four-channel boundary variation `SMSren2`** and `SMscalarcutoffresponse`; `SMresponsematrix`, canonical-source factor and Hessians, one-point functions (also from Eq. (9a) nonlinearly on the NR family), PBH kernel and two-point normalization |
| `NRH05_SM_Charges_Action.wl` | SM 2–3 labels `SMCPSform`–`SMgamma2value` (SM (47)–(59)) | Noether surface potential; Brown–Henneaux cocycle; **`SMCPSresult` with every Θ̂-component vanishing and W₁ dropping out through `SMNRchargecancellation`**; `SMCPSalgebra` and its footnote; the charge Poisson algebra; `SMgamma2`, `SMgamma2flux`, cutoff and renormalized value |
| `NRH06_SM_Worldsheet.wl` | SM 5 labels `SMRfirstorder`–`SMBRSTfusion` (SM (74)–(92)) | Riemannian first-order form and `SMceff`; SNC kernel and `SMdygconstraints`/`SMdygGO`; `SMGO`, `SMvertex`, `SMlongstringE`, `SMdeltaL`; `SMWZWlevel`; Fradkin–Tseytlin radial weight, marginal roots and the central charge `SMBRSTcentral`; the gauge-obstruction derivation `SMWgaugeobstruction`; the fusion exponent and resonances `SMBRSTfusion` |
| `NRH07_SM_KillingSpinors.wl` | SM 6 labels `SMupliftblocks`–`SMhairyKS` (SM (93)–(119)) | uplift checks for both factors and exact ten-dimensional probes; `SMRlocaliso`; `SMexactiso` and `SMweighteddilaton`; `SMNRlocalstabilizer`; **the exact non-Riemannian double vielbein `SMvielbein`/`SMvielbeincheck`**; Clifford/Majorana/barred-algebra package; `SMkillingspinor`, `SMreducedDirac`, `SMinternalprojectors` (including the S³ spinor integrability), `SMspinorcountchain`, and `SMhairyKS`; **the Killing-spinor bilinear of `SMsusyclosure` is an exact vacuum isometry `SMexactiso`**; Hill-system and monodromy checks |
| `NRH08_SM_BoundaryCandidate.wl` | SM 7 labels `SMcandidateD`–`SMcandidatetrivialL` (SM (120)–(130)) | projector selection and candidate action; non-abelian gauge invariance; commuting Witt symmetries; boundary isometry ⇔ chirality; chiral fermionic transformations with exact Grassmann arithmetic; the extra chiral transformation; the equation-of-motion redundancies |

Coverage is keyed to LaTeX labels rather than snapshot-dependent equation numbers.
`MANUSCRIPT_MAP.md` gives the current number of every label used by the suite and pins
the exact manuscript source by SHA-256.  `EQUATION_LEDGER.md` walks through **every
numbered display of the paper in order** — Letter (1)–(16) and SM (1)–(130) — and names,
for each one, the check(s) that verify it, or states that it is a definition or a cited
statement with nothing to compute.

The suite currently comprises **247 exact checks**; the recorded reference execution
passes 247/247.  See `REFERENCE_RUN.md` for the neutral environment, command,
result, and reproducibility limitations.

## Running

Requires Mathematica (tested with 13.2) or the free Wolfram Engine.

```
wolframscript -file NRH00_RunAll.wl
```

runs everything (a few minutes) and exits nonzero if any check fails.  Each
`NRH0k_*.wl` file is self-contained (it loads `NRH01_DFT_Tools.wl` from its own
directory) and can be run alone the same way, or opened and evaluated in the Mathematica
front end, where the cell markers render titles and explanations as a notebook.

**Notebooks.**  Every file is also provided as a double-clickable notebook with the
identical content: download the whole `mathematica/` folder, open `NRH00_RunAll.nb`
(or any section notebook, e.g. `NRH03_Letter_NonRiemannian.nb`) in Mathematica, and use
*Evaluation → Evaluate Notebook*.  The section notebooks locate `NRH01_DFT_Tools.wl` in
the same directory, so keep the folder together.  The `.nb` files are generated from the
`.wl` sources, which remain the canonical, diff-able versions; the two always carry the
same checks.

## How to read the files

* Each `.wl` file is organized like a short paper: markers such as `(* ::Title:: *)`,
  `(* ::Section:: *)`, `(* ::Text:: *)` render as headings and explanatory prose when
  the file is opened in Mathematica.  The `.nb` notebooks show the same content with
  the formatting already applied.
* Every verification is a single call of the small framework loaded from
  `NRH01_DFT_Tools.wl`:
  - ``NRH`CheckZero["label", expr]`` proves that `expr` (a scalar, a list, or a whole
    matrix) is **identically zero** as a symbolic expression;
  - ``NRH`Check["label", statement]`` records a structural true/false claim
    (a rank, a set equality, a "contains no ..." statement).
  Each check prints `[PASS]` or `[FAIL]`; each file ends with an `n/n` summary, and
  `NRH00_RunAll.wl` prints the grand total.  Any failure makes a command-line run exit
  with a nonzero code.
* "Exact" means: arbitrary functions stay arbitrary.  When a check says it holds "for
  arbitrary chiral L±(x±)", the functions `Lp[xp]`, `Lm[xm]` are never specialized —
  the zero is an identity of symbolic algebra, not a numerical coincidence.
* "Modulo total derivatives": densities that sit under a boundary integral are compared
  with an **Euler–Lagrange test** — a density is a total derivative if and only if all
  of its variational (Euler–Lagrange) derivatives vanish.  This is how statements like
  "the cocycle is a total derivative" are decided exactly; in two boundary variables the
  built-in `VariationalD` of the VariationalMethods package is used.
* "Nonlinear" versus "linearized": the response dictionary Eq. (9a) is evaluated in two
  independent ways — through the quadratic action and the linearized solution (SM 1),
  and directly from the full DFT connection of each exact saddle (the momentum A^y_MN
  projected on the fixed boundary frame).  Both give the same one-point functions.

## Notation and conventions (following the paper)

* **Coordinates.**  Boundary lightcone x^± = (t ± lφ)/√2 with φ ~ φ + 2π; y is the
  holographic radial coordinate and the boundary sits at y → ∞; l is the AdS₃ radius.
  Three interchangeable radial variables appear in the code, chosen to keep every
  computation rational: `u` = e^{2y/l} (Riemannian side), `z` = e^{−2y/l} = 1/u
  (charge falloffs), and the non-Riemannian variable `ch` = χ with
  `T` = e^{χ/(2√2)} rationalizing all hyperbolic functions.
* **Doubled indices.**  The fixed coordinate order is
  x^M = (x̃₊, x̃₋, ỹ; x⁺, x⁻, y): the three dual ("winding") coordinates first, then
  the three physical ones.  All 6×6 matrices use this order (4×4 on the boundary,
  20×20 in the ten-dimensional probes).  The section condition sets all dual
  derivatives to zero, which is what the derivative operators in the code implement.
* **Frame (flat) indices.**  The three-dimensional lightcone frame is (⊕, ⊖, y) with
  flat metric η = ((0,−1,0),(−1,0,0),(0,0,1)); the barred frame has η̄ = −η.  In code
  the symbols ⊕ and ⊖ are written `p`(lus)/`op` and `m`(inus)/`om`: for instance
  `hpm` is the fluctuation component h_{⊕⊖̄} and `Kmp` is the response K_{⊖⊕̄}.
* **State data.**  L₊(x⁺) and L₋(x⁻) are the two chiral Banados (stress-tensor)
  functions.  W₀(x⁺,x⁻) and W₁(x⁺,x⁻) are the two radial modes of the hair
  W = W₀ + e^{−2y/l}W₁ + (inhomogeneous terms): W₀ is the non-normalizable marginal
  source (locally pure gauge), W₁ the normalizable "soft hair".
* **Derived variables.**  ψ± := L±^{−1/2} (Hill variables), Π = L₊L₋,
  q = e^{−2y/l}√(Π/2), χ = 2√2 arctanh q.  The Hill equation (l²/2)s″ = L s governs
  the Riemannian Killing spinors; A± = ψ±″/ψ± is its potential data.

## Glossary of abbreviations

| Abbreviation | Meaning |
|---|---|
| DFT | Double Field Theory: the metric, Kalb–Ramond field B, and dilaton unified as one object on doubled coordinates, with T-duality manifest |
| O(D,D), J | the T-duality group and its invariant metric J_MN (off-diagonal unit blocks in our basis) |
| H_MN | the generalized metric — the O(D,D) tensor packaging (g, B); it obeys the constraint H J H = J |
| d | the DFT dilaton: e^{−2d} = √(−g) e^{−2φ} on Riemannian backgrounds |
| P, P̄ | the projectors (J ± H)/2; their mixed components carry the physical fluctuations |
| EDFE | Einstein Double Field Equations G_MN = T_MN; on these saddles G_MN = 2 l^{−2} J_MN |
| Γ (semi-covariant connection) | the torsionless DFT analogue of the Christoffel symbols, determined by H and d |
| S_ABCD, S_AB, S₍₀₎ | the semi-covariant Riemann, Ricci, and scalar curvatures of DFT |
| ĥL (generalized Lie derivative) | the DFT gauge transformation combining diffeomorphisms and B-field gauge shifts |
| R / NR | the two exact saddles: the Riemannian (Banados/BTZ) branch and the everywhere non-Riemannian branch |
| type (n, n̄) | the Morand–Park classification of non-Riemannian generalized metrics; type (1,1) is the Gomis–Ooguri geometry |
| GO | Gomis–Ooguri: the non-relativistic string theory describing the long-string sector |
| SNC | string Newton–Cartan geometry; τ± are its clock one-forms, m_μ^± its mass gauge fields |
| BTZ | the Banados–Teitelboim–Zanelli black hole (constant L± > 0) |
| BH | Brown–Henneaux: the asymptotic Virasoro symmetry of AdS₃ with central charge c = 3l/2G |
| PBH | Penrose–Brown–Henneaux: the radial completion of a boundary diffeomorphism into the bulk |
| FG | Fefferman–Graham gauge for the radial expansion |
| GKPW | the Gubser–Klebanov–Polyakov–Witten source/partition-function dictionary |
| Ward identities | the boundary conservation laws of the response tensor, Eq. (10) of the Letter |
| KS | Killing spinor; SDFT = supersymmetric DFT; MW = Majorana–Weyl spinor conditions |
| Hill equation | (l²/2)s″ = L(x)s: the periodic-coefficient ODE controlling the Riemannian Killing spinors and their monodromy |
| DYG | the doubled-yet-gauged worldsheet action, which couples the string directly to H_MN |
| FT | the Fradkin–Tseytlin worldsheet dilaton coupling (it improves the stress tensor by −l^{−1}∂²y) |
| OPE | operator product expansion (worldsheet CFT) |
| c_eff | the effective Gomis–Ooguri speed of light, c_eff² = 2F = 2(e^{2y/l} + L₊L₋e^{−2y/l}) |
| C-bracket | the DFT bracket of doubled gauge parameters (the analogue of the Lie bracket) |
| PRRS | Park–Rey–Rim–Sakatani: the DFT covariant-phase-space construction of surface charges (arXiv:1507.07545) |
| K^{AB}, Θ, k_ξ | Noether surface potential, symplectic potential, and the surface-charge one-form built from them |
| BT | Barnich–Troessaert: the adjustment of brackets for field-dependent asymptotic parameters |
| Witt algebra | the centerless Virasoro algebra of chiral reparametrizations |
| GF cocycle | Gelfand–Fuchs: the Virasoro two-cocycle ∮ε₁∂³ε₂, whose cocycle condition is the Jacobi identity of the centrally extended algebra |
| JW | Jordan–Wigner: the matrix realization of anticommuting (Grassmann) generators used for the fermionic checks |
| μ (mu-dichotomy) | the endpoint parameter of the on-shell action: μ = L₊L₋ on the R branch, μ = L₊L₋/2 on the NR branch |
| k = l²/α′ | the flux (WZW) level fixed by the S³ NS–NS flux; α′ is the string tension parameter, w the winding number |
| SM (n) / Eq. (n) | equation numbers in the source build pinned by `MANUSCRIPT_MAP.md`; use the quoted LaTeX labels as stable identifiers |

## Code symbol dictionary

| In the code | In the paper |
|---|---|
| `Lp[xp]`, `Lm[xm]` | the chiral Banados data L₊(x⁺), L₋(x⁻) |
| `ep[xp]`, `em[xm]` (also `e1f`, `e2f`, `e1g`, `e2g`, `al`) | asymptotic-symmetry parameters ε⁺(x⁺), ε⁻(x⁻) (and the off-shell PBH test function α⁺) |
| `W0[xp,xm]`, `W1[xp,xm]` | the hair modes W₀, W₁ |
| `psip[xp]`, `psim[xm]` | ψ± = L±^{−1/2} |
| `u`, `z`, `ch`, `T`, `LT` | e^{2y/l}, e^{−2y/l}, χ, e^{χ/(2√2)}, log T (kept as an independent symbol) |
| `hpp`, `hpm`, `hmp`, `hmm`, `dd` | the fixed-frame fluctuations h_{⊕⊕̄}, h_{⊕⊖̄}, h_{⊖⊕̄}, h_{⊖⊖̄} and δd |
| `Kpp`, `Kpm`, `Kmp`, `Kmm`, `T0` | the boundary responses K_{⊕⊕̄}, K_{⊕⊖̄}, K_{⊖⊕̄}, K_{⊖⊖̄} and T₍₀₎ |
| `Hinf`, `dinf`, `H0` | the boundary data H^∞, d^∞ = −y/l and the induced 4×4 boundary generalized metric H⁽⁰⁾ |
| `Vinf`, `Vbinf`, `V0`, `Vb0`, `eta3`/`etab3` | the aligned double-vielbeins V^∞, V̄^∞, their boundary restrictions, and the flat metrics η, η̄ |
| `GammaDFT`, `RiemannR4`, `RicciS`, `ScalarS0`, `EinsteinG` | Γ_CAB, R_CDAB, S_AB, S₍₀₎, G_AB |
| `GenLieH`, `GenLieD` | ĥL_ξ H_MN and ĥL_ξ d |
| `Gamma2Density`, `GammaBVector`, `MomentumAK` | the Γ²-Lagrangian density, its boundary vector B^M, and the generalized-metric momentum A^K_MN (frame components A^K_{pq̄} of SM (12)) |
| `VexU`, `VbexU` (NRH07) | the exact non-Riemannian double vielbein of SM (106) with raised local indices |
| `E1`, `E2`, `Cmat`, `XM` (NRH07) | two vacuum Killing spinors, the Majorana conjugation C = iσ₂, and their bilinear X^M |
| `H4`, `Lq4`, `ThetaY` (NRH04) | the four-channel coset family, its quadratic Γ² density, and the radial boundary momentum at the cutoff |
| `NoetherK`, `KhatComp`, `ThetaHat`, `kPlus`/`kMinus` | K^{AB}, K̂^{AB}, e^{−2d}Θ̂^A, and the finite charge one-forms k^{∓y}[ε^±] |
| `CBracket` | the DFT C-bracket of doubled vectors |
| `xiPlus`, `xiMinus`, `xiIso`, `xiPBH` | the asymptotic generators ξ[ε±], the exact vacuum isometries ξ[v±, ω±], and the PBH vector ξ[α⁺] |
| `tauP`, `tauM`, `esig` | the SNC clock forms τ± and e^σ = √(L₊/L₋) |
| `FF`, `rho`, `GG`/`Gp` | the radial source F, the kernel ρ(χ), and the profile G(χ) with its derivative |
| `Gam[a]`, `GamBar[a]`, `Gam11`, `BB10` | the ten-dimensional Γ^{p̂}, the barred Γ̄^{p̂̄} = Γ^{p̂}Γ¹¹, the chirality Γ¹¹, and the Majorana intertwiner 𝖡₁₀ |
| `th[i]` | the Jordan–Wigner Grassmann generators θᵢ |
| `phiP`/`phiM`, `psiP`/`psiM`, `Ap`/`Am`, `Lam` | the boundary-candidate fields φ^±, ψ^±, the internal gauge field 𝔸_±, and the gauge parameter Λ |
| ``NRH`CheckZero``, ``NRH`Check`` | "this expression is identically zero" / "this structural statement is true" |

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
* **Radial frame convention.**  The suite uses the radial-frame sign convention of the
  manuscript's SM (14) (radial columns (+,+) for V and (+,−) for V̄); with it the toolbox
  reproduces the displayed vacuum spin-connection components, the Killing spinor
  E = (√2 f, l∂₊f)ᵀ with eigenvalue +1/(√2 l), and the reduced matrix diag(1,−1) of SM (116).

## Reference execution

The executable identities are accompanied by the environment and result record in
`REFERENCE_RUN.md`.  The Python scripts in `checks/` provide a separate implementation
of many of the same formulas.
