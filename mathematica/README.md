# Mathematica verification

Exact symbolic checks for *Non-Riemannian Holography: Long Strings and Soft Hair*.
The `.wl` files are canonical; each `.nb` notebook contains the same input expressions.
The files follow the order of the manuscript: one file for the Letter, then one file per
Supplemental Material section (SM 3 is split into its linearized-dynamics, renormalization
and charge parts). Every check prints the manuscript label(s) it verifies, so a label from
[EQUATION_LEDGER.md](EQUATION_LEDGER.md) can be located in the sources and in the run log by
text search.

```sh
wolframscript -file NRH00_RunAll.wl
```

Or keep this directory together, open `NRH00_RunAll.nb`, and choose
**Evaluation > Evaluate Notebook**. Section files also run independently.
A failed check returns a nonzero command-line exit status. The full run takes about five
minutes with Mathematica 13.2; see [REFERENCE_RUN.md](REFERENCE_RUN.md).

| File | Manuscript part | Checks |
|---|---|---:|
| NRH00 | runner | — |
| NRH01 | shared tools: O(3,3) metric, torsionless connection, curvatures, Γ² density and flux, radial momenta, generalized Lie derivative, both exact saddles as series in u = e^{−2y/l}, frame variation, projected linearized EDFE | — |
| NRH02 | Letter (1)–(22): both saddles, boundary data and frames, EDFE, Ward identities, asymptotic symmetries and the Virasoro / NR transformation laws, one-point matrices, particular two-point kernels, worldsheet vertex | 49 |
| NRH03 | SM 1: Γ² action density and flux on both saddles, cutoff value, on-shell value | 13 |
| NRH04 | SM 2: exact radial branches (Hill data, χ equation, exact hair profile), W₀ as a B-field shift, constant-L radial flow | 20 |
| NRH05 | SM 3.1–3.4: projectors and mixed fluctuation, boundary frames and sources, exact linearized EDFE (normal form, Cauchy constraints, constraint propagation, Taylor recursion), the component operators E^(0), E^(1),NR, ΔE^(1), the coupled hierarchy solved order by order, logarithmic coefficients, stress and type constraints | 57 |
| NRH06 | SM 3.5–3.8: Γ² variation and momenta, dictionary, second variation, background and linearized momenta, one-points, state fluctuations, ordered cutoff bilinear, derivative counterterms, finite bilinear and rows, contact rows, Lorentzian convention, Ward operators, nonchiral completion, kernels, conditional R completion | 67 |
| NRH07 | SM 3.9: covariant phase-space charges, Virasoro cocycle, NR charge cancellation and algebra | 29 |
| NRH08 | SM 4: first-order worldsheet, SNC clocks, Gomis–Ooguri limit, long-string energy, radial vertex operator, BRST weights and fusion | 32 |
| NRH09 | SM 5: ten-dimensional uplift and flux, Clifford algebra, vacuum and Riemannian Killing spinors, one-sided non-Riemannian jet system | 55 |
| NRH10 | SM 6: boundary candidate action and its bosonic and Grassmann symmetries | 14 |

[EQUATION_LEDGER.md](EQUATION_LEDGER.md) states the coverage of every numbered display;
[MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) pins the source hash and the label-to-number map;
[REFERENCE_RUN.md](REFERENCE_RUN.md) records the execution.

## Methods

The checks implement the computations as the manuscript states them.

- Backgrounds are the exact Bañados and everywhere non-Riemannian saddles with arbitrary
  chiral L±(x±) and, on the NR branch, arbitrary W₁(x⁺,x⁻); W₂ = −(l²/4)L₊′L₋′ is taken from
  the exact hair profile. Nothing is evaluated on sample data.
- The linearized EDFE are obtained by direct linearization of the exact curvature on each
  saddle, projected on the saddle frames (SM 3.3). The near-boundary solutions follow from the
  coupled hierarchy (SM 3.4), solved sequentially in the log-branch coefficients; the responses
  R±, H_s and the zero mode c_s are kept as free functions, and the stress and type constraints
  are read off from the remaining equations.
- Two-point functions are computed from the ordered second variation of the on-shell action
  (SM 3.8) with two independent bulk solutions inserted in the two slots. They are never
  obtained by differentiating one-point functions. The counterterm test integrates by parts to
  the canonical form Σ j₁,I q_I[j₂] and imposes the stress constraints in both slots, so
  "modulo total tangential derivatives" is explicit; a negative control confirms that the
  divergences remain without the counterterm.
- The nonlocal kernels follow from the Ward operators acting on the manuscript's Lorentzian
  inverse derivative ∂∓⁻¹δ² → 1/(2πiΔ±), with the time-ordering factor 1/i.
- The Killing-spinor jet system on the one-sided hairy background is derived from the printed
  frame and the semi-covariant connection, not transcribed.

## Conventions and symbols

The doubled order is `(dual x+, dual x-, dual y; x+, x-, y)`; dual derivatives vanish.
`JJ` is the O(3,3) metric. Antisymmetrization has unit weight. The connection obeys
`Gamma^B_BA = -2 partial_A d`; compatibility, trace and torsion identities are checked on
both saddles. `eta3` is the null-frame metric and `etab3 = -eta3`.

| Code | Manuscript quantity |
|---|---|
| `u` | e^{2y/l} in the Letter, SM 1 and SM 2 files (`xs` carries d_y = (2u/l) d_u) |
| `z`, `yy` | e^{−2y/l} (the manuscript's u) and the explicit radial coordinate y in the SM 3 files; `DyZY` is d_y = d_yy − (2z/l) d_z |
| `Lp[xp]`, `Lm[xm]`, `W0`, `W1[xp,xm]` | chiral data and hair modes; `psip`, `psim` are ψ± = L±^{−1/2} (`NRLpsi` maps L± to ψ±) |
| `RiemannianSaddleExact[]`, `NonRiemannianSaddleExact[W]` | exact saddles with frames V, V̄, dilaton and (R) metric and B field |
| `SaddleSeries[sd, n]` | the saddle truncated to order z^n (exact: nothing lowers the z order) |
| `hpp, hpm, hmp, hmm, dd` | h_{⊕⊕̄}, h_{⊕⊖̄}, h_{⊖⊕̄}, h_{⊖⊖̄} and δd as functions of (xp, xm, yy) |
| `a, b, r, c, v` (slot k: `a1`, `a2`, …) | the sources (h^(0)_{⊖⊖̄}, h^(0)_{⊕⊕̄}, h^(0)_{⊖⊕̄}, h^(0)_{⊕⊖̄}, δd^(0)) of SM (105) |
| `Rp, Rm, H, cs` | the responses h^(2)_{⊕⊕̄}, h^(2)_{⊖⊖̄}, h^(2)_{⊕⊖̄} and the zero mode c_s |
| `LinearizedEDFEComponents` | E_{pq̄}, E₀ of SM (41) through the requested z order |
| `MomentumProjected` | A^y_{pq̄} and the same-chirality projections of SM (87b) |
| `GenLieH`, `GenLieD` | generalized Lie derivative of H and d |
| `chiq`, `ch`, `esig` | χ = 2√2 arctanh q, χ as a symbol, e^σ = ψ₋/ψ₊ |

`T = exp(chi/(2 sqrt(2)))` simplifies hyperbolic identities on the exact NR saddle.
Constant-L symbols denote explicit constant-profile sectors.

## Scope

The general variation formulas of SM 3.5 are not separately checked (the linearized EDFE are
computed directly), the symmetry of the action-response matrix is a stated property, and the
ten-dimensional Killing-spinor statements marked **S** in the ledger are cited. The undetermined
interior responses, the hair self-response and the finite local terms induced by the cutoff
subtraction are carried as free data, as in the manuscript; the suite does not determine complete
Green functions. A passing reduced jet system is not a construction of nonzero supercharges.
