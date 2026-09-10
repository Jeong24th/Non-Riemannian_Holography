# Mathematica verification

Exact symbolic checks for *Non-Riemannian Hair in Long-String Holography*.
The `.wl` files are canonical; `.nb` files contain the same input expressions.

```sh
wolframscript -file NRH00_RunAll.wl
```

Or keep this directory together, open `NRH00_RunAll.nb`, and choose
**Evaluation > Evaluate Notebook**. Section files also run independently.
A failed check returns a nonzero command-line exit status.

| File | Checks |
|---|---|
| NRH00 | runner |
| NRH01 | DFT connection, curvature, action, momenta, shared backgrounds |
| NRH02 | Riemannian saddle, boundary geometry, Ward laws, one-points |
| NRH03 | non-Riemannian saddle, radial equation, hair, asymptotic laws |
| NRH04 | variations, linearized equations, response normalization, particular two-point kernels |
| NRH05 | covariant surface charges, C-bracket, cocycle, on-shell action |
| NRH06 | first-order worldsheet, GO limit, radial weights and marginality |
| NRH07 | ten-dimensional bosonic uplifts, frames, Clifford algebra, Killing-spinor reductions |
| NRH08 | boundary candidate action and its bosonic/Grassmann symmetries |

[EQUATION_LEDGER.md](EQUATION_LEDGER.md) states the coverage of each current equation.
[MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) pins the source and equation numbers.
[REFERENCE_RUN.md](REFERENCE_RUN.md) records actual execution results.

## Conventions and symbols

The doubled order is `(dual x+, dual x-, dual y; x+, x-, y)`; dual derivatives vanish.
`JJ` is the O(3,3) metric. Antisymmetrization has unit weight. The connection obeys
`Gamma^B_BA = -2 partial_A d`; compatibility and curvature identities are checked.
`eta3` is the null-frame metric and `etab3 = -eta3`.

| Code | Manuscript quantity |
|---|---|
| `gR`, `bR`, `HR`, `dR` | Riemannian metric, B field, generalized metric and dilaton |
| `HNR`, `dNR` | non-Riemannian generalized metric and dilaton |
| `Hinf`, `Vinf`, `Vbinf` | common limiting generalized metric and double frame |
| `Lp[xp]`, `Lm[xm]`, `W0[xp,xm]`, `W1[xp,xm]` | chiral data and hair modes |
| `psip`, `psim` | manuscript psi = L^(-1/2); no extra LpPsi/LpP aliases |
| `u`, `z` | exp(2y/l), exp(-2y/l) = 1/u |
| `ch`, `chU`, `esig` | chi, chi expressed in u, exp(sigma) |
| `hmat`, `hpp`, `hpm`, `hmp`, `hmm` | flat tangential fluctuation and its components |
| `MomentumCore`, `MomentumAK` | unprojected mathcal A, its mixed doubled projection |
| `aR`, `aNR`, `mR`, `mNR` | particular kernels A_R, A_NR, M_R, M_NR |
| `centralCharge` | c = 3l/(2G) |

`RiemannianMetric/B/D`, `NonRiemannianH`, and `NRBoundaryH/D` share background
formulas across files. Independent frame reconstructions and curvature contractions
remain separate checks. Constant L symbols denote explicit constant-profile sectors;
10-dimensional uplift matrices and radial expansions differ from the exact
three-dimensional fields. `T = exp(chi/(2 sqrt(2)))` simplifies hyperbolic identities.
State data are never replaced by numerical sample values.

## Response scope

NRH04 differentiates both bulk frames in the momentum while holding the independent
variation argument h1 fixed. It checks the vacuum and both finite-data backgrounds
through exp(-2y/l), logarithmic finite coefficients, and Ward-derived particular
kernels A/M. Background coefficients are evaluated at the first point. The code
uses M directly, without additional B_R/C_NR aliases.

The 5x5 matrix is an organizational definition. Its Ward-undetermined remainder,
including hair self-response, interior/state/zero-mode prescriptions and local contacts,
is not computed. These checks do not determine complete Green functions or finite local
counterterms. The ledger also states the restricted scope of the SUSY checks.
