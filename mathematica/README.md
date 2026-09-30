# Mathematica verification

Symbolic checks for *Non-Riemannian Holography: Long Strings and Soft Hair*.
The `.wl` files are canonical; each `.nb` notebook contains the same input
expressions. The files follow the manuscript's broad organization, with SM3
split into linearized dynamics, renormalization, and charges.

## Version and execution

Current documentation targets source SHA-256
`57960AAD42F2A9150C0958D1C253F08B0939B76C0E646EE250A7D824D6319B68`
(2026-09-30): Letter (1)-(21), 23 numbered displays; SM1-SM129, 130 displays;
153 total. SM3 now has six subsections and 49 displays.

The recorded Mathematica 13.2.1 run passed 336/336 checks on 2026-09-17 against
an earlier source. The current update preserves executable algebra while
updating documentation, comments, and check descriptions. No Wolfram runtime
was available for a new run. Historical success does not imply verification
of every revised or newly added formula.

From this directory, run:

```sh
wolframscript -file NRH00_RunAll.wl
```

Or keep the directory together, open `NRH00_RunAll.nb`, and choose
**Evaluation > Evaluate Notebook**. Section files run independently; a failed
check returns a nonzero command-line exit status. The recorded full run took
303.6 seconds. See [REFERENCE_RUN.md](REFERENCE_RUN.md).

| File | Current correspondence and retained calculations | Checks in recorded run |
|---|---|---:|
| NRH00 | Runner | -- |
| NRH01 | Shared tools: O(3,3) metric, connection, curvatures, Gamma-squared density and flux, momenta, generalized Lie derivative, saddle series, frame variation, projected linearized EDFE | -- |
| NRH02 | Letter (1)-(21): saddle and boundary data, field equations, boundary Ward components, symmetry transformations, one-point matrices, particular kernels, worldsheet vertex; no general off-shell Codazzi check | 49 |
| NRH03 | SM1: action density, radial flux, regulated and on-shell saddle values | 13 |
| NRH04 | SM2: Hill data, radial equation and hair profile, W0 shift, constant-profile radial flow | 20 |
| NRH05 | SM3.1 and SM3.3-3.5: fluctuation algebra, frames, direct curvature linearization and near-boundary hierarchy; retains expanded component operators, constraints and auxiliary recursion checks | 57 |
| NRH06 | SM3.2-3.5: action variation, momenta, ordered asymptotic bilinear, quadratic counterterms, finite response rows and particular kernels; retains the earlier expanded calculation | 67 |
| NRH07 | SM3.6: surface-charge one-form, R and NR charge components, cocycles and charge algebra | 29 |
| NRH08 | SM4: first-order worldsheet, SNC clocks, Gomis-Ooguri limit, static long-string energy, radial vertex, linearized weights and fusion | 32 |
| NRH09 | SM5: specified uplift sectors, Clifford algebra, vacuum and R Killing-spinor ingredients, reduced one-sided NR jet system | 55 |
| NRH10 | SM6: classical candidate action, bosonic non-Abelian tests and specified Grassmann symmetry tests | 14 |

[EQUATION_LEDGER.md](EQUATION_LEDGER.md) maps every current numbered display to
its actual public coverage. [MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) pins the source
and label-to-number map. [SM_REVISION_NOTES.md](SM_REVISION_NOTES.md) explains
label migration. A legacy label in a check name may refer to a retained
intermediate calculation whose display was removed from the compact manuscript;
it is not automatically a current equation or a claim of complete verification.

## Methods

- General-background calculations retain arbitrary chiral L-plus/L-minus and
  arbitrary NR hair W1. Constant-profile, vacuum, and one-sided restrictions
  are stated separately. The derivative-dependent W2 coefficient follows from
  the exact hair profile.
- The linearized EDFE are obtained by direct linearization of the scalar and
  projected curvature, followed by projection on the saddle frames. The
  coupled near-boundary hierarchy determines the logarithmic coefficients and
  stress/type constraints while retaining the undetermined responses.
- The two-point calculation uses the ordered second variation with two
  independent bulk solutions and extracts source linear response with the
  manuscript normalization. It is the retained expanded asymptotic calculation;
  the current compact metric Hessian has no direct public equality check.
- Counterterm cancellation is tested after integration by parts, with the
  stress constraints imposed in both variation slots. A negative control
  confirms that the divergences remain without the counterterm. This does
  not establish a universal nonlinear intrinsic counterterm completion.
- Particular nonlocal kernels follow from Ward operators acting on the joint
  Lorentzian inverse derivative, with the time-ordering factor 1/i. This does
  not determine all homogeneous responses or a complete Green function.
- The reduced one-sided Killing-spinor jet system is derived from its frame
  and semi-covariant connection. Separate Clifford-algebra checks do not make
  it a verification of every equation in the current real ten-dimensional basis.

## Conventions and symbols

The doubled order is `(dual x+, dual x-, dual y; x+, x-, y)`; dual derivatives
vanish. `JJ` is the O(3,3) metric. Antisymmetrization has unit weight. The
connection obeys `Gamma^B_BA = -2 partial_A d`; compatibility, trace and torsion
identities are checked on the stated saddle series. `eta3` is the null-frame
metric and `etab3 = -eta3`.

| Code | Manuscript quantity |
|---|---|
| `u` | e^(2y/l) in NRH02-NRH04; the radial derivative is `(2u/l) D[...,u]` |
| `z`, `yy` | e^(-2y/l), the current manuscript's u, and the explicit radial coordinate y in SM3 calculations; `DyZY = d_yy - (2z/l) d_z` |
| `Lp[xp]`, `Lm[xm]`, `W0`, `W1[xp,xm]` | Chiral data and hair modes; `psip`, `psim` are psi-plus/minus = L-plus/minus^(-1/2) |
| `RiemannianSaddleExact[]`, `NonRiemannianSaddleExact[W]` | Exact backgrounds with double vielbeins and dilaton, plus the R metric and B field |
| `SaddleSeries[sd,n]` | Series truncation to order z^n; the radial derivative does not lower the z order |
| `hpp,hpm,hmp,hmm,dd` | h_(plus,bar-plus), h_(plus,bar-minus), h_(minus,bar-plus), h_(minus,bar-minus), delta d |
| `a,b,r,c,v` (slot k: `a1,a2,...`) | Leading sources h_(minus,bar-minus), h_(plus,bar-plus), h_(minus,bar-plus), h_(plus,bar-minus), delta d, in the order of Letter (19) |
| `r` / `r0` | Type-changing source h^(0)_(minus,bar-plus); `c` / `c0` is the fourth, W0-source channel |
| `Rp,Rm,H,cs` | Full stress responses h^(2)_(plus,bar-plus), h^(2)_(minus,bar-minus), hair response H_s, and constant mode c_s; on R, c_s = -4 f0 |
| `LinearizedEDFEComponents` | Direct projected curvature variations through the requested z order; not a separate implementation of the universal tensor Box |
| `MomentumProjected` | Mixed and same-chirality projections of radial A^y |
| `GenLieH`, `GenLieD` | Generalized Lie derivatives of H and d |
| `chiq`, `ch`, `esig` | chi = 2 sqrt(2) arctanh(q), symbolic chi, and exp(sigma) |
| `GG[ch]`, `Gp[ch]` in NRH02/NRH04 | Manuscript radial calligraphic I(chi) and its first derivative, respectively; labels `MainGprofile` and `NRGprofile` are unchanged |
| Module-local `GG[c]` in NRH05 | The first derivative of calligraphic I; its series is integrated to reconstruct the profile |

The radial helper identifiers are retained. Bare `I` in Wolfram Language is the
protected imaginary unit and is not used as a replacement profile name. Newton's
`G`, the Einstein tensor and the two-point matrix retain their distinct meanings.
`T = exp(chi/(2 sqrt(2)))` simplifies hyperbolic identities on the NR saddle.

## Limits of current coverage

The ledger distinguishes direct checks, partial checks, definitions, cited
results and uncovered formulas. In particular:

- Letter (13), the general off-shell Codazzi identity in FG gauge, is not
  checked by evaluating exact saddles or by the constant-boundary Ward test.
- The universal tensor Box includes derivative, curvature and connection terms
  whose combination is fully covariant. Direct component linearization does
  not verify the full Box identities or off-shell curvature-variation formulas.
- SM34's exact radial metric variation and SM35's compact polarized Hessian
  are not directly tested. Tests retaining the `SMsecondvariation` label
  concern coset-coordinate auxiliaries and the earlier asymptotic bilinear.
- The fixed-flux K3 solution and interior H_R matching in SM43-SM44 are not
  computed by the suite. Old conditional normalization tests do not establish
  the specified-vacuum hair correlator in SM48. General H_s remains free.
- The new compact a/b/omega and integrated NR coefficient expressions have
  supporting checks of the earlier component hierarchy and differential
  constraints, but no direct public bridge for every rewritten display.
- Worldsheet checks have their stated radial, near-boundary or constant-profile
  scope. A longitudinal contraction check does not establish nonsingularity
  or exact marginality of the dressed W1 vertex.
- The uplift and spinor tests cover the specified sectors, including the reduced
  one-sided NR system. They do not construct nonzero fermionic charges or prove
  all current real-basis and global counting claims. Grassmann candidate-action
  checks use the stated Abelian/zero-connection realization; separate bosonic
  tests retain non-Abelian covariant derivatives.

The boundary candidate is classical. Quantum conformal invariance, anomalies,
and a complete holographic operator identification are not computed here.
