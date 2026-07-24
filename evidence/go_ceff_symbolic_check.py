#!/usr/bin/env python3
"""GO-CEFF-001 independent verification (Claude), hard-fail + negative controls.

Conventions (manuscript `Rfields`, baseline 241E1754...77A5, lines 679-688):
  theta^+ = e^{y/l} dx^+ - e^{-y/l} L_-(x^-) dx^-
  theta^- = e^{y/l} dx^- - e^{-y/l} L_+(x^+) dx^+
  ds^2    = dy^2 - 2 theta^+ theta^-          (symmetrized product)
  B_(2)   = F dx^- ^ dx^+ ,  F := e^{2y/l} + e^{-2y/l} L_+ L_-,  phi = phi_0
Coordinate order (y, x^+, x^-) = indices (0, 1, 2); B = (1/2) B_{mu nu} dx^mu ^ dx^nu
so B[2,1] = F = -B[1,2].

Claims verified:
  C-a  g_{+-} = -F exactly  (so metric and flux are aligned component-wise);
  C-b  exactly ONE chirality pairing cancels: (g+B)[2,1] = 0 with survivor
       (g+B)[1,2] = -2F, and mirrored for (g-B); which pairing survives is an
       orientation choice, pinned by matching the free system of (SMGO)
       [multipliers enforce dbar x^+ = 0 = d x^-]: the big term must sit on
       dbar x^+ d x^-.
  C-c  exact first-order rewriting at EVERY radius (completion of square,
       unit Jacobian):
         -2F dbar x^+ dx^- = beta dbar x^+ + betabar dx^- + (1/(2F)) beta betabar
       at stationarity;  identically
         beta dbar x^+ + betabar dx^- + (1/(2F)) beta betabar
           = (1/(2F)) (beta + 2F dx^-)(betabar + 2F dbar x^+) - 2F dbar x^+ dx^-.
  C-d  c_eff^2 := 2F: boundary limit 1/(2F) -> (1/2) e^{-2y/l} -> 0 (GO point,
       no tuning); vacuum L+L- = 0 gives 1/(2F) = (1/2) e^{-2y/l} exactly;
       at the constant-L horizon e^{4y/l} = L+L-: 2F = 4 sqrt(L+L-) != 0
       (c_eff finite and nonzero at the horizon; only e^{-2d} = E^2 - Pi/E^2
       vanishes there).
  C-e  near the boundary the Banados couplings 2L_+ dx^+ dbar x^+ +
       2L_- dx^- dbar x^- multiply derivatives that vanish on the free-GO
       constraint surface (dbar x^+ = 0 = d x^-) - classical decoupling,
       mirroring the NR-side chi->0 statement.

Negative controls:
  N-1  sign-flipped interaction -(1/(2F)) beta betabar breaks the round trip;
  N-2  B -> -B swaps which pairing cancels (orientation selectivity);
  N-3  claiming the interaction is exactly (1/2) e^{-2y/l} fails for L+L- != 0;
  N-4  the direction claim c_eff = e^{-y/l} implies interaction -> infinity at
       the boundary; the verified interaction -> 0.  Direction is c_eff = e^{+y/l}.
Run with --selftest-fail to verify the harness catches an injected failure.
"""
import sys
import sympy as sp

selftest = "--selftest-fail" in sys.argv
fails = []


def check(name, cond):
    ok = bool(cond)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        fails.append(name)


y, l, xp, xm, w = sp.symbols("y l x_p x_m w", positive=True)
Lp = sp.Function("L_p")(xp)
Lm = sp.Function("L_m")(xm)
E = sp.exp(y / l)
F = E**2 + Lp * Lm / E**2

# --- metric from ds^2 = dy^2 - 2 theta+ theta-  -----------------------------
tp = sp.Matrix([0, E, -Lm / E])       # theta^+ components on (y, x+, x-)
tm = sp.Matrix([0, -Lp / E, E])       # theta^- components
T = tp * tm.T
g = -(T + T.T)
g[0, 0] = 1

check("g_{++} = 2 L_+", sp.simplify(g[1, 1] - 2 * Lp) == 0)
check("g_{--} = 2 L_-", sp.simplify(g[2, 2] - 2 * Lm) == 0)
check("C-a: g_{+-} = -F exactly", sp.simplify(g[1, 2] + F) == 0)
check("no y-x mixing: g_{y+} = g_{y-} = 0", g[0, 1] == 0 and g[0, 2] == 0)

# --- B field ----------------------------------------------------------------
B = sp.zeros(3, 3)
B[2, 1] = F
B[1, 2] = -F

Mp = g + B     # combination multiplying  dx^mu dbar x^nu  (orientation s=+1)
Mm = g - B     # orientation s=-1

check("C-b: (g+B)[2,1] = 0 exactly (all y, arbitrary chiral L)",
      sp.simplify(Mp[2, 1]) == 0)
check("C-b: survivor (g+B)[1,2] = -2F", sp.simplify(Mp[1, 2] + 2 * F) == 0)
check("C-b mirror: (g-B)[1,2] = 0", sp.simplify(Mm[1, 2]) == 0)
check("C-b mirror: survivor (g-B)[2,1] = -2F", sp.simplify(Mm[2, 1] + 2 * F) == 0)
# Orientation pinning: in the manuscript's conventions the free limit must be
# (SMGO) with constraints dbar x^+ = 0 = d x^-, i.e. the multipliers pair as
# beta dbar x^+ + betabar d x^-; hence the BIG term must multiply dbar x^+ d x^-,
# which is the s = -1 slot [2,1] above.  Recorded, not derived here.

# --- C-c: exact first-order rewriting (completion of square) ----------------
beta, betab, dbxp, dxm = sp.symbols("beta betabar dbxp dxm")
L1 = beta * dbxp + betab * dxm + (1 / (2 * F)) * beta * betab
square = (1 / (2 * F)) * (beta + 2 * F * dxm) * (betab + 2 * F * dbxp) \
         - 2 * F * dbxp * dxm
check("C-c: completion-of-square identity (exact, every radius)",
      sp.simplify(sp.expand(L1 - square)) == 0)

sol = sp.solve([sp.diff(L1, beta), sp.diff(L1, betab)], [beta, betab], dict=True)[0]
L1_onshell = sp.simplify(L1.subs(sol))
check("C-c: round trip: on-shell L1 = -2F dbar x^+ dx^-",
      sp.simplify(L1_onshell + 2 * F * dbxp * dxm) == 0)

# --- C-d: c_eff^2 = 2F ------------------------------------------------------
coeff = sp.simplify(1 / (2 * F))
check("C-d: boundary limit of interaction is 0  (c_eff -> infinity)",
      sp.limit(coeff.subs([(Lp, sp.Symbol("a", positive=True)),
                           (Lm, sp.Symbol("b", positive=True))]), y, sp.oo) == 0)
check("C-d: vacuum (L+L- = 0) coefficient = (1/2) e^{-2y/l} exactly",
      sp.simplify(coeff.subs(Lp, 0) - sp.exp(-2 * y / l) / 2) == 0)
a, b = sp.symbols("a b", positive=True)
Fc = (E**2 + a * b / E**2)
hor = sp.solve(sp.Eq(E**4, a * b), y)  # e^{4y/l} = L+L-
c2_at_horizon = sp.simplify((2 * Fc).subs(y, sp.Rational(1, 4) * l * sp.log(a * b)))
check("C-d: c_eff^2 at the horizon = 4 sqrt(L+L-) != 0",
      sp.simplify(c2_at_horizon - 4 * sp.sqrt(a * b)) == 0)
em2d = E**2 - a * b / E**2
check("C-d: e^{-2d} DOES vanish at the horizon (contrast)",
      sp.simplify(em2d.subs(y, sp.Rational(1, 4) * l * sp.log(a * b))) == 0)

# --- C-e: Banados couplings multiply constrained derivatives ----------------
# The free-GO constraint surface is dbxp = 0 (and dxm = 0); the L couplings are
# 2 L_+ (dx^+)(dbar x^+) + 2 L_- (dx^-)(dbar x^-): each term contains one
# constrained factor (dbar x^+ or d x^-), so both vanish classically there.
dxp, dbxm = sp.symbols("dxp dbxm")
banados = 2 * Lp * dxp * dbxp + 2 * Lm * dxm * dbxm
check("C-e: Banados couplings vanish on the free-GO constraint surface",
      sp.simplify(banados.subs([(dbxp, 0), (dxm, 0)])) == 0)

# --- Negative controls ------------------------------------------------------
L1_bad = beta * dbxp + betab * dxm - (1 / (2 * F)) * beta * betab
solb = sp.solve([sp.diff(L1_bad, beta), sp.diff(L1_bad, betab)],
                [beta, betab], dict=True)[0]
bad_onshell = sp.simplify(L1_bad.subs(solb))
check("N-1: sign-flipped interaction FAILS the round trip (control)",
      sp.simplify(bad_onshell + 2 * F * dbxp * dxm) != 0)

Mp_flip = g - B  # B -> -B inside (g+B)
check("N-2: B -> -B swaps the cancelling pairing (control)",
      sp.simplify(Mp_flip[2, 1]) != 0 and sp.simplify(Mp_flip[1, 2]) == 0)

check("N-3: pure e^{-2y/l}/2 claim FAILS at finite y for L+L- != 0 (control)",
      sp.simplify(coeff - sp.exp(-2 * y / l) / 2) != 0)

wrong_dir = 1 / (2 * sp.exp(-2 * y / l))  # interaction if c_eff were e^{-y/l}
check("N-4: c_eff = e^{-y/l} would send the interaction to infinity (control)",
      sp.limit(wrong_dir, y, sp.oo) == sp.oo)

if selftest:
    check("injected failure (selftest)", sp.simplify(F - E**2) == 0)

print()
if fails:
    print("FAILED CHECKS: " + "; ".join(fails), file=sys.stderr)
    sys.exit(1)
print("ALL CHECKS PASSED")

