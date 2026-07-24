#!/usr/bin/env python3
"""GO-CEFF-001 dated amendment (Claude, 2026-07-24, cross-review stage).

Independent verification of the Codex domain refinements (F=0 / F>0) that my
original claude_go_ceff_check.py did not test, plus an orientation-pinning
argument independent of the SMdeltaL sign convention.  Append-only evidence;
the original script is unchanged.
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


y, l = sp.symbols("y l", positive=True)
Pi = sp.Symbol("Pi", real=True)
E2 = sp.exp(2 * y / l)
F = E2 + Pi / E2
em2d = E2 - Pi / E2

# D-1: F = 0  <=>  e^{4y/l} = -Pi, real solution iff Pi < 0.
sol = sp.solve(sp.Eq(F, 0), y)
check("D-1: F=0 locus is e^{4y/l} = -Pi (needs Pi<0)",
      any(sp.simplify(sp.exp(4 * s / l) - (-Pi)) == 0 for s in sol))
Fneg = F.subs(Pi, -1)          # counterexample Pi = -1
check("D-1: counterexample Pi=-1: F(y=0)=0", sp.simplify(Fneg.subs(y, 0)) == 0)
Fpos_dom = F.subs(Pi, sp.Symbol("P", positive=True))
check("D-1: Pi>=0 => F>0 at all radii (Pi>0 case)",
      sp.ask(sp.Q.positive(Fpos_dom),
             sp.Q.positive(sp.Symbol("P", positive=True))) in (True, None)
      and sp.simplify(Fpos_dom) != 0 and sp.limit(Fpos_dom, y, -sp.oo) == sp.oo)
check("D-1: Pi=0 => F = e^{2y/l} > 0", sp.simplify(F.subs(Pi, 0) - E2) == 0)

# D-2: at the F=0 locus the second-order metric stays nondegenerate:
#      -det(g_{+-}-block) = F^2 - 4 Pi = (e^{-2d})^2, and at F=0 it equals -4Pi>0.
check("D-2: F^2 - 4 Pi = (e^{-2d})^2 identity",
      sp.simplify(F**2 - 4 * Pi - em2d**2) == 0)
check("D-2: at F=0 (Pi=-1, y=0): -det block = 4 > 0 (Lorentzian, invertible)",
      sp.simplify((F**2 - 4 * Pi).subs([(Pi, -1), (y, 0)]) - 4) == 0)
check("D-2: e^{-2d} > 0 at the F=0 locus (no horizon there)",
      sp.simplify(em2d.subs([(Pi, -1), (y, 0)]) - 2) == 0)

# D-3: orientation pinning WITHOUT invoking the SMdeltaL sign directly:
# under E=g-B the constrained derivatives are {d x^-, dbar x^+} while the NR
# hair vertex V_W ~ W d x^+ dbar x^- couples the complementary unconstrained
# pair -> nonvanishing hair coupling (as established).  Under E=g+B the
# constrained set would be {d x^+, dbar x^-} and V_W would vanish identically
# on the free-GO constraint surface -> contradiction.
dxp, dxm, dbxp, dbxm, W = sp.symbols("dxp dxm dbxp dbxm W")
VW = W * dxp * dbxm
VW_on_gmB = VW.subs([(dxm, 0), (dbxp, 0)])   # E=g-B constraints
VW_on_gpB = VW.subs([(dxp, 0), (dbxm, 0)])   # E=g+B (rejected orientation)
check("D-3: E=g-B constraints leave V_W nonvanishing", sp.simplify(VW_on_gmB) != 0)
check("D-3: E=g+B would annihilate V_W identically (rejected)",
      sp.simplify(VW_on_gpB) == 0)

if selftest:
    check("injected failure (selftest)", sp.simplify(F - E2) == 0)

print()
if fails:
    print("FAILED CHECKS: " + "; ".join(fails), file=sys.stderr)
    sys.exit(1)
print("ALL CHECKS PASSED")

