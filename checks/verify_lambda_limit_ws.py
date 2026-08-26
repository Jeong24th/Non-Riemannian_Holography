#!/usr/bin/env python3
"""Hard-fail regression test: lambda-limit worldsheet reduction and the
(SMdeltaL) sign dictionary (HCKim notes Sec. 6 vs NR_Holography.tex SM).

Provenance: session-agent verification (2026-07-23) of HCKim notes Sec. 6,
shared and hardened 2026-07-24 after the Codex cross-review.  Codex's
independent transformation chain (identity map fails; d_K = db_M forced;
A_K = -A_M) reproduces the same conclusions as checks 9-11 below.

Guarded facts (conventions of the notes: (d, db) = (d_sigma, d_sigmabar),
sigma-model L = (1/2 pi alpha')(g+B)_{mu nu} dx^mu dbx^nu, Banados family
with B = (e^{2y/l}+L+L- e^{-2y/l}) dx^- ^ dx^+):

  1) (6.10) is exactly the (g+B) sigma model                      -> True
  2) the (g-B) convention does NOT reproduce (6.10)               -> False
  3) (6.11): y -> y + l*lam, L+- -> e^lam L+- scaling             -> True
  4) (6.12): J Jbar rewriting identity                            -> True
  5) (6.14): integrating out b, bbar reproduces -2 e^{2lam} J Jb  -> True
  6) (6.15): lam -> infty limit of L_aux                          -> True
  7) (6.16): assembled limit Lagrangian                           -> True
  8) (6.17)|_{W1 = 4 L+ L-} == (6.16)                             -> True
  9) orientation-flipped (6.17) == paper L0+DeltaL with
     Wcal = -W1 e^{-2y/l}                                         -> False
 10) orientation-flipped (6.17) == paper L0+DeltaL with
     Wcal = +W1 e^{-2y/l}                                         -> True
 11) unflipped: his DeltaL = -(paper DeltaL)                      -> True

Checks 9/10 are the sign negative/positive controls: under the UNIQUE
dictionary (target labels x^+- are bulk-pinned, so only the worldsheet
orientation flip d <-> db relates the two constraint displays; beta =
e^{y/l} bbar, betabar = e^{y/l} b), the relative sign of (SMdeltaL) is
"+".  The ABSOLUTE worldsheet sign remains withheld (the manuscript does
not fix its orientation), and finite-W1 exact marginality is withheld.

Exit codes: 0 = ALL CHECKS PASSED; 1 = at least one check failed;
2 = missing dependency.  Usage:
    python checks/verify_lambda_limit_ws.py [--selftest-fail] [--strict-pin]
(--selftest-fail deliberately flips the expectation of check 10 to
demonstrate that the harness fails hard; it must exit 1.  --strict-pin
turns a SymPy version different from the pin into a dependency error,
exit 2, instead of the default warning.)

A stdlib-only core of the decisive sign checks (9-11) is provided in
checks/verify_lambda_limit_sign_core.py for Pythons without sympy.
"""

from __future__ import annotations

import sys

PINNED_SYMPY = "1.14.0"


def _configure_streams() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass


_configure_streams()

try:
    import sympy as sp
except ImportError:
    print(f"[DEPENDENCY ERROR] SymPy {PINNED_SYMPY} is required. "
          "Install per requirements-verification.txt "
          "(isolated environment only).", file=sys.stderr)
    sys.exit(2)

FAILURES: list[str] = []


def check(name: str, actual, expected) -> None:
    ok = (actual == expected)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: actual={actual} expected={expected}")
    if not ok:
        FAILURES.append(name)


def main(selftest_fail: bool = False, strict_pin: bool = False) -> None:
    print(f"sys.executable : {sys.executable}")
    print(f"python version : {sys.version.split()[0]}")
    print(f"sympy version  : {sp.__version__}"
          + ("" if sp.__version__ == PINNED_SYMPY
             else f"  [WARNING: pinned version is {PINNED_SYMPY}]"))
    if strict_pin and sp.__version__ != PINNED_SYMPY:
        print(f"[DEPENDENCY ERROR] SymPy {PINNED_SYMPY} is required "
              f"(--strict-pin; found {sp.__version__}). Install per "
              "requirements-verification.txt "
              "(isolated environment only).", file=sys.stderr)
        sys.exit(2)
    if selftest_fail:
        print("[SELFTEST] deliberately flipping the expectation of check 10")

    l, t = sp.symbols("l t", positive=True)          # t = e^{-lambda}
    y = sp.symbols("y", real=True)
    Lp, Lm = sp.symbols("L_p L_m", positive=True)
    E = sp.exp(y / l)
    dxp, dxm, dy = sp.symbols("dx_p dx_m dy")        # d x^+, d x^-, d y
    bxp, bxm, by = sp.symbols("bx_p bx_m by")        # db x^+, db x^-, db y
    dX = sp.Matrix([dxp, dxm, dy])
    bX = sp.Matrix([bxp, bxm, by])

    b = E**2 + Lp * Lm * E**-2
    g = sp.Matrix([[2 * Lp, -b, 0], [-b, 2 * Lm, 0], [0, 0, 1]])
    B = sp.Matrix([[0, -b, 0], [b, 0, 0], [0, 0, 0]])   # B_{-+} = +b

    L_gplusB = sp.expand((dX.T * (g + B) * bX)[0, 0])
    L_gminusB = sp.expand((dX.T * (g - B) * bX)[0, 0])
    L610 = sp.expand(dy * by + 2 * Lp * dxp * bxp + 2 * Lm * dxm * bxm
                     - 2 * (E**2 + Lm * Lp * E**-2) * dxp * bxm)
    check("(6.10) g+B identification", sp.simplify(L_gplusB - L610) == 0, True)
    check("(6.10) g-B identification", sp.simplify(L_gminusB - L610) == 0, False)

    L_scaled = L_gplusB.subs({E: E / t, Lp: Lp / t, Lm: Lm / t},
                             simultaneous=True)
    L611 = sp.expand(dy * by + (2 * Lp / t) * dxp * bxp
                     + (2 * Lm / t) * dxm * bxm
                     - 2 * (E**2 / t**2 + Lm * Lp * E**-2) * dxp * bxm)
    check("(6.11) lambda-scaled Lagrangian",
          sp.simplify(sp.expand(L_scaled) - L611) == 0, True)

    Jc = E * dxp - t * E**-1 * Lm * dxm
    Jbc = E * bxm - t * E**-1 * Lp * bxp
    L612 = sp.expand(dy * by - 2 * Jc * Jbc / t**2
                     - 2 * E**-2 * Lp * Lm * (dxp * bxm - dxm * bxp))
    check("(6.12) J Jbar identity", sp.simplify(L612 - L611) == 0, True)

    bb, bbb = sp.symbols("b bbar")
    Laux = bbb * Jc + bb * Jbc + sp.Rational(1, 2) * t**2 * bb * bbb
    solb = sp.solve(sp.diff(Laux, bbb), bb)[0]
    solbb = sp.solve(sp.diff(Laux, bb), bbb)[0]
    Laux_os = sp.simplify(Laux.subs({bb: solb, bbb: solbb}))
    check("(6.14) auxiliary-field on-shell identity",
          sp.simplify(Laux_os + 2 * Jc * Jbc / t**2) == 0, True)

    Laux_lim = sp.limit(Laux, t, 0)
    check("(6.15) lambda->infty multiplier limit",
          sp.simplify(Laux_lim - (E * bbb * dxp + E * bb * bxm)) == 0, True)

    L616 = sp.expand(dy * by + E * bbb * dxp + E * bb * bxm
                     - 2 * E**-2 * Lp * Lm * (dxp * bxm - dxm * bxp))
    Lfull_lim = sp.limit(sp.expand(dy * by + Laux
                                   - 2 * E**-2 * Lp * Lm
                                   * (dxp * bxm - dxm * bxp)), t, 0)
    check("(6.16) assembled limit Lagrangian",
          sp.simplify(Lfull_lim - L616) == 0, True)

    W1 = sp.symbols("W_1")
    L617 = sp.expand(dy * by + E * bbb * dxp + E * bb * bxm
                     - sp.Rational(1, 2) * E**-2 * W1 * (dxp * bxm - dxm * bxp))
    check("(6.17) W1 = 4 L+ L- substitution",
          sp.simplify(L617.subs(W1, 4 * Lp * Lm) - L616) == 0, True)

    beta, brbeta, cW = sp.symbols("beta brbeta cW")
    Lpaper = sp.expand(dy * by + beta * bxp + brbeta * dxm
                       + sp.Rational(1, 2) * cW * (dxp * bxm - bxp * dxm))
    flip = {dxp: bxp, bxp: dxp, dxm: bxm, bxm: dxm, dy: by, by: dy}
    L617_flipped = L617.subs(flip, simultaneous=True)
    subs_minus = {beta: E * bbb, brbeta: E * bb, cW: -W1 * E**-2}
    subs_plus = {beta: E * bbb, brbeta: E * bb, cW: +W1 * E**-2}
    check("flipped, Wcal = -W1 e^{-2y/l} identification",
          sp.simplify(L617_flipped - Lpaper.subs(subs_minus)) == 0, False)
    expected10 = False if selftest_fail else True
    check("flipped, Wcal = +W1 e^{-2y/l} identification",
          sp.simplify(L617_flipped - Lpaper.subs(subs_plus)) == 0, expected10)

    his_dL = -sp.Rational(1, 2) * E**-2 * W1 * (dxp * bxm - dxm * bxp)
    paper_dL = sp.Rational(1, 2) * (W1 * E**-2) * (dxp * bxm - bxp * dxm)
    check("unflipped: his DeltaL = -(paper DeltaL)",
          sp.simplify(his_dL + paper_dL) == 0, True)

    if FAILURES:
        print("FAILED CHECKS: " + "; ".join(FAILURES))
        sys.exit(1)
    print("ALL CHECKS PASSED")
    sys.exit(0)


if __name__ == "__main__":
    main(selftest_fail="--selftest-fail" in sys.argv[1:],
         strict_pin="--strict-pin" in sys.argv[1:])

