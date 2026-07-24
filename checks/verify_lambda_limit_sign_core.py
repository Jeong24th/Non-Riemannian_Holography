#!/usr/bin/env python3
"""Stdlib-only core negative/positive controls for the (SMdeltaL) sign
dictionary (no sympy required).

This isolates the DECISIVE orientation/sign content of
calculations/verify_lambda_limit_ws.py (checks 9-11 there) in exact
integer polynomial arithmetic over commuting monomials, so it runs on
any Python >= 3.8 without third-party packages.

Setting: near-boundary (chi -> 0) limit.  Fields as commuting symbols
dx+, dx-, bx+, bx- (d = d_sigma, b = d_sigmabar), multipliers bb ("b"),
bbb ("bbar"), weight w standing for W1 e^{-2y/l} (opaque positive
coefficient), and e for e^{y/l}.

Kim (6.17), flat limit, in (1/2 pi alpha') units:
    L_K  =  e*bbb*dx+  +  e*bb*bx-  -  (w/2)*(dx+*bx- - dx-*bx+)
Paper (SMGO)+(SMdeltaL), same units, with Wcal = s*w (s = +-1):
    L_P(s) =  beta*bx+  +  betab*dx-  +  (s*w/2)*(dx+*bx- - bx+*dx-)

Guarded facts:
  1) identity map: multiplier monomial sets do not match           -> True
  2) orientation flip (d<->b) maps Kim's multiplier terms onto the
     paper's with beta = e*bbar, betabar = e*b                     -> True
  3) flipped Kim == paper with s = +1                              -> True
  4) flipped Kim == paper with s = -1                              -> False
  5) unflipped: Kim DeltaL = -(paper DeltaL with s = +1)           -> True

Exit codes: 0 = ALL CHECKS PASSED; 1 = failure.
Usage: python calculations/verify_lambda_limit_sign_core.py [--selftest-fail]
"""

from __future__ import annotations

import sys
from fractions import Fraction


def _configure_streams() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass


_configure_streams()

FAILURES: list[str] = []


def check(name: str, actual, expected) -> None:
    ok = (actual == expected)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: actual={actual} expected={expected}")
    if not ok:
        FAILURES.append(name)


# --- exact polynomials: dict {frozenset-of-(symbol,power) -> Fraction} ---
def mono(*syms: str):
    key: dict[str, int] = {}
    for s in syms:
        key[s] = key.get(s, 0) + 1
    return frozenset(key.items())


def poly(*terms):
    """terms: (coeff, sym, sym, ...)"""
    out: dict = {}
    for term in terms:
        c = Fraction(term[0])
        k = mono(*term[1:])
        out[k] = out.get(k, Fraction(0)) + c
        if out[k] == 0:
            del out[k]
    return out


def padd(p, q, cq=Fraction(1)):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, Fraction(0)) + cq * c
        if out[k] == 0:
            del out[k]
    return out


def pmap(p, symmap):
    """apply a symbol substitution (symbol -> symbol) to every monomial"""
    out: dict = {}
    for k, c in p.items():
        newk: dict[str, int] = {}
        for s, n in k:
            s2 = symmap.get(s, s)
            newk[s2] = newk.get(s2, 0) + n
        kk = frozenset(newk.items())
        out[kk] = out.get(kk, Fraction(0)) + c
        if out[kk] == 0:
            del out[kk]
    return out


def main(selftest_fail: bool = False) -> None:
    print(f"sys.executable : {sys.executable}")
    print(f"python version : {sys.version.split()[0]}")
    print("dependencies   : stdlib only")
    if selftest_fail:
        print("[SELFTEST] deliberately flipping the expectation of check 3")

    half = Fraction(1, 2)

    # Kim (6.17) flat limit
    L_K = poly((1, "e", "bbar", "dx+"), (1, "e", "b", "bx-"),
               (-half, "w", "dx+", "bx-"), (half, "w", "dx-", "bx+"))

    def L_P(s: int):
        return poly((1, "beta", "bx+"), (1, "betab", "dx-"),
                    (Fraction(s) * half, "w", "dx+", "bx-"),
                    (-Fraction(s) * half, "w", "bx+", "dx-"))

    # dictionary: beta = e*bbar, betab = e*b  (applied by renaming in L_P)
    dict_sub = {"beta": None, "betab": None}  # placeholder; done via pmap on L_K side

    # 1) identity map: Kim multiplier monomials {e*bbar*dx+, e*b*bx-} vs
    #    paper's {beta*bx+, betab*dx-} with beta=e*bbar, betab=e*b, i.e.
    #    {e*bbar*bx+, e*b*dx-}: no overlap without a flip.
    kim_mult = {mono("e", "bbar", "dx+"), mono("e", "b", "bx-")}
    paper_mult = {mono("e", "bbar", "bx+"), mono("e", "b", "dx-")}
    check("identity map fails (multiplier monomials disjoint)",
          len(kim_mult & paper_mult) == 0, True)

    # orientation flip d <-> b on coordinates
    flip = {"dx+": "bx+", "bx+": "dx+", "dx-": "bx-", "bx-": "dx-"}
    L_K_flip = pmap(L_K, flip)

    # 2) flipped multiplier terms match the paper's exactly
    flip_mult = {k for k in L_K_flip if any(s in ("b", "bbar") for s, _ in k)}
    check("flip maps multipliers onto beta=e*bbar, betab=e*b structure",
          flip_mult == paper_mult, True)

    # paper polynomials with beta, betab replaced by e*bbar, e*b
    def L_P_sub(s: int):
        return poly((1, "e", "bbar", "bx+"), (1, "e", "b", "dx-"),
                    (Fraction(s) * half, "w", "dx+", "bx-"),
                    (-Fraction(s) * half, "w", "bx+", "dx-"))

    expected3 = False if selftest_fail else True
    check("flipped Kim == paper with s = +1",
          padd(L_K_flip, L_P_sub(+1), Fraction(-1)) == {}, expected3)
    check("flipped Kim == paper with s = -1",
          padd(L_K_flip, L_P_sub(-1), Fraction(-1)) == {}, False)

    # 5) unflipped DeltaL comparison
    dL_K = poly((-half, "w", "dx+", "bx-"), (half, "w", "dx-", "bx+"))
    dL_P = poly((half, "w", "dx+", "bx-"), (-half, "w", "bx+", "dx-"))
    check("unflipped: Kim DeltaL = -(paper DeltaL)",
          padd(dL_K, dL_P) == {}, True)

    if FAILURES:
        print("FAILED CHECKS: " + "; ".join(FAILURES))
        sys.exit(1)
    print("ALL CHECKS PASSED")
    sys.exit(0)


if __name__ == "__main__":
    main(selftest_fail="--selftest-fail" in sys.argv[1:])

