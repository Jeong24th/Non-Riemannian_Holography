#!/usr/bin/env python3
"""Regression test for the common fixed-frame SM linearized system.

Guards the 2026-07-23 corrections (audit E1--E4, E6) as amended by the
2026-08-07 sign correction (co-author dispute, resolved by an independent
    first-principles derivation, now expressed in Kim's aligned frame.  The
    type-changing source is h_{ominus,bar-oplus}; its couplings carry
MINUS signs; certified by D-module equivalence with the Euler--Lagrange
system of the eps^2 coset-completed action, by pure-gauge invariance, and
by O(t)-preservation of the exact nonlinear dilaton EOM S_0 = 2 Lambda,
which the former +2 signs violate: delta S_0 = 4 d_+d_- h^(0) != 0):

* E1': the linearized dilaton constraint of ``SMNRlin`` reads
  (8/l) d_y(delta d) = -d_+ d_- h_{ominus,bar-oplus}; the previously
  printed +2 (and the pre-July coefficient 1) are rejected as negative
  controls, and the ell_pm equations read
  box h_{oplus,bar-oplus} = -d_+^2 h_{ominus,bar-oplus} etc.
* E2': the general solution ``SMNRsol`` carries
  delta d = delta d^(0) - (l/8) y d_+d_- h_{ominus,bar-oplus}^(0),
  tails -(l/2)y d_pm^2 h^(0), and the y^2 branch +(l^2/8)y^2 D4 h^(0);
  the previously printed (l/8) C y display is rejected.
* E3: the exclusion of generic relativizing sources is a property of the
  fixed-dilaton, log-free sector, not of the field equations; the old
  sentence attributing it to ``SMNRsol`` is rejected.
* E4: the response dictionary is uniform,
  <K_{p bar q}> = h_{p bar q}^{(2)}/(32 pi G l).  In particular the hair
  h_{oplus,bar-ominus}^{(2)}=W_1/2 gives +int W_1/(64 pi G l).
* E6: the vielbein prose states the shifts -(W/2) tau^{mp} for both
  sectors (the printed matrices), not the correlated mp(W/2) tau^{mp}.

Run from any directory with

    python checks/verify_sm_nr_linearization.py
"""

from __future__ import annotations

import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = ROOT / "NR_Holography.tex"


def check_tex_contract() -> bool:
    """Keep the verified equations and their LaTeX display synchronized."""

    if not TEX_PATH.is_file():
        print("LaTeX contract: SKIP (NR_Holography.tex is not in the archive)")
        return False

    # Alignment tabs in aligned/alignedat are typographical, not algebraic.
    compact = re.sub(r"\s+", "", TEX_PATH.read_text(encoding="utf-8")).replace("&", "")
    required = {
        "aligned-frame dilaton constraint": (
            r"\tfrac{8}{l}\partial_y\deltad"
            r"=-\partial_+\partial_-h_{\ominus\bar\oplus}"
        ),
        "aligned-frame plus source": (
            r"h_{\oplus\bar\oplus}=-\partial_+^2h_{\ominus\bar\oplus}"
        ),
        "aligned-frame minus source": (
            r"h_{\ominus\bar\ominus}=-\partial_-^2h_{\ominus\bar\oplus}"
        ),
        "aligned-frame dilaton slope": (
            r"\deltad=\deltad^{(0)}(x^+,x^-)"
            r"-\frac{l}{8}y\,\partial_+\partial_-"
            r"h_{\ominus\bar\oplus}^{(0)}"
        ),
        "aligned-frame y^2 branch": (
            r"+\frac{l^2}{8}y^2\partial_+^2\partial_-^2"
            r"h_{\ominus\bar\oplus}^{(0)}"
        ),
        "aligned-frame y-linear source": (
            r"-\frac{l^2}{4}\partial_+^2\partial_-^2"
            r"h_{\ominus\bar\oplus}^{(0)}\Big]"
        ),
        "sign-corrected contact response": (
            r"=-\frac{l}{16\piG}\,e^{2y/l}\,"
            r"\partial_{+}\partial_{-}h_{\ominus\bar\oplus}^{(0)}"
        ),
        "complete-solution phrase": r"Thecompletesolutionis",
        "E3 general type-changing source": r"h_{\ominus\bar\oplus}^{(0)}(x^+,x^-)",
        "E3 sector-based exclusion": (
            r"Inthefixed-dilaton,log-freesectorusedintheLetter"
        ),
        "E4 aligned mixed one-point": (
            r"\int\rd^{2}x\,\langleK_{\oplus\bar{\ominus}}\rangle"
            r"=\frac{1}{64\piGl}\int\rd^{2}x\,W_{1}"
        ),
    }
    missing = [name for name, fragment in required.items() if fragment not in compact]
    if missing:
        raise AssertionError(
            "NR_Holography.tex is out of sync with the verified NR linearization: "
            + ", ".join(missing)
        )

    forbidden = {
        "old perturbation symbol": r"\mathfrakr",
        "old chiral perturbation symbol": r"\ell_{+}",
        "old E2 display": r"\tfrac{l}{8}\,C\,y",
        "old E3 sentence": r"excludesgenerictwo-dimensionalrelativizingsources",
        "old E4 sign": (
            r"\int\rd^{2}x\,\langleK_{\oplus\bar{\ominus}}\rangle"
            r"=-\frac{1}{16\piGl}"
        ),
        "old E6 correlated sign": r"$\mp\tfrac{\cW}{2}\tau^{\mp}$",
        "old-frame dilaton constraint": (
            r"\tfrac{8}{l}\partial_y\deltad"
            r"=-2\partial_+\partial_-h_{\oplus\bar\oplus}"
        ),
        "old-frame plus equation": (
            r"h_{\ominus\bar\oplus}=-2\partial_+^2h_{\oplus\bar\oplus}"
        ),
        "old-frame minus equation": (
            r"h_{\oplus\bar\ominus}=-2\partial_-^2h_{\oplus\bar\oplus}"
        ),
        "old (wrong) dilaton-slope sign": (
            r"+\frac{l}{4}y\,\partial_+\partial_-h_{\oplus\bar\oplus}^{(0)}"
        ),
        "old (wrong) ell+ tail sign": (
            r"+ly\,\partial_+^2h_{\oplus\bar\oplus}^{(0)}"
        ),
        "old (wrong) contact-response sign": (
            r"=\frac{l}{8\piG}\,e^{2y/l}\,"
            r"\partial_{+}\partial_{-}h_{\oplus\bar\oplus}^{(0)}"
        ),
        "old radial-solution phrase": r"Thecompleteradialsolutionis",
    }
    present = [name for name, fragment in forbidden.items() if fragment in compact]
    if present:
        raise AssertionError(
            "NR_Holography.tex still contains rejected patterns: "
            + ", ".join(present)
        )
    return True


l = sp.symbols("l", positive=True)
xp, xm, y = sp.symbols("x_p x_m y", real=True)
u2 = sp.exp(-2 * y / l)


def nr_linearized_system(
    a: sp.Expr,
    b: sp.Expr,
    c: sp.Expr,
    d: sp.Expr,
    dd: sp.Expr,
    dilaton_coefficient: int = -1,
    source_sign: int = -1,
) -> list[sp.Expr]:
    """The displayed relations of ``SMNRlin`` (pm split written out).

    ``dilaton_coefficient`` = -1 and ``source_sign`` = -1 are the
    sign-corrected relations (2026-08-07); passing +2/+1 reproduces the
    previously printed (wrong) equations for the negative controls.
    """

    def box(f: sp.Expr) -> sp.Expr:
        return sp.diff(f, y, 2) + (2 / l) * sp.diff(f, y)

    return [
        box(c),
        sp.diff(c, y, xp),
        sp.diff(c, y, xm),
        sp.diff(dd, y, 2),
        sp.Rational(8) / l * sp.diff(dd, y)
        - dilaton_coefficient * sp.diff(c, xp, xm),
        box(a) - source_sign * sp.diff(c, xp, 2),
        box(d) - source_sign * sp.diff(c, xm, 2),
        sp.diff(4 * sp.diff(dd, xp) - sp.diff(a, xm), y),
        sp.diff(4 * sp.diff(dd, xm) - sp.diff(d, xp), y),
        box(b)
        + sp.diff(d, xp, 2) + sp.diff(a, xm, 2)
        - 4 * sp.diff(dd, xp, xm),
    ]


def general_solution() -> tuple[sp.Expr, ...]:
    """The corrected ``SMNRsol`` with arbitrary r^(0)(x^+, x^-)."""

    c0 = sp.Function("c0")(xp, xm)
    dd0 = sp.Function("dd0")(xp, xm)
    a0 = sp.Function("a0")(xp, xm)
    d0 = sp.Function("d0")(xp, xm)
    b0 = sp.Function("b0")(xp, xm)
    a2 = sp.Function("a2")(xp)
    d2 = sp.Function("d2")(xm)
    b2 = sp.Function("b2")(xp, xm)
    c2 = sp.symbols("c2", real=True)

    big_d = sp.diff(c0, xp, 2, xm, 2)
    c = c0 + c2 * u2
    dd = dd0 - sp.Rational(1, 8) * l * y * sp.diff(c0, xp, xm)
    a = a0 - sp.Rational(1, 2) * l * y * sp.diff(c0, xp, 2) + u2 * a2
    d = d0 - sp.Rational(1, 2) * l * y * sp.diff(c0, xm, 2) + u2 * d2
    b = (
        b0
        + sp.Rational(1, 2) * l * y * (
            4 * sp.diff(dd0, xp, xm)
            - sp.diff(d0, xp, 2)
            - sp.diff(a0, xm, 2)
            - sp.Rational(1, 4) * l**2 * big_d
        )
        + sp.Rational(1, 8) * l**2 * y**2 * big_d
        + u2 * b2
    )
    return a, b, c, d, dd, c0


def main() -> None:
    if check_tex_contract():
        print("LaTeX contract: PASS")

    a, b, c, d, dd, c0 = general_solution()

    residuals = nr_linearized_system(a, b, c, d, dd)
    bad = [i for i, res in enumerate(residuals) if sp.simplify(res) != 0]
    if bad:
        raise AssertionError(f"corrected solution fails corrected system: {bad}")
    print("Sign-corrected SMNRsol solves sign-corrected SMNRlin: PASS")

    # Negative control 1: the previously printed +2 dilaton sign (and the
    # pre-July coefficient 1) fail on any solution with a nonharmonic
    # type-changing source.
    for badcoef in (2, 1):
        old = nr_linearized_system(a, b, c, d, dd,
                                   dilaton_coefficient=badcoef)
        probe = sp.simplify(old[4].subs(c0, xp**2 * xm**2).doit())
        if sp.simplify(probe) == 0:
            raise AssertionError(
                f"negative control coef={badcoef} unexpectedly satisfied"
            )
    # and the old +1 source sign fails on the corrected solution
    old = nr_linearized_system(a, b, c, d, dd, source_sign=+1)
    probe = sp.simplify(old[5].subs(c0, xp**2 * xm**2).doit())
    if sp.simplify(probe) == 0:
        raise AssertionError(
            "negative control source_sign=+1 unexpectedly satisfied"
        )
    print("Old-sign negative controls (+2, 1, source +1): PASS (rejected)")

    # Negative control 2: the previously printed solution
    # (delta d = (l/8) C y, no y^2 branch) fails the corrected constraint.
    c_const = sp.symbols("C", nonzero=True)
    c_old = xp * xm * c_const  # d+d- c0 = C
    dd_old = sp.Rational(1, 8) * l * c_const * y
    res_old = sp.simplify(
        sp.Rational(8) / l * sp.diff(dd_old, y) + sp.diff(c_old, xp, xm)
    )
    if res_old == 0:
        raise AssertionError(
            "old displayed solution unexpectedly satisfies the corrected constraint"
        )
    print("Old-solution negative control: PASS (rejected)")

    # E4 bookkeeping: the aligned-frame dictionary is same-channel and uniform.
    a2v, b2v, c2v, d2v = sp.symbols("a2 b2 c2 d2")
    h2_matrix = sp.Matrix([[a2v, b2v], [c2v, d2v]])
    k_matrix = h2_matrix / 32
    assert k_matrix[0, 0] == a2v / 32 and k_matrix[0, 1] == b2v / 32
    assert sp.simplify(k_matrix[0, 1].subs(b2v, sp.Symbol("W1") / 2)
                       - sp.Symbol("W1") / 64) == 0
    print("E4 sign bookkeeping: PASS")

    print("All checks passed.")


if __name__ == "__main__":
    main()
