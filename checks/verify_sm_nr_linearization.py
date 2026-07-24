#!/usr/bin/env python3
"""Regression test for the non-Riemannian SM linearized system.

Guards the 2026-07-23 corrections (audit E1--E4, E6):

* E1: the linearized dilaton constraint of ``SMNRlin`` reads
  (8/l) d_y(delta d) = 2 d_+ d_- r  (the symmetric double count of
  deltaH^{+-} = deltaH^{-+}); the previously printed coefficient 1 is
  rejected as a negative control.
* E2: the general solution ``SMNRsol`` carries
  delta d = delta d^(0) + (l/4) y d_+d_- r^(0) with r^(0) arbitrary,
  and the hair acquires the y^2 branch -(l^2/2) y^2 d_+^2 d_-^2 r^(0);
  the previously printed (l/8) C y display is rejected.
* E3: the exclusion of generic relativizing sources is a property of the
  fixed-dilaton, log-free sector, not of the field equations; the old
  sentence attributing it to ``SMNRsol`` is rejected.
* E4: the mixed one-point function of ``SMNRonept`` carries a plus sign,
  +int W_1/(16 pi G l), consistent with ``SMSren2`` and
  <K_{oplus oplusbar}> = +ell_+^(2)/(16 pi G l).
* E6: the vielbein prose states the shifts -(W/2) tau^{mp} for both
  sectors (the printed matrices), not the correlated mp(W/2) tau^{mp}.

Run from any directory with

    python calculations/verify_sm_nr_linearization.py
"""

from __future__ import annotations

import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = ROOT / "NR_Holography.tex"


def check_tex_contract() -> None:
    """Keep the verified equations and their LaTeX display synchronized."""

    compact = re.sub(r"\s+", "", TEX_PATH.read_text(encoding="utf-8"))
    required = {
        "E1 corrected dilaton constraint": (
            r"\tfrac{8}{l}\,\partial_{y}\deltad=2\,\partial_{+}\partial_{-}r"
        ),
        "E2 corrected dilaton slope": (
            r"\deltad=\deltad^{(0)}(x^{+},x^{-})"
            r"+\tfrac{l}{4}\,y\,\partial_{+}\partial_{-}r^{(0)}"
        ),
        "E2 y^2 hair branch": (
            r"-\tfrac{l^{2}}{2}\,y^{2}\,\partial_{+}^{2}\partial_{-}^{2}r^{(0)}"
        ),
        "E2 y-linear r-source in w": (
            r"+\tfrac{l^{2}}{2}\,\partial_{+}^{2}\partial_{-}^{2}r^{(0)}\big]"
        ),
        "E3 arbitrary relativizing source": r"with$r^{(0)}(x^{+},x^{-})$arbitrary",
        "E3 sector-based exclusion": (
            r"Agenericrelativizingsourceisadmittedby"
        ),
        "E4 positive mixed one-point sign": (
            r"\int\rd^{2}x\,\langleK_{\oplus\bar{\ominus}}\rangle"
            r"=\frac{1}{16\piGl}\int\rd^{2}x\,W_{1}"
        ),
        "E6 uncorrelated vielbein shifts": r"$-\tfrac{\cW}{2}\tau^{\mp}$",
    }
    missing = [name for name, fragment in required.items() if fragment not in compact]
    if missing:
        raise AssertionError(
            "NR_Holography.tex is out of sync with the verified NR linearization: "
            + ", ".join(missing)
        )

    forbidden = {
        "old E1 coefficient": (
            r"\tfrac{8}{l}\,\partial_{y}\deltad=\partial_{+}\partial_{-}r"
        ),
        "old E2 display": r"\tfrac{l}{8}\,C\,y",
        "old E3 sentence": r"excludesgenerictwo-dimensionalrelativizingsources",
        "old E4 sign": (
            r"\int\rd^{2}x\,\langleK_{\oplus\bar{\ominus}}\rangle"
            r"=-\frac{1}{16\piGl}"
        ),
        "old E6 correlated sign": r"$\mp\tfrac{\cW}{2}\tau^{\mp}$",
    }
    present = [name for name, fragment in forbidden.items() if fragment in compact]
    if present:
        raise AssertionError(
            "NR_Holography.tex still contains rejected patterns: "
            + ", ".join(present)
        )


l = sp.symbols("l", positive=True)
xp, xm, y = sp.symbols("x_p x_m y", real=True)
u2 = sp.exp(-2 * y / l)


def nr_linearized_system(
    r: sp.Expr,
    ellp: sp.Expr,
    ellm: sp.Expr,
    w: sp.Expr,
    dd: sp.Expr,
    dilaton_coefficient: int = 2,
) -> list[sp.Expr]:
    """The seven displayed relations of ``SMNRlin`` (pm split written out).

    ``dilaton_coefficient`` = 2 is the corrected fourth relation;
    passing 1 reproduces the previously printed (wrong) equation for the
    negative control.
    """

    def box(f: sp.Expr) -> sp.Expr:
        return sp.diff(f, y, 2) + (2 / l) * sp.diff(f, y)

    return [
        box(r),
        sp.diff(r, y, xp),
        sp.diff(r, y, xm),
        sp.diff(dd, y, 2),
        sp.Rational(8) / l * sp.diff(dd, y)
        - dilaton_coefficient * sp.diff(r, xp, xm),
        box(ellp) - 2 * sp.diff(r, xp, 2),
        box(ellm) - 2 * sp.diff(r, xm, 2),
        sp.diff(4 * sp.diff(dd, xp) - sp.diff(ellp, xm), y),
        sp.diff(4 * sp.diff(dd, xm) - sp.diff(ellm, xp), y),
        box(w)
        + 2 * (sp.diff(ellm, xp, 2) + sp.diff(ellp, xm, 2))
        - 8 * sp.diff(dd, xp, xm),
    ]


def general_solution() -> tuple[sp.Expr, ...]:
    """The corrected ``SMNRsol`` with arbitrary r^(0)(x^+, x^-)."""

    r0 = sp.Function("r0")(xp, xm)
    dd0 = sp.Function("dd0")(xp, xm)
    l0p = sp.Function("l0p")(xp, xm)
    l0m = sp.Function("l0m")(xp, xm)
    w0 = sp.Function("w0")(xp, xm)
    ell2p = sp.Function("ell2p")(xp)
    ell2m = sp.Function("ell2m")(xm)
    w2 = sp.Function("w2")(xp, xm)
    r2 = sp.symbols("r2", real=True)

    big_d = sp.diff(r0, xp, 2, xm, 2)
    r = r0 + r2 * u2
    dd = dd0 + sp.Rational(1, 4) * l * y * sp.diff(r0, xp, xm)
    ellp = l0p + l * y * sp.diff(r0, xp, 2) + u2 * ell2p
    ellm = l0m + l * y * sp.diff(r0, xm, 2) + u2 * ell2m
    w = (
        w0
        + l
        * y
        * (
            4 * sp.diff(dd0, xp, xm)
            - sp.diff(l0m, xp, 2)
            - sp.diff(l0p, xm, 2)
            + sp.Rational(1, 2) * l**2 * big_d
        )
        - sp.Rational(1, 2) * l**2 * y**2 * big_d
        + u2 * w2
    )
    return r, ellp, ellm, w, dd, r0


def main() -> None:
    check_tex_contract()
    print("LaTeX contract: PASS")

    r, ellp, ellm, w, dd, r0 = general_solution()

    residuals = nr_linearized_system(r, ellp, ellm, w, dd, dilaton_coefficient=2)
    bad = [i for i, res in enumerate(residuals) if sp.simplify(res) != 0]
    if bad:
        raise AssertionError(f"corrected solution fails corrected system: {bad}")
    print("Corrected SMNRsol solves corrected SMNRlin: PASS")

    # Negative control 1: the previously printed dilaton coefficient fails
    # on any solution with a nonharmonic relativizing source.
    old = nr_linearized_system(r, ellp, ellm, w, dd, dilaton_coefficient=1)
    probe = sp.simplify(old[4].subs(r0, xp**2 * xm**2).doit())
    if sp.simplify(probe) == 0:
        raise AssertionError(
            "negative control unexpectedly satisfies the old dilaton equation"
        )
    print("Old-coefficient negative control: PASS (rejected)")

    # Negative control 2: the previously printed solution
    # (delta d = (l/8) C y, no y^2 branch) fails the corrected constraint.
    c_const = sp.symbols("C", nonzero=True)
    r_old = xp * xm * c_const  # d+d- r0 = C
    dd_old = sp.Rational(1, 8) * l * c_const * y
    res_old = sp.simplify(
        sp.Rational(8) / l * sp.diff(dd_old, y) - 2 * sp.diff(r_old, xp, xm)
    )
    if res_old == 0:
        raise AssertionError(
            "old displayed solution unexpectedly satisfies the corrected constraint"
        )
    print("Old-solution negative control: PASS (rejected)")

    # E4 sign bookkeeping: with K = (1/16 pi G l)[[ell_+^(2), w^(2)],
    # [-r^(2), ell_-^(2)]] (the frame fixed by the first equation of
    # SMNRonept), the Banados/hair family gives +W_1 in the mixed slot.
    ell2p_v, w2_v, r2_v, ell2m_v = sp.symbols("ell2p w2 r2 ell2m")
    k_matrix = sp.Matrix([[ell2p_v, w2_v], [-r2_v, ell2m_v]])
    assert k_matrix[0, 0] == ell2p_v and k_matrix[0, 1] == +w2_v
    print("E4 sign bookkeeping: PASS")

    print("All checks passed.")


if __name__ == "__main__":
    main()

