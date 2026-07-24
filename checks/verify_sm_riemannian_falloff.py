#!/usr/bin/env python3
"""Regression test for the Riemannian SM linearized falloff.

This test is deliberately tied to both the equations and their LaTeX
presentation.  It rejects the previously printed one-sided logarithmic
falloff, verifies that the corrected falloff satisfies all ten components of
``SMlinEDFE``, and checks that ``NR_Holography.tex`` contains the corrected
signs, logarithms, and response constraint.

Run from any directory with

    python calculations/verify_sm_riemannian_falloff.py
"""

from __future__ import annotations

import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = ROOT / "NR_Holography.tex"


def check_tex_contract() -> None:
    """Keep the verified equations and the displayed falloff synchronized."""

    compact = re.sub(r"\s+", "", TEX_PATH.read_text(encoding="utf-8"))
    required = {
        "correct h_minus radial sign": (
            r"\partial_{y}\big(2h_{+}+h_{-}-4\deltad\big)"
        ),
        "correct constraint equation": (
            r"+4h_{+}+2h_{-}-8\deltad\big]"
        ),
        "opposite h_mixed logarithm": (
            r"h_{\ominus\bar{\oplus}}=4\deltad^{(0)}"
            r"+e^{-2y/l}\big[h^{(2)}_{\ominus\bar{\oplus}}"
            r"+\tfrac{4y}{l}\,\deltad^{(2)}\big]+\ldots"
        ),
        "mixed response relation": (
            r"h^{(2)}_{\oplus\bar{\ominus}}"
            r"+h^{(2)}_{\ominus\bar{\oplus}}=4\deltad^{(2)}"
        ),
    }
    missing = [name for name, fragment in required.items() if fragment not in compact]
    if missing:
        raise AssertionError(
            "NR_Holography.tex is out of sync with the verified falloff: "
            + ", ".join(missing)
        )


l = sp.symbols("l", positive=True)
xp, xm, y = sp.symbols("x_p x_m y", real=True)
u2 = sp.exp(-2 * y / l)


def field(name: str) -> sp.Expr:
    return sp.Function(name)(xp, xm)


hpp0, hmm0, dd0 = field("hpp0"), field("hmm0"), field("dd0")
hpp2, hmm2 = field("hpp2"), field("hmm2")
hpm2, hmp2, dd2 = field("hpm2"), field("hmp2"), field("dd2")


def sm_linearized_system(
    hpp: sp.Expr,
    hpm: sp.Expr,
    hmp: sp.Expr,
    hmm: sp.Expr,
    ddf: sp.Expr,
) -> list[sp.Expr]:
    """The nine displayed EDFEs, with the final +/- equation split in two."""

    h_plus = hpm + hmp
    h_minus = hpm - hmp
    dy = lambda f, k=1: sp.diff(f, *([y] * k))
    dp = lambda f, k=1: sp.diff(f, *([xp] * k))
    dm = lambda f, k=1: sp.diff(f, *([xm] * k))

    equations = [
        l**2 * dy(ddf, 2) + l * dy(hpm) + 2 * hpm,
        dy(hpp, 2)
        + 2 / l * dy(hpp)
        + u2 * dp(h_plus - 4 * ddf, 2),
        dy(hmm, 2)
        + 2 / l * dy(hmm)
        + u2 * dm(h_plus - 4 * ddf, 2),
        dy(h_plus, 2)
        + 2 / l * dy(2 * h_minus + h_plus + 4 * ddf)
        + 8 / l**2 * (h_plus + h_minus)
        + 2
        * u2
        * (dm(hpp, 2) + dp(hmm, 2) - 4 * sp.diff(ddf, xp, xm)),
        dy(h_minus, 2)
        + 2 / l * dy(2 * h_plus + h_minus - 4 * ddf),
        2 * dy(ddf, 2)
        + 8 / l * dy(ddf)
        + u2
        * (
            dm(hpp, 2)
            + dp(hmm, 2)
            + sp.diff(h_plus - 8 * ddf, xp, xm)
        ),
    ]

    bracket = 4 * l * dy(ddf) - l / 2 * dy(h_plus) + h_plus - 4 * ddf
    equations.extend(
        [
            l * sp.diff(dy(hmm), xp) - sp.diff(bracket, xm),
            l * sp.diff(dy(hpp), xm) - sp.diff(bracket, xp),
        ]
    )

    final_bracket = l * dy(h_minus) + 4 * h_plus + 2 * h_minus - 8 * ddf
    equations.extend([sp.diff(final_bracket, xp), sp.diff(final_bracket, xm)])
    return [sp.expand(eq.doit()) for eq in equations]


def reduce_on_shell(expr: sp.Expr) -> sp.Expr:
    """Apply the response, trace, and Ward relations at this order."""

    trace_response = l**2 / 8 * (
        sp.diff(hpp0, xm, 2)
        + sp.diff(hmm0, xp, 2)
        - 4 * sp.diff(dd0, xp, xm)
    )
    expr = expr.subs(hmp2, 4 * dd2 - hpm2).doit()

    def apply_ward_identity(atom: sp.Expr) -> sp.Expr:
        if not isinstance(atom, sp.Derivative):
            return atom

        counts = dict(atom.variable_count)
        if atom.expr == hmm2 and counts.get(xp, 0) >= 1:
            counts[xp] -= 1
            remaining = [var for var, count in counts.items() for _ in range(count)]
            base = 2 * sp.diff(dd2, xm)
            return sp.diff(base, *remaining) if remaining else base

        if atom.expr == hpp2 and counts.get(xm, 0) >= 1:
            counts[xm] -= 1
            remaining = [var for var, count in counts.items() for _ in range(count)]
            base = 2 * sp.diff(dd2, xp)
            return sp.diff(base, *remaining) if remaining else base

        return atom

    # Repeat because replacing a differentiated response may expose another
    # derivative node after SymPy expands the expression.
    for _ in range(4):
        expr = expr.replace(
            lambda atom: isinstance(atom, sp.Derivative), apply_ward_identity
        ).doit()
        expr = sp.expand(expr)

    expr = expr.subs(dd2, trace_response).doit()
    return sp.simplify(sp.expand(expr))


def main() -> None:
    check_tex_contract()

    hpp = hpp0 + u2 * hpp2
    hmm = hmm0 + u2 * hmm2
    hpm = u2 * (hpm2 - 4 * y / l * dd2)
    ddf = dd0 + u2 * dd2

    # Negative control: this is the falloff that was previously printed.
    hmp_missing_log = 4 * dd0 + u2 * hmp2
    missing_log_residuals = sm_linearized_system(
        hpp, hpm, hmp_missing_log, hmm, ddf
    )
    if reduce_on_shell(missing_log_residuals[4]) == 0:
        raise AssertionError("negative control unexpectedly satisfies the h_- equation")

    # Correct falloff: the mixed logarithms have equal magnitude and opposite sign.
    hmp_corrected = 4 * dd0 + u2 * (hmp2 + 4 * y / l * dd2)
    corrected_residuals = sm_linearized_system(hpp, hpm, hmp_corrected, hmm, ddf)
    reduced = [reduce_on_shell(eq) for eq in corrected_residuals]
    failures = [(index + 1, value) for index, value in enumerate(reduced) if value != 0]
    if failures:
        details = "; ".join(f"T{index}={value}" for index, value in failures)
        raise AssertionError("corrected falloff failed: " + details)

    print("LaTeX contract: PASS")
    print("Missing-log negative control: PASS (rejected)")
    print("Corrected falloff: PASS (10/10 EDFEs)")


if __name__ == "__main__":
    main()

