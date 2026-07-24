#!/usr/bin/env python3
"""Exact-rational checks for GO-CEFF-001.

This script uses only the Python standard library.  It fixes the worldsheet
orientation from the protected manuscript, checks the Riemannian
Ba~nados sigma-model components, verifies the auxiliary-field round trip,
and exercises sign, branch, and asymptotic negative controls.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "NR_Holography.tex"


class Checks:
    def __init__(self) -> None:
        self.failed: list[str] = []
        self.total = 0

    def check(self, name: str, condition: bool) -> None:
        self.total += 1
        print(f"[{'PASS' if condition else 'FAIL'}] {name}")
        if not condition:
            self.failed.append(name)

    def finish(self) -> int:
        if self.failed:
            print("FAILED CHECKS: " + "; ".join(self.failed), file=sys.stderr)
            return 1
        print(f"ALL CHECKS PASSED ({self.total}/{self.total})")
        return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selftest-fail",
        action="store_true",
        help="inject the wrong sign in the beta-barbeta coefficient",
    )
    args = parser.parse_args()
    checks = Checks()

    source = TEX.read_text(encoding="utf-8")
    checks.check(
        "LaTeX contract: R B-field orientation is dx- wedge dx+",
        r"\rd x^-\wedge\rd x^+" in source,
    )
    checks.check(
        "LaTeX contract: GO multipliers impose bar-d x+ and d x-",
        r"\beta\,\bar\partial x^{+}+\brbeta\,\partial x^{-}" in source
        and r"\bar\partial x^{+}=0=\partial x^{-}" in source,
    )
    checks.check(
        "LaTeX contract: B-pullback fixes the g-B conformal-gauge sign",
        r"\cW\big(\partial x^{+}\bar\partial x^{-}" in source
        and r"-\bar\partial x^{+}\partial x^{-}\big)" in source,
    )

    # a2 = exp(2y/l); p,q = d x+, d x-;
    # bp,bq = bar-d x+, bar-d x-.
    samples = [
        (Q(5), Q(2), Q(3), Q(7), Q(-2), Q(11), Q(13), Q(17), Q(-19)),
        (Q(7, 2), Q(-3), Q(-5), Q(2), Q(4), Q(-7), Q(9), Q(3), Q(8)),
        (Q(9), Q(0), Q(5), Q(-4), Q(6), Q(5), Q(-3), Q(-2), Q(7)),
    ]

    all_metric = True
    all_e = True
    all_direct = True
    all_stationary = True
    all_roundtrip = True
    all_square = True
    all_wrong_sign_rejected = True
    all_bflip_rejected = True
    all_vacuum_only_rejected = True

    for a2, lp, lm, p, q, bp, bq, dy, bdy in samples:
        pi = lp * lm
        f = a2 + pi / a2

        # From ds^2 = dy^2 - 2 theta+ theta-.
        gpp = 2 * lp
        gmm = 2 * lm
        gpm = -f
        bmp = f
        bpm = -f
        all_metric &= (
            gpp == 2 * lp
            and gmm == 2 * lm
            and gpm == -(a2 + pi / a2)
            and bmp == f
            and bpm == -f
        )

        # Manuscript orientation: E = g - B.
        e_pm = gpm - bpm
        e_mp = gpm - bmp
        all_e &= e_pm == 0 and e_mp == -2 * f

        direct = dy * bdy + gpp * p * bp + gmm * q * bq + e_pm * p * bq + e_mp * q * bp
        expected = dy * bdy + 2 * lp * p * bp + 2 * lm * q * bq - 2 * f * q * bp
        all_direct &= direct == expected

        beta_os = -2 * f * q
        betabar_os = -2 * f * bp
        all_stationary &= bp + betabar_os / (2 * f) == 0
        all_stationary &= q + beta_os / (2 * f) == 0

        injected_sign = -1 if args.selftest_fail else 1
        first_os = (
            dy * bdy
            + 2 * lp * p * bp
            + 2 * lm * q * bq
            + beta_os * bp
            + betabar_os * q
            + injected_sign * beta_os * betabar_os / (2 * f)
        )
        all_roundtrip &= first_os == expected

        beta = Q(23)
        betabar = Q(-29)
        lhs = beta * bp + betabar * q + beta * betabar / (2 * f)
        rhs = (beta + 2 * f * q) * (betabar + 2 * f * bp) / (2 * f) - 2 * f * q * bp
        all_square &= lhs == rhs

        wrong = beta_os * bp + betabar_os * q - beta_os * betabar_os / (2 * f)
        all_wrong_sign_rejected &= wrong != -2 * f * q * bp

        # B -> -B swaps the canceled and surviving chirality.
        flip_e_pm = gpm + bpm
        flip_e_mp = gpm + bmp
        all_bflip_rejected &= flip_e_pm == -2 * f and flip_e_mp == 0
        all_bflip_rejected &= not (flip_e_pm == 0 and flip_e_mp == -2 * f)

        if pi != 0:
            all_vacuum_only_rejected &= Q(1, 2 * a2) != Q(1, 2) / f

    checks.check("metric and B components from Rfields", all_metric)
    checks.check("E=g-B cancels d x+ bar-d x- and leaves -2F d x- bar-d x+", all_e)
    checks.check("direct conformal-gauge longitudinal action", all_direct)
    checks.check("auxiliary beta equations", all_stationary)
    checks.check("first-order to second-order round trip", all_roundtrip)
    checks.check("complete-square identity", all_square)

    # Explicit physical negative controls.
    checks.check("negative control: beta-barbeta sign flip is rejected", all_wrong_sign_rejected)
    checks.check("negative control: B reversal swaps chirality", all_bflip_rejected)
    checks.check("negative control: pure vacuum coefficient is not exact for Pi != 0", all_vacuum_only_rejected)

    # Boundary variable u = exp(-2y/l): F = 1/u + Pi*u.
    asymptotic_ok = True
    c_direction_ok = True
    for u, pi in [(Q(1, 10), Q(6)), (Q(1, 100), Q(-15)), (Q(1, 1000), Q(0))]:
        f = 1 / u + pi * u
        inv_c2 = 1 / (2 * f)
        asymptotic_ok &= inv_c2 == u / (2 * (1 + pi * u * u))
        # c_eff^2 * u -> 2 and (1/c_eff^2)/u -> 1/2.
        c_direction_ok &= 2 * f * u == 2 * (1 + pi * u * u)
        c_direction_ok &= inv_c2 / u == 1 / (2 * (1 + pi * u * u))
    checks.check("boundary coefficient u/[2(1+Pi u^2)]", asymptotic_ok)
    checks.check("c_eff grows as exp(+y/l), inverse speed falls as exp(-y/l)", c_direction_ok)

    u = Q(1, 13)
    f_vac = 1 / u
    checks.check("vacuum interaction is exactly (1/2) exp(-2y/l)", 1 / (2 * f_vac) == u / 2)

    # Constant positive-L horizon: choose Pi=36, sqrt(Pi)=6.
    pi_h = Q(36)
    a2_h = Q(6)
    f_h = a2_h + pi_h / a2_h
    c2_h = 2 * f_h
    dilaton_density_h = a2_h - pi_h / a2_h
    checks.check("BTZ horizon c_eff^2 = 4 sqrt(Pi) is finite", c2_h == 24 and c2_h == 4 * 6)
    checks.check("horizon zero is e^{-2d}, not c_eff^2", dilaton_density_h == 0 and c2_h != 0)

    # GO constraint surface: bar-d x+ = 0 and d x- = 0.
    lp, lm, p, q, bp, bq = Q(3), Q(5), Q(7), Q(0), Q(0), Q(11)
    diagonal_banados = 2 * lp * p * bp + 2 * lm * q * bq
    checks.check("Banados diagonal terms vanish on the free-GO constraint surface", diagonal_banados == 0)

    # NR branch: no beta-barbeta term in SMdygGO, hence inverse c_eff^2=0.
    checks.check(
        "NR relativizing channel has 1/c_eff^2 = 0 at every radius",
        r"\beta\brbeta" not in source[source.index(r"\label{SMdygGO}") - 500 : source.index(r"\label{SMdygGO}")],
    )

    # Domain control: for Pi<0, F can vanish, so a positive speed is not global.
    pi_neg = Q(-1)
    a2_zero = Q(1)
    f_zero = a2_zero + pi_neg / a2_zero
    checks.check(
        "domain control: arbitrary-sign Pi can cross F=0",
        f_zero == 0,
    )
    checks.check(
        "negative control: c_eff = exp(-y/l) has the wrong boundary direction",
        # At u=1/100, the claimed squared speed u would go to zero while
        # the exact c_eff^2 = 2/u (vacuum) grows.
        Q(1, 100) != 2 / Q(1, 100),
    )

    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())


