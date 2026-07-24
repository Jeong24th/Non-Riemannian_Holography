#!/usr/bin/env python3
"""Hard-fail scope checks for the FRAMING-001 Codex prose proposal."""

from __future__ import annotations

import argparse
import sys


TITLE = "Long Strings and Non-Riemannian Hair"

ABSTRACT = r"""
Long strings, macroscopic winding strings near the NS--NS $AdS_3$ boundary,
elude supergravity because their Gomis--Ooguri worldsheet has no Riemannian
target metric. Double Field Theory supplies a bulk description. At vanishing
marginal source, the same type-$(1,1)$ boundary generalized metric admits two
exact saddles: the Riemannian Ba\~nados family and an everywhere
non-Riemannian branch. The doubled-yet-gauged string action couples directly
to a field $\cW$ on the latter: its boundary mode $W_0$ sources an exactly
marginal null vertex, whereas its normalizable mode $W_1$ is soft
winding-charge hair probed by long strings. The saddles share the
Virasoro-sector one-point functions but realize inequivalent finite-radius
branches and, respectively, Brown--Henneaux Virasoro and centerless Witt
charge algebras.
"""

INTRO = r"""
At $W_0=0$, in a fixed boundary representative, the saddles are called two
saddle-resolved radial/RG trajectories only in this branchwise classical
sense; whether their finite-cutoff functionals define distinct Wilsonian
flows or competing Euclidean phases remains open. At the boundary, $W_0$
sources the exactly marginal vertex, while $W_1$ is the normalizable
response. The non-Riemannian charges reduce, after quotienting zero-charge
$B$-gauge directions, to centerless Witt algebras. Arbitrary smooth
Riemannian data preserve sixteen local real polarizations, with global
survival fixed by Hill monodromy and the common spin structure, whereas
generic one-sided non-Riemannian families preserve complementary
four-polarization $\Spin(1,9)$ or $\Spin(9,1)$ sectors.
"""

SM_SECTIONS = (
    "Linearized holographic response",
    "Covariant charges and asymptotic algebras",
    "Renormalized on-shell action",
    "Radial branches and infrared structure",
    "Worldsheet realizations",
    "Ten-dimensional uplift, exact isometries, and supersymmetry",
)


def evaluate(abstract: str, intro: str) -> list[tuple[str, bool]]:
    abstract = " ".join(abstract.split())
    intro = " ".join(intro.split())
    return [
        ("title has no doubled space", "  " not in TITLE),
        ("common source is qualified", "vanishing marginal source" in abstract),
        ("W0 is the boundary source", "boundary mode $W_0$" in abstract),
        ("W1 is normalizable", "normalizable mode $W_1$" in abstract),
        ("charge contrast is present", "centerless Witt" in abstract),
        ("branchwise scope is present", "branchwise classical sense" in intro),
        ("Wilsonian/Euclidean questions remain open", "remains open" in intro),
        ("zero-charge quotient is present", "quotienting zero-charge" in intro),
        ("R SUSY local/global split is present", "sixteen local" in intro and "global" in intro),
        ("NR SUSY is one-sided", "one-sided non-Riemannian" in intro),
        ("official type-II labels are present", r"\Spin(1,9)" in intro and r"\Spin(9,1)" in intro),
        ("six SM topics", len(SM_SECTIONS) == 6),
        ("charge heading names algebra", "asymptotic algebras" in SM_SECTIONS[1]),
        ("worldsheet heading says realization", SM_SECTIONS[4] == "Worldsheet realizations"),
        ("no unqualified two-flow claim", "two distinct Wilsonian flows" not in abstract + intro),
        ("no W1 marginality claim", "$W_1$ is exactly marginal" not in abstract + intro),
        ("no broad same-correlators claim", "same correlators" not in abstract + intro),
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selftest-fail", action="store_true")
    args = parser.parse_args()

    abstract = ABSTRACT
    intro = INTRO
    if args.selftest_fail:
        intro = intro.replace("one-sided non-Riemannian", "non-Riemannian")

    checks = evaluate(abstract, intro)
    failed = []
    for name, passed in checks:
        print(f"[{'PASS' if passed else 'FAIL'}] {name}")
        if not passed:
            failed.append(name)

    if failed:
        print("FAILED CHECKS: " + "; ".join(failed), file=sys.stderr)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

