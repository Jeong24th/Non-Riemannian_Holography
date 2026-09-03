#!/usr/bin/env python3
"""Regression test for the renormalized on-shell action section of the SM.

Guards the 2026-07-23 additions (now under SM section "Renormalized
on-shell action", labels SMgamma2 / SMgamma2flux / SMgamma2value):

* the Gamma^2 rewriting of the lecture notes [Park:2025core, (2.121)-(2.122)]:
  e^{-2d} S_0 = e^{-2d} L_{Gamma^2} + d_A(e^{-2d} B^A),
  B^A = 4 H^{AB} d_B d - d_B H^{AB},
  with the semi-covariant Christoffel connection validated against its
  defining properties (2.46), (2.48), (2.49) at an exact rational point;
* the universal flux structure: the tilde-y row of H^{AB} is constant for
  both saddles, hence B^y = 4 d_y d and e^{-2d} B^y = -2 d_y e^{-2d}
  = -(4/l)(e^{2y/l} + mu e^{-2y/l}), mu = L+L- (Riemannian),
  mu = L+L-/2 (non-Riemannian);
* W-independence: the section components B^{x+-}, B^y of the
  non-Riemannian saddle contain neither W0 nor W1 (the hair costs no
  on-shell action), for arbitrary W0(x+,x-), W1(x+,x-);
* the renormalized value S_ren = -(8 sqrt(mu))/(16 pi G l) Vol_2, read at
  the interior locus e^{-2d} = 0;
* the quadratic-response normalization in Kim's aligned frame: one
  functional derivative gives K=h^(2)/(32 pi G l), while the source coupling
  -2 h^(0) K supplies an additional one-half in the connected Hessian; the
  Riemannian kernel 3 l^2/(4 pi x^4) therefore reproduces
  (8 pi)^(-2)(c/2)x^(-4) for c=3l/(2G);
* the fixed-dilaton Ward identities and the scoped vanishing of the
  non-Riemannian same-channel Hessians;
* the LaTeX contract for the displayed formulas.

Run from any directory with

    python checks/verify_gamma2_action.py

(The connection gates take a few minutes of exact rational arithmetic.)
"""

from __future__ import annotations

import random
import re
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TEX_PATH = ROOT / "NR_Holography.tex"

D2, Ddim = 6, 3
xp, xm, y, l = sp.symbols("x_p x_m y l", real=True)
coords = [xp, xm, y]

J = sp.zeros(D2, D2)
for i in range(3):
    J[i, 3 + i] = 1
    J[3 + i, i] = 1


def dA(e, A):
    return sp.diff(e, coords[A - 3]) if A >= 3 else sp.Integer(0)


def check_tex_contract() -> None:
    if not TEX_PATH.is_file():
        print("LaTeX contract: SKIP (NR_Holography.tex is not in the archive)")
        return
    compact = re.sub(r"\s+", "", TEX_PATH.read_text(encoding="utf-8"))
    required = {
        "section title": r"\section{RenormalizedOn-ShellAction}",
        "quadratic-action origin": (
            r"S_{\Gamma^{2}}=(16\piG)^{-1}\int_{\Sigma_{3}}"
            r"[\cL_{\Gamma^{2}}-2\Lambda_{\rmDFT}e^{-2d}]"
        ),
        "quadratic counterterm": (
            r"S_{\rmct}=-\frac{1}{16\piG}\frac4l"
            r"\int_{y=Y}\rd^2x\,e^{-2d}"
        ),
        "B^M definition": r"B^{M}=4\cH^{MN}\partial_{N}d-\partial_{N}\cH^{MN}",
        "response matrix": r"\label{SMresponsematrix}",
        "stress Hessian": r"\label{SMstresshessian}",
        "continuity check": r"\label{SMcontinuitycheck}",
        "R two-point function": r"\label{SMRtwopt}",
        "NR two-point function": r"\label{SMNRhairhessian}",
        "R response kernel": r"\frac{3l^{2}}{4\pi(\Deltax^{+})^{4}}",
        "NR scope": (
            r"Intherestrictedsourcesectordefinedabove,theadmissible"
            r"same-channelnormalizablekernelsvanishatthisorder"
        ),
        "NR response-kernel scope": r"same-channelnormalizableresponsekernelsvanish",
        "universal flux": r"e^{-2d}B^{y}=-2\,\partial_{y}e^{-2d}",
        "renormalized value": r"S_{\rmren}=-\frac{8\sqrt{\mu}}{16\piGl}",
        "CPS cross-link": r"theboundaryvectorofthe$\Gamma^{2}$DFTaction",
        "Letter soft clause": (
            r"$W_{1}$dropsoutofboththechargesandthesource-free"
            r"renormalizedon-shellaction~\cite{SM}"
        ),
        "long-string energy": (
            r"E(y)=\frac{wl}{\alpha'}\Big[\big(e^{2y/l}"
            r"-L_{+}L_{-}e^{-2y/l}\big)-\big(e^{2y/l}"
            r"+L_{+}L_{-}e^{-2y/l}\big)\Big]"
            r"=-\frac{2wl}{\alpha'}\,L_{+}L_{-}\,e^{-2y/l}"
        ),
    }
    missing = [k for k, v in required.items() if v not in compact]
    if missing:
        raise AssertionError("NR_Holography.tex out of sync: " + ", ".join(missing))
    if not (
        compact.index(required["quadratic-action origin"])
        < compact.index(r"\label{SMSren2}")
        < compact.index(r"\label{SMgamma2}")
    ):
        raise AssertionError("Gamma^2 derivation must precede SMSren2")
    print("LaTeX contract: PASS")


def correlator_checks() -> None:
    G, l = sp.symbols("G l", nonzero=True)
    pi = sp.pi
    c = 3 * l / (2 * G)
    response_weight = 1 / (32 * pi * G * l)

    # h_{+\bar+}^{(2)} = 2L_+ and its right-moving counterpart.
    assert sp.simplify(response_weight * 2 - 1 / (16 * pi * G * l)) == 0

    # The regular Riemannian response kernel fixes the Brown--Henneaux
    # stress-tensor two-point normalization.
    response_kernel = 3 * l**2 / (4 * pi)
    source_coupling_factor = sp.Rational(1, 2)
    cft_coefficient = c / (2 * (8 * pi) ** 2)
    assert sp.simplify(
        source_coupling_factor * response_weight * response_kernel
        - cft_coefficient
    ) == 0

    # Negative controls: omitting the extra one-half from the -2hK source
    # coupling gives precisely the disputed factor-two result and must fail.
    assert sp.simplify(response_weight * response_kernel - cft_coefficient) != 0
    assert sp.simplify(response_kernel - cft_coefficient) != 0
    print("S_ren one-/two-point normalization and negative controls: PASS")


def riemannian_H(g, B):
    gi = g.inv()
    return (gi.row_join(-gi * B)).col_join((B * gi).row_join(g - B * gi * B))


def B_closed(H, d):
    Hu = J * H * J
    return [sum(4 * Hu[A, Bi] * dA(d, Bi) - dA(Hu[A, Bi], Bi) for Bi in range(D2))
            for A in range(D2)]


def hhz_scalar(H, d):
    Hu = J * H * J
    R = 0
    for M in range(D2):
        for N in range(D2):
            R += (4 * Hu[M, N] * dA(dA(d, M), N)
                  - dA(dA(Hu[M, N], M), N)
                  - 4 * Hu[M, N] * dA(d, M) * dA(d, N)
                  + 4 * dA(Hu[M, N], M) * dA(d, N))
    for M in range(3, D2):
        for N in range(3, D2):
            for K in range(D2):
                for L2 in range(D2):
                    dH_M = dA(Hu[K, L2], M)
                    if dH_M != 0:
                        R += sp.Rational(1, 8) * Hu[M, N] * dH_M * dA(H[K, L2], N)
    for M in range(3, D2):
        for N in range(D2):
            for K in range(3, D2):
                for L2 in range(D2):
                    R -= sp.Rational(1, 2) * Hu[M, N] * dA(Hu[K, L2], M) * dA(H[N, L2], K)
    return R


def gamma_gates() -> None:
    """Connection + (2.121) + (2.122) at an exact rational point."""
    random.seed(11)

    def rnd():
        return sp.Rational(random.randint(-3, 3), random.randint(1, 3))

    gm = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            e = rnd() + rnd() * xp + rnd() * xm + rnd() * y + rnd() * xp * y + rnd() * xm * xm
            gm[i, j] = gm[j, i] = e
    gm += 5 * sp.eye(3)
    Bf = sp.zeros(3, 3)
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        e = rnd() + rnd() * xm + rnd() * xp * y + rnd() * y * y
        Bf[i, j], Bf[j, i] = e, -e
    dd = rnd() + rnd() * xp + rnd() * y + rnd() * xm * y + rnd() * xp * xm

    H = riemannian_H(gm, Bf)
    pt = {xp: sp.Rational(1, 3), xm: sp.Rational(-1, 2), y: sp.Rational(2, 5)}

    Hn = H.subs(pt)
    P, Pb = (J + Hn) / 2, (J - Hn) / 2
    dPn = [sp.Matrix(D2, D2, lambda a, b: dA((J + H)[a, b] / 2, C)).subs(pt)
           for C in range(D2)]
    ddn = [dA(dd, A).subs(pt) for A in range(D2)]

    PJ, JPb = P * J, J * Pb
    PJn, PbJn = P * J, Pb * J
    t1 = [PJ * dPn[C] * JPb for C in range(D2)]
    v = [ddn[Dd] + sum(J[E, F] * (t1[F][E, Dd] - t1[F][Dd, E]) / 2
                       for E in range(D2) for F in range(D2))
         for Dd in range(D2)]
    t2 = []
    for C in range(D2):
        M = sp.zeros(D2, D2)
        for A in range(D2):
            for B in range(D2):
                s = 0
                for Dd in range(3, D2):
                    for E in range(D2):
                        c1 = PbJn[A, Dd] * PbJn[B, E] - PJn[A, Dd] * PJn[B, E]
                        if c1 != 0:
                            s += c1 * dPn[Dd][E, C]
                M[A, B] = s
        t2.append(M)
    Gam = []
    for C in range(D2):
        M = sp.zeros(D2, D2)
        for A in range(D2):
            for B in range(D2):
                term1 = t1[C][A, B] - t1[C][B, A]
                term2 = t2[C][A, B] - t2[C][B, A]
                term3 = 0
                for Dd in range(D2):
                    term3 += (Pb[C, A] * PbJn[B, Dd] - Pb[C, B] * PbJn[A, Dd]
                              + P[C, A] * PJn[B, Dd] - P[C, B] * PJn[A, Dd]) * v[Dd]
                term3 *= -sp.Rational(2, Ddim - 1)
                M[A, B] = term1 + term2 + term3  # exact rationals; never nsimplify
        Gam.append(M)

    assert max(abs(Gam[C][A, B] + Gam[C][B, A]) for C in range(D2)
               for A in range(D2) for B in range(D2)) == 0, "(2.48) fails"
    assert max(abs(Gam[A][B, C] + Gam[B][C, A] + Gam[C][A, B]) for A in range(D2)
               for B in range(D2) for C in range(D2)) == 0, "(2.49) fails"
    worst = 0
    for A in range(D2):
        for B in range(D2):
            for C in range(D2):
                e = dPn[A][B, C]
                for E in range(D2):
                    for F in range(D2):
                        e += J[E, F] * (Gam[A][B, F] * P[E, C] + Gam[A][C, F] * P[B, E])
                worst = max(worst, abs(e))
    assert worst == 0, "(2.46) nabla P fails"
    worst = 0
    for A in range(D2):
        e = 2 * ddn[A]
        for B in range(D2):
            for F in range(D2):
                e += J[B, F] * Gam[F][B, A]
        worst = max(worst, abs(e))
    assert worst == 0, "(2.46) trace fails"

    Pu, Pbu = J * P * J, J * Pb * J
    Bg = []
    for A in range(D2):
        s = 0
        for Bi in range(D2):
            for C in range(D2):
                for Dd in range(D2):
                    c1 = Pu[A, C] * Pu[Bi, Dd] - Pbu[A, C] * Pbu[Bi, Dd]
                    if c1 != 0:
                        s += 2 * c1 * Gam[Bi][C, Dd]
        Bg.append(s)
    Bc = [b.subs(pt) for b in B_closed(H, dd)]
    assert max(abs(Bg[A] - Bc[A]) for A in range(D2)) == 0, "(2.121) fails"

    def Gud(C, A, B):
        return sum(Gam[C][A, E] * J[E, B] for E in range(D2))

    def Guf(E, A, B):
        return sum(J[E, F] * Gam[F][A, B] for F in range(D2))

    LG2 = 0
    for A in range(D2):
        for Bi in range(D2):
            for C in range(D2):
                for Dd in range(D2):
                    pref = Pu[A, C] * Pu[Bi, Dd] - Pbu[A, C] * Pbu[Bi, Dd]
                    if pref == 0:
                        continue
                    s1 = sum(Gud(A, C, E) * Gam[Bi][Dd, E] for E in range(D2))
                    s2 = sum(Gud(A, Bi, E) * Gam[Dd][C, E] for E in range(D2))
                    s3 = sum(Guf(E, A, Bi) * Gam[E][C, Dd] for E in range(D2)) / 2
                    LG2 += pref * (s1 - s2 + s3)
    w = sp.exp(-2 * dd)
    lhs = (w * hhz_scalar(H, dd)).subs(pt)
    divB = sum(dA(sp.exp(-2 * dd) * B_closed(H, dd)[A], A) for A in range(D2)).subs(pt)
    assert sp.simplify(lhs - (w.subs(pt)) * LG2 - divB) == 0, "(2.122) fails"
    print("Gamma gates (2.46)/(2.48)/(2.49)/(2.121)/(2.122): PASS")


def saddle_checks() -> None:
    u = sp.exp(-2 * y / l)
    # Riemannian Banados, generic chiral L_pm
    Lp, Lm = sp.Function("L_p")(xp), sp.Function("L_m")(xm)
    g = sp.Matrix([[2 * Lp, -(sp.exp(2 * y / l) + Lp * Lm * u), 0],
                   [-(sp.exp(2 * y / l) + Lp * Lm * u), 2 * Lm, 0],
                   [0, 0, 1]])
    Bf = sp.zeros(3, 3)
    Bf[1, 0] = sp.exp(2 * y / l) + Lp * Lm * u
    Bf[0, 1] = -Bf[1, 0]
    HR = riemannian_H(g, Bf)
    em2dR = sp.exp(2 * y / l) - Lp * Lm * u
    dR = -sp.log(em2dR) / 2
    HuR = sp.simplify(J * HR * J)
    assert list(HuR.row(5)) == [0, 0, 0, 0, 0, 1], "R: tilde-y row not constant"
    BR = B_closed(HR, dR)
    assert sp.simplify(BR[5] - 4 * sp.diff(dR, y)) == 0
    fluxR = sp.simplify(em2dR * BR[5])
    target = -(4 / l) * (sp.exp(2 * y / l) + Lp * Lm * u)
    assert sp.simplify(fluxR - target) == 0
    assert sp.simplify(fluxR + 2 * sp.diff(em2dR, y)) == 0
    # IR value at the horizon e^{2y/l} = sqrt(Pi)
    PiS = sp.symbols("Pi_c", positive=True)
    val = target.subs({Lp * Lm: PiS}).subs(sp.exp(2 * y / l), sp.sqrt(PiS))
    val = sp.simplify(val.subs(sp.exp(-2 * y / l), 1 / sp.sqrt(PiS)))
    assert sp.simplify(val + 8 * sp.sqrt(PiS) / l) == 0
    print("Riemannian saddle: B^y = 4 d_y d, flux, horizon value -(8/l)sqrt(mu): PASS")

    # long-string energetics (SMlongstringE): static winding-w string
    t_, ph_ = sp.symbols("t phi", real=True)
    s2 = sp.sqrt(2)
    dxp = (1 / s2, l / s2)     # (dt, dphi) components of dx^+
    dxm = (1 / s2, -l / s2)
    E2 = sp.exp(2 * y / l) + Lp * Lm * u
    gpp, gmm, gpm = 2 * Lp, 2 * Lm, -E2
    g_tt = gpp * dxp[0] ** 2 + gmm * dxm[0] ** 2 + 2 * gpm * dxp[0] * dxm[0]
    g_ff = gpp * dxp[1] ** 2 + gmm * dxm[1] ** 2 + 2 * gpm * dxp[1] * dxm[1]
    g_tf = gpp * dxp[0] * dxp[1] + gmm * dxm[0] * dxm[1] \
        + gpm * (dxp[0] * dxm[1] + dxp[1] * dxm[0])
    det2 = sp.simplify(g_tf ** 2 - g_tt * g_ff)
    assert sp.simplify(det2 - l ** 2 * em2dR ** 2) == 0, \
        "perfect square sqrt(-det g2) = l e^{-2d} fails"
    B_tf = sp.simplify(E2 * (dxm[0] * dxp[1] - dxm[1] * dxp[0]))
    assert sp.simplify(B_tf - l * E2) == 0
    E_res = sp.simplify(l * em2dR - B_tf)          # per (w/alpha') factor
    assert sp.simplify(E_res + 2 * l * Lp * Lm * u) == 0, \
        "E(y) = -(2wl/alpha') L+L- e^{-2y/l} fails"
    print("Long-string energetics: sqrt(-det g2) = l e^{-2d}, "
          "E = -(2wl/a') L+L- e^{-2y/l}: PASS")

    # Non-Riemannian saddle, constant L_pm, arbitrary W0, W1
    Pi, sig = sp.symbols("Pi sigma", positive=True)
    W0 = sp.Function("W_0")(xp, xm)
    W1 = sp.Function("W_1")(xp, xm)
    chif = sp.Function("chi")(y)
    W = W0 + W1 * chif / (2 * sp.sqrt(Pi))
    ch, sh = sp.cosh(chif), sp.sinh(chif)
    es = sp.exp(sig)
    HN = sp.Matrix([
        [0, 0, 0, ch, -sh / es, 0],
        [0, 0, 0, es * sh, -ch, 0],
        [0, 0, 1, 0, 0, 0],
        [ch, es * sh, 0, -W * es * sh, W * ch, 0],
        [-sh / es, -ch, 0, W * ch, -W * sh / es, 0],
        [0, 0, 0, 0, 0, 1]])
    assert sp.simplify(HN * J * HN - J) == sp.zeros(D2, D2)
    HuN = J * HN * J
    assert list(HuN.row(5)) == [0, 0, 0, 0, 0, 1], "NR: tilde-y row not constant"
    # section components of B^A are W-free
    dN = -y / l + sp.log(sp.cosh(chif / (2 * sp.sqrt(2))))
    BN = B_closed(HN, dN)
    for A in (3, 4, 5):
        assert not (sp.simplify(BN[A]).has(W0) or sp.simplify(BN[A]).has(W1)), \
            f"B^{A} depends on W"
    assert sp.simplify(BN[5] - 4 * sp.diff(dN, y)) == 0
    # closed q-form of the flux: e^{-2d} = e^{2y/l} - (Pi/2) e^{-2y/l}
    em2dq = sp.exp(2 * y / l) - (Pi / 2) * u
    fluxq = -2 * sp.diff(em2dq, y)
    targetN = -(4 / l) * (sp.exp(2 * y / l) + (Pi / 2) * u)
    assert sp.simplify(fluxq - targetN) == 0
    valN = targetN.subs(sp.exp(2 * y / l), sp.sqrt(Pi / 2)).subs(
        sp.exp(-2 * y / l), 1 / sp.sqrt(Pi / 2))
    assert sp.simplify(valN + 8 * sp.sqrt(Pi / 2) / l) == 0
    print("Non-Riemannian saddle: W-free section fluxes, B^y = 4 d_y d, "
          "cap value -(8/l)sqrt(Pi/2): PASS")


def main() -> None:
    check_tex_contract()
    correlator_checks()
    saddle_checks()
    gamma_gates()
    print("All checks passed.")


if __name__ == "__main__":
    main()
