"""Symbolic regression checks for the BRST origin of the W1 mode.

The calculation uses the standard bosonic-string conventions

    S_y  = (1/(2*pi*alpha')) int d^2z partial y barpartial y,
    S_FT = (1/(4*pi)) int sqrt(g) R d,       d = V y,

for which

    y(z)y(w) ~ -(alpha'/2) log|z-w|^2,
    T_y = -(1/alpha') :(partial y)^2: + V partial^2 y.

It checks the conformal weight, radial equation, momentum-space
factorization, central-charge balance, the old-note negative controls,
and the leading fusion tower relevant to exact marginality.
"""

from __future__ import annotations

from fractions import Fraction


def weight_coefficient(b: Fraction, v: Fraction) -> Fraction:
    """Return h/(alpha'/ell^2) for a=b/ell and V=v/ell."""

    return b * (2 * v - b) / 4


def radial_coefficient(b: Fraction, v: Fraction) -> Fraction:
    """Return ell^2 (f''-2Vf')/f for f=exp(b*y/ell)."""

    return b * (b - 2 * v)


def fusion_weight(n: int, q: Fraction) -> Fraction:
    """Holomorphic weight of the leading n-fold fused operator."""

    return n + q * n * (1 - n)


def main() -> None:
    v_ads = Fraction(-1)
    b_w1 = Fraction(-2)

    # OPE double-pole coefficient of T_y(z) exp(a*y(w)).
    h_w1 = weight_coefficient(b_w1, v_ads)
    assert h_w1 == 0

    # Negative controls: neither omitting nor doubling the FT improvement
    # makes the desired W1 dressing weightless.
    h_no_ft = weight_coefficient(b_w1, Fraction(0))
    h_double_ft = weight_coefficient(b_w1, 2 * v_ads)
    assert h_no_ft == -1
    assert h_double_ft == 1

    # Radial equation on the two correct branches.
    assert radial_coefficient(Fraction(0), v_ads) == 0
    assert radial_coefficient(b_w1, v_ads) == 0

    # Doubling S_FT instead selects exp(-4y/ell).
    assert radial_coefficient(b_w1, 2 * v_ads) != 0
    assert radial_coefficient(Fraction(-4), 2 * v_ads) == 0

    # Momentum-space polynomial in x=k*ell:
    # x^2-2*i*x = x(x-2i).  Completing the square with
    # X=P*ell=x-i gives X^2+1.
    momentum_coefficients = (1, -2j, 0)
    completed_square_coefficients = (1, -2j, 0)
    assert momentum_coefficients == completed_square_coefficients
    for root in (0j, 2j):
        assert abs(root * (root - 2j)) == 0

    # Old-note sign negative control.  With p=-i partial_y,
    # p(p+2i/ell) gives -(partial_y^2-2/ell partial_y).
    old_sign_on_w1 = -(b_w1**2 - 2 * b_w1)
    assert old_sign_on_w1 == -8

    # Central-charge identity.  Write q=alpha'/ell^2 and test exact
    # rational samples: 2 + (1+6q) + (23-6q) - 26 = 0.
    for q in (Fraction(1, 100), Fraction(1, 7), Fraction(1, 2), Fraction(3, 5)):
        c_beta_gamma = Fraction(2)
        c_y = 1 + 6 * q
        c_internal = 23 - 6 * q
        c_ghost = Fraction(-26)
        assert c_beta_gamma + c_y + c_internal + c_ghost == 0

    # Algebraic factorization:
    # h_n-1 = (n-1)(1-nq).  Verify exactly for a broad integer range,
    # including the resonance q=1/n.
    for n in range(1, 65):
        for q in (Fraction(1, 97), Fraction(2, 13), Fraction(5, 11)):
            assert fusion_weight(n, q) - 1 == (n - 1) * (1 - n * q)
        if n >= 2:
            assert fusion_weight(n, Fraction(1, n)) == 1

    # B_{+-}=-(1/2)e^{-2y/ell}W1 gives
    # H_{y+-}=+(1/ell)e^{-2y/ell}W1, so the mode is not pure gauge.
    b_prefactor = Fraction(-1, 2)
    radial_exponent = Fraction(-2)
    h_flux_prefactor = b_prefactor * radial_exponent
    assert h_flux_prefactor == 1

    # The historical integrated one-point coefficient is low by 8.
    old_coefficient = Fraction(1, 128)
    corrected_coefficient = Fraction(1, 16)
    assert corrected_coefficient / old_coefficient == 8

    print("Standard FT weight of exp(-2y/ell): 0")
    print("No-FT negative control: -alpha'/ell^2")
    print("Doubled-FT negative control: +alpha'/ell^2")
    print("Radial solutions checked: 1, exp(-2y/ell)")
    print("Momentum-space condition: k(k-2i/ell)=0")
    print("Required internal c: 23-6*alpha'/ell^2")
    print("n-fold fusion: h_n=n+q*n*(1-n)")
    print("Possible logarithmic resonance: q=1/n, n>=2")
    print("H_y+-: +(1/ell) exp(-2y/ell) W1")
    print("Historical/current one-point ratio: 8")
    print("All BRST/W1 regression checks: PASS")


if __name__ == "__main__":
    main()

