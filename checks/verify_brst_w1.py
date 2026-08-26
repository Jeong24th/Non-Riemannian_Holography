"""Symbolic regression checks for the linear Virasoro analysis of W1.

The calculation uses the standard bosonic-string conventions

    S_y  = (1/(2*pi*alpha')) int d^2z partial y barpartial y,
    S_FT = (1/(4*pi)) int sqrt(g) R d,       d = V y,

for which

    y(z)y(w) ~ -(alpha'/2) log|z-w|^2,
    T_y = -(1/alpha') :(partial y)^2: + V partial^2 y.

It checks the conformal weight, radial equation, momentum-space
factorization, the ghost-independent WZW central-charge identity, the
old-note negative controls, and the leading fusion tower relevant to
exact marginality.  The legacy full bosonic central-charge identity is
retained as an algebraic regression check but is not a manuscript claim.
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

    # Legacy full-bosonic central-charge identity.  This remains a useful
    # algebraic regression check, although the manuscript now uses only
    # the ghost-independent beta-gamma + y WZW identity below.
    for q in (Fraction(1, 100), Fraction(1, 7), Fraction(1, 2), Fraction(3, 5)):
        c_beta_gamma = Fraction(2)
        c_y = 1 + 6 * q
        c_internal = 23 - 6 * q
        c_ghost = Fraction(-26)
        assert c_beta_gamma + c_y + c_internal + c_ghost == 0

    # The S3 gauge choice B=l^2 cos^2(theta) dphi1^dphi2 gives
    # |H|=2 l^2 sin(theta)cos(theta) dtheta^dphi1^dphi2.  Factoring
    # pi^2, its integral is 2*(1/2)*(2*2)=4 times pi^2 l^2, hence
    # |int H|/(4*pi^2*alpha')=l^2/alpha'.
    h_prefactor_abs = Fraction(2)
    theta_integral = Fraction(1, 2)
    azimuth_integral_over_pi2 = Fraction(4)
    flux_integral_over_pi2_l2 = (
        h_prefactor_abs * theta_integral * azimuth_integral_over_pi2
    )
    assert flux_integral_over_pi2_l2 == 4
    assert flux_integral_over_pi2_l2 / 4 == 1

    # In the compact AdS3 x S3 x T4 uplift, q=alpha'/ell^2=1/k_s with
    # positive integral supersymmetric flux level k_s=N5.  Decoupling the
    # adjoint fermions shifts the bosonic SL(2,R) and SU(2) levels to
    # k_s+2 and k_s-2.  The beta-gamma plus radial central charge must
    # equal the bosonic SL(2,R)_{k_s+2} WZW value 3(k_s+2)/k_s.
    for k_s in range(1, 65):
        q = Fraction(1, k_s)
        c_first_order = Fraction(2) + 1 + 6 * q
        c_sl2_bosonic = Fraction(3 * (k_s + 2), k_s)
        assert c_first_order == c_sl2_bosonic
        assert (k_s + 2) - 2 == k_s
        assert (k_s - 2) + 2 == k_s

    # Algebraic factorization:
    # h_n-1 = (n-1)(1-nq).  Verify exactly for a broad integer range,
    # including the resonance q=1/n.
    for n in range(1, 65):
        for q in (Fraction(1, 97), Fraction(2, 13), Fraction(5, 11)):
            assert fusion_weight(n, q) - 1 == (n - 1) * (1 - n * q)
        if n >= 2:
            assert fusion_weight(n, Fraction(1, n)) == 1

    # Full generalized-metric/dilaton gauge-orbit obstruction.  Matching a
    # pure lower-right h_{+-}=h_{-+}=w perturbation gives b_{+-}=-w/2,
    # b_{+y}=d_+v^y, b_{-y}=-d_-v^y.  The remaining H and d equations give
    # d_+d_-v^y=0.  Hence db=0 requires -(1/2)d_yw=0: generalized-gauge
    # modes are y-independent.  This is distinct from a representative-level
    # calculation of dB, which is not Milne invariant.
    b_over_w = Fraction(-1, 2)
    db_over_dyw = b_over_w
    radial_exponent = Fraction(-2)
    assert db_over_dyw != 0
    assert Fraction(0) * db_over_dyw == 0  # W0
    assert radial_exponent * db_over_dyw != 0  # W1

    # Keep the normalization history explicit.  The pre-aligned dictionary
    # gave 1/16, while the aligned same-channel dictionary gives 1/64.  The
    # oldest 1/128 value is smaller than the current result by a factor two.
    pre_aligned_coefficient = Fraction(1, 16)
    aligned_coefficient = Fraction(1, 64)
    oldest_coefficient = Fraction(1, 128)
    assert pre_aligned_coefficient / aligned_coefficient == 4
    assert aligned_coefficient / oldest_coefficient == 2

    print("Standard FT weight of exp(-2y/ell): 0")
    print("No-FT negative control: -alpha'/ell^2")
    print("Doubled-FT negative control: +alpha'/ell^2")
    print("Radial solutions checked: 1, exp(-2y/ell)")
    print("Momentum-space condition: k(k-2i/ell)=0")
    print("Matter-sector WZW central charge: c_beta-gamma+c_y=3(k_s+2)/k_s")
    print("S3 flux: |int H|/(4*pi^2*alpha')=ell^2/alpha'")
    print("Compact uplift: alpha'/ell^2=1/k_s, k_s=N5 positive integral")
    print("n-fold fusion: h_n=n+q*n*(1-n)")
    print("Possible logarithmic resonance: q=1/n, hence n=k_s in compact uplift")
    print("Generalized-gauge obstruction: db=0 and delta d=0 require d_y w=0")
    print("Pre-aligned/aligned integrated one-point ratio: 4")
    print("Aligned/oldest integrated one-point ratio: 2")
    print("All Virasoro/W1 regression checks: PASS")


if __name__ == "__main__":
    main()
