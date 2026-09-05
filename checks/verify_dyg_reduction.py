#!/usr/bin/env python3
"""Hard-fail regression test: doubled-yet-gauged string on the exact
non-Riemannian background -> curved first-order (GO-type) action with the
hair coupling; boundary limit -> GO + (SMvertex).

Provenance: session verification 2026-07-24 (scratchpad dyg_reduction.py),
promoted to checks/ after a cross-review of the pending SM
worldsheet-section restructuring (see SM_DYG_WORLDSHEET_REVISION_PLAN.md
and SM_DYG_RESTRUCTURE_DRAFT_CLAUDE.md).  The restructure itself was
applied to the manuscript (P-004 record); the label-scoped LaTeX
contracts required by the draft's apply checklist were added 2026-08-10
(author instruction), so NR_Holography.tex is now READ (never modified)
to pin the (SMdyg)/(SMdygconstraints)/(SMdygGO)/(SMGO)/(SMvertex)
coefficient strings inside their named displays; every symbolic
reduction check remains independent of the manuscript.

Setup (Lorentzian conformal gauge; commuting placeholders dx*, bx* for the
worldsheet derivatives d = d_sigma, db = d_sigmabar; overall 1/(2 pi
alpha') suppressed; doubled ordering (tilde_+, tilde_-, tilde_y; +, -, y)):

  L = -(1/2)[ Dx^A H_AB Dbx^B + Dbx^A H_AB Dx^B ]
      + s ( a_mu dbx^mu - ab_mu dx^mu ),        s = +1 anchored below,

  Dx^A = (d tx_mu - a_mu ; d x^mu),  Dbx^A = (db tx_mu - ab_mu ; db x^mu),

with the exact NR blocks of (NRHcompact): upper-left H^{mu nu} =
diag(0,0,1), mixed block M = Y (x) tau+ - Yb (x) tau-, and the displayed
representative's lower-right block W (tau+ tau- + tau- tau+) + e_y e_y,
where (Morand-Park type-(1,1) data)

  tau+ = (cosh(chi/2), -e^{-sigma} sinh(chi/2)),
  tau- = (-e^{sigma} sinh(chi/2), cosh(chi/2)),
  Y    = (cosh(chi/2),  e^{sigma} sinh(chi/2)),
  Yb   = (e^{-sigma} sinh(chi/2), cosh(chi/2)).

Guarded facts:
  - Morand-Park duality Y.tau+ = Yb.tau- = 1, Y.tau- = Yb.tau+ = 0, and
    longitudinal completeness Y(x)tau+ + Yb(x)tau- = 1.
  - Upper-block degeneracy is precisely longitudinal: H^{mu nu} tau^pm_nu
    = 0 (two-dimensional kernel), while H^{ty ty} = 1 is invertible ---
    so the transverse auxiliary components are Gaussian and the
    longitudinal ones are Lagrange multipliers (NOT "H^{mu nu} = 0").
  - Mixed and lower-right blocks equal the displayed Morand-Park-form
    representative.  Local Milne shifts change these component data while
    leaving the full generalized metric fixed.
  - a_pm variation coefficients are proportional to Y^mu (tau+ . dbx) and
    -Yb^mu (tau- . dx): each lightcone doublet enters through a single
    combination (the orthogonal one is the gauged direction), enforcing
    the curved chirality constraints tau+ . dbx = 0 = tau- . dx.
  - Integrating out the Gaussian pair (a_y, ab_y) yields the Riemannian
    kinetic term: -2 dy dby in the raw normalization (whence the overall
    L_phys = -L/2), i.e. canonical dy dby after normalization; with
    multipliers beta = -(a.Y), betab = +(ab.Yb) (gauge a_- = ab_+ = 0),
      L_phys = dy dby + beta (tau+ . dbx) + betab (tau- . dx)
               + (W/2)[ (tau+ . dx)(tau- . dbx) + (tau- . dx)(tau+ . dbx) ],
    i.e. the lower-right W block descends SYMMETRICALLY, valid at every radius on the
    exact non-Riemannian saddle (curved first-order action; the free GO CFT
    statements apply only in the boundary limit).
  - The unit-Jacobian shift beta -> beta - (W/2)(tau- . dx) removes the
    doubly-constrained bilinear; chi -> 0 then gives exactly
    (SMGO) + (W/2) dx+ dbx-  ==  (SMGO) + (SMvertex)  [(1/4 pi alpha') W].
  - Consistency with the B-pullback route: the finite-chi clock-frame
    identity tau+ wedge tau- = dx+ wedge dx- first identifies the clock
    antisymmetric form with the actual coordinate pullback.  Then
    sym - pullback = W (tau- . dx)(tau+ . dbx), a product of the two
    constraints (multiplier-absorbable).  The antisymmetric form
    (SMdeltaL) belongs to the pullback route, not to the direct reduction.
  - Orientation: s = -1 exchanges the constraint chiralities; the
    ABSOLUTE worldsheet sign remains withheld (manuscript orientation
    unfixed), only relative statements are guarded here.

Exit codes: 0 = ALL CHECKS PASSED; 1 = at least one check failed;
2 = missing dependency.  Usage:
    python checks/verify_dyg_reduction.py [--selftest-fail] [--strict-pin]
(--selftest-fail deliberately flips the expectation of the chi->0 GO+W
check to demonstrate hard failure; it must exit 1.  --strict-pin turns a
SymPy version different from the pin into a dependency error, exit 2.)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

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


BE_RE = re.compile(r"\\be(?![A-Za-z])")
EE_RE = re.compile(r"\\ee(?![A-Za-z])")


def _labelled_equation(text: str, label: str) -> str:
    """Label-scoped ``\\be ... \\label{label} ... \\ee`` extraction.

    Same contract semantics as verify_response_dictionary.py (C-033): a
    fragment relocated to a different equation, or a renamed/missing
    label, must fail; ``(?![A-Za-z])`` keeps ``\\begin`` from matching
    the project display macro ``\\be``.
    """
    token = rf"\label{{{label}}}"
    if text.count(token) != 1:
        raise AssertionError(f"label {label!r} must occur exactly once")
    label_pos = text.index(token)
    starts = [m.start() for m in BE_RE.finditer(text, 0, label_pos)]
    if not starts:
        raise AssertionError(f"no display start before label {label!r}")
    start = starts[-1]
    if EE_RE.search(text, start, label_pos):
        raise AssertionError(f"label {label!r} is outside its display")
    if BE_RE.search(text, start + 3, label_pos):
        raise AssertionError(f"nested display start before label {label!r}")
    end_match = EE_RE.search(text, label_pos)
    if not end_match:
        raise AssertionError(f"no display end after label {label!r}")
    return text[start:end_match.end()]


def tex_contracts() -> None:
    """Pin the applied-restructure displays to this file's derivations.

    Derivation artifacts (rule 7(vi)): every fragment below is the
    LaTeX image of a quantity derived by the symbolic checks in main()
    -- (SMdyg) normalization 1/(4 pi alpha') anchors the reduction
    input; (SMdygconstraints) are the multiplier constraints; (SMdygGO)
    is L_phys after the unit-Jacobian shift; (SMGO)+(SMvertex) are its
    chi->0 limit with the inherited (1/4 pi alpha') W coefficient.
    """
    tex_path = Path(__file__).resolve().parents[1] / "NR_Holography.tex"
    fixture = r"\be RIGHT \label{target}\ee \be NEEDLE \label{other}\ee"
    check("parser control: wrong-equation fragment invisible",
          "NEEDLE" in _labelled_equation(fixture, "target"), False)
    missing_rejected = False
    try:
        _labelled_equation(fixture, "missing")
    except AssertionError:
        missing_rejected = True
    check("parser control: missing label rejected", missing_rejected, True)
    if not tex_path.is_file():
        print("LaTeX contracts: SKIP (NR_Holography.tex is not in the archive)")
        return
    tex = tex_path.read_text(encoding="utf-8")
    contracts = {
        "SMdyg": [
            r"S_{\rmdyg}=\frac{1}{4\pi\alpha'}",
            r"-\epsilon^{\alpha\beta}D_{\alpha}x^{M}\cA_{\betaM}",
        ],
        "SMdygconstraints": [
            r"\tau^{+}_{\mu}\bar\partialx^{\mu}=0",
            r"\tau^{-}_{\mu}\partialx^{\mu}=0",
        ],
        "SMdygGO": [
            r"L_{\rmNR}=\frac{1}{2\pi\alpha'}\Big(\partialy\,\bar\partialy",
            r"+\half\cW\,\tau^{+}_{\mu}\partialx^{\mu}\,"
            r"\tau^{-}_{\nu}\bar\partialx^{\nu}",
        ],
        "SMGO": [
            r"L_{\rmGO}=\frac{1}{2\pi\alpha'}\big(\beta\,\bar\partialx^{+}"
            r"+\brbeta\,\partialx^{-}+\partialy\,\bar\partialy\big)",
        ],
        "SMvertex": [
            r"V_{\cW}=\frac{1}{4\pi\alpha'}\,\cW(x^{+},x^{-})\,"
            r"\partialx^{+}\bar\partialx^{-}",
        ],
    }
    for label, fragments in contracts.items():
        try:
            block = re.sub(r"\s+", "", _labelled_equation(tex, label))
        except AssertionError as exc:
            check(f"[{label}] labelled display located ({exc})", False, True)
            continue
        for fragment in fragments:
            check(f"[{label}] pins {fragment[:44]}", fragment in block, True)


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
        print("[SELFTEST] deliberately flipping the expectation of the "
              "chi->0 GO+W check")

    chi, sg, W, s = sp.symbols("chi sigma W s")
    ch2, sh2 = sp.cosh(chi / 2), sp.sinh(chi / 2)
    es = sp.exp(sg)

    # Morand-Park type-(1,1) longitudinal data
    X = sp.Matrix([ch2, -sh2 / es])       # tau^+
    Xb = sp.Matrix([-es * sh2, ch2])      # tau^-
    Y = sp.Matrix([ch2, es * sh2])
    Yb = sp.Matrix([sh2 / es, ch2])

    dual = all(sp.simplify(v) == 0 for v in [
        (Y.T * X)[0, 0] - 1, (Yb.T * Xb)[0, 0] - 1,
        (Y.T * Xb)[0, 0], (Yb.T * X)[0, 0]])
    check("Morand-Park duality Y.tau+=Yb.tau-=1, Y.tau-=Yb.tau+=0",
          dual, True)
    check("longitudinal completeness Y(x)tau+ + Yb(x)tau- = 1",
          sp.simplify(Y * X.T + Yb * Xb.T - sp.eye(2)) == sp.zeros(2, 2),
          True)
    check("clock-frame orientation det(tau+,tau-) = +1",
          sp.simplify(sp.Matrix.hstack(X, Xb).det()) == 1, True)

    # upper-left block H^{mu nu} = diag(0,0,1): longitudinal kernel only
    UL = sp.diag(0, 0, 1)
    tauP3 = sp.Matrix([X[0], X[1], 0])
    tauM3 = sp.Matrix([Xb[0], Xb[1], 0])
    ey3 = sp.Matrix([0, 0, 1])
    check("upper block: H tau+ = 0 = H tau- (longitudinal kernel)",
          sp.simplify(UL * tauP3) == sp.zeros(3, 1)
          and sp.simplify(UL * tauM3) == sp.zeros(3, 1), True)
    check("upper block: H e_y = e_y (transverse invertible, NOT H=0)",
          sp.simplify(UL * ey3 - ey3) == sp.zeros(3, 1), True)

    # blocks of (NRHcompact)
    M2 = sp.simplify(Y * X.T - Yb * Xb.T)
    M2paper = sp.Matrix([[sp.cosh(chi), -sp.sinh(chi) / es],
                         [es * sp.sinh(chi), -sp.cosh(chi)]])
    check("mixed block of (NRHcompact) = Y(x)tau+ - Yb(x)tau-",
          sp.simplify(M2 - M2paper) == sp.zeros(2, 2), True)
    K2 = sp.simplify(W * (X * Xb.T + Xb * X.T))
    K2paper = sp.Matrix([[-W * es * sp.sinh(chi), W * sp.cosh(chi)],
                         [W * sp.cosh(chi), -W * sp.sinh(chi) / es]])
    check("displayed lower-right W block = W(tau+ tau- + tau- tau+)",
          sp.simplify(K2 - K2paper) == sp.zeros(2, 2), True)

    # worldsheet fields (commuting placeholders)
    dxp, dxm, dy = sp.symbols("dx_p dx_m dy")       # d x^+, d x^-, d y
    bxp, bxm, by = sp.symbols("bx_p bx_m by")       # db x^+, db x^-, db y
    dtp, dtm, dty = sp.symbols("dt_p dt_m dt_y")    # d tx_mu
    btp, btm, bty = sp.symbols("bt_p bt_m bt_y")    # db tx_mu
    ap_, am_, ay_ = sp.symbols("a_p a_m a_y")       # d-components of a_mu
    bp_, bm_, byy_ = sp.symbols("ab_p ab_m ab_y")   # db-components of a_mu

    dx = sp.Matrix([dxp, dxm, dy])
    bx = sp.Matrix([bxp, bxm, by])
    dt = sp.Matrix([dtp, dtm, dty])
    bt = sp.Matrix([btp, btm, bty])
    av = sp.Matrix([ap_, am_, ay_])
    abv = sp.Matrix([bp_, bm_, byy_])

    M3 = sp.zeros(3, 3)
    M3[0:2, 0:2] = M2
    LR = sp.zeros(3, 3)
    LR[0:2, 0:2] = K2
    LR[2, 2] = 1

    Dt, Dbt = dt - av, bt - abv
    half = sp.Rational(1, 2)
    Ldd = -half * ((Dt.T * UL * Dbt)[0, 0] + (Dbt.T * UL * Dt)[0, 0]
                   + (Dt.T * M3 * bx)[0, 0] + (bx.T * M3.T * Dt)[0, 0]
                   + (dx.T * M3.T * Dbt)[0, 0] + (Dbt.T * M3 * dx)[0, 0]
                   + (dx.T * LR * bx)[0, 0] + (bx.T * LR * dx)[0, 0])
    Ltop = s * ((av.T * bx)[0, 0] - (abv.T * dx)[0, 0])
    L = sp.expand(Ldd + Ltop)

    tauP_dx = sp.simplify((X.T * sp.Matrix([dxp, dxm]))[0, 0])
    tauM_dx = sp.simplify((Xb.T * sp.Matrix([dxp, dxm]))[0, 0])
    tauP_bx = sp.simplify((X.T * sp.Matrix([bxp, bxm]))[0, 0])
    tauM_bx = sp.simplify((Xb.T * sp.Matrix([bxp, bxm]))[0, 0])
    clock_wedge = sp.expand(tauP_dx * tauM_bx - tauM_dx * tauP_bx)
    coordinate_wedge = dxp * bxm - dxm * bxp
    check("finite-chi wedge bridge: tau+ wedge tau- = dx+ wedge dx-",
          sp.simplify(clock_wedge - coordinate_wedge) == 0, True)
    check("negative control: opposite clock-frame orientation",
          sp.simplify(clock_wedge + coordinate_wedge) == 0, False)

    # longitudinal multipliers: single surviving combination per doublet
    ok_a = (sp.simplify(sp.diff(L, ap_).subs(s, 1) - 2 * Y[0] * tauP_bx) == 0
            and sp.simplify(sp.diff(L, am_).subs(s, 1)
                            - 2 * Y[1] * tauP_bx) == 0)
    check("a_pm variation = 2 Y^mu (tau+ . dbx)  [multiplier, not Gaussian]",
          ok_a, True)
    ok_ab = (sp.simplify(sp.diff(L, bp_).subs(s, 1)
                         + 2 * Yb[0] * tauM_dx) == 0
             and sp.simplify(sp.diff(L, bm_).subs(s, 1)
                             + 2 * Yb[1] * tauM_dx) == 0)
    check("ab_pm variation = -2 Yb^mu (tau- . dx)  [multiplier, not Gaussian]",
          ok_ab, True)

    # orientation flip: s = -1 exchanges the constraint chiralities
    ok_flip = (sp.simplify(sp.diff(L, ap_).subs(s, -1)
                           + 2 * Yb[0] * tauM_bx) == 0
               and sp.simplify(sp.diff(L, bp_).subs(s, -1)
                               - 2 * Y[0] * tauP_dx) == 0)
    check("s=-1 exchanges constraint chiralities (orientation caveat)",
          ok_flip, True)

    # transverse pair (a_y, ab_y): Gaussian, cross-determined
    solaby = sp.solve(sp.diff(L, ay_).subs(s, 1), byy_)[0]
    solay = sp.solve(sp.diff(L, byy_).subs(s, 1), ay_)[0]
    Ly = sp.expand(L.subs(s, 1).subs({ap_: 0, am_: 0, bp_: 0, bm_: 0,
                                      ay_: solay, byy_: solaby}))
    for v in (dtp, dtm, dty, btp, btm, bty):
        Ly = Ly.subs(v, 0)
    target_y = -2 * dy * by - W * (tauP_dx * tauM_bx + tauM_dx * tauP_bx)
    check("a_y Gaussian: y+W sector = -2 dy dby - W[sym]  "
          "(pre-normalization; fixes L_phys = -L/2)",
          sp.simplify(Ly - target_y) == 0, True)

    # physical normalization and multiplier gauge
    Lfull = sp.expand(L.subs(s, 1).subs({ay_: solay, byy_: solaby}))
    for v in (dtp, dtm, dty, btp, btm, bty):
        Lfull = Lfull.subs(v, 0)
    Lphys = sp.expand(-Lfull / 2)

    beta, betab = sp.symbols("beta betab")
    Lmult = sp.simplify(Lphys.subs({am_: 0, bp_: 0}).subs(
        {ap_: -beta / Y[0], bm_: betab / Yb[1]}))
    target_sym = (dy * by + beta * tauP_bx + betab * tauM_dx
                  + (W / 2) * (tauP_dx * tauM_bx + tauM_dx * tauP_bx))
    check("L_phys = dy dby + beta(tau+.dbx) + betab(tau-.dx) + (W/2)[sym]",
          sp.simplify(Lmult - target_sym) == 0, True)
    target_sym_wrong = (dy * by + beta * tauP_bx + betab * tauM_dx
                        - (W / 2) * (tauP_dx * tauM_bx + tauM_dx * tauP_bx))
    check("negative control: same with -(W/2)[sym]",
          sp.simplify(Lmult - target_sym_wrong) == 0, False)

    # unit-Jacobian multiplier shift removes the doubly-constrained half
    Lshift = sp.simplify(Lmult.subs(beta, beta - (W / 2) * tauM_dx))
    target_shift = (dy * by + beta * tauP_bx + betab * tauM_dx
                    + (W / 2) * tauP_dx * tauM_bx)
    check("shift beta -> beta - (W/2)(tau-.dx) leaves (W/2)(tau+.dx)(tau-.dbx)",
          sp.simplify(Lshift - target_shift) == 0, True)
    Lshift_wrong = sp.simplify(Lmult.subs(beta, beta + (W / 2) * tauM_dx))
    check("negative control: opposite shift sign",
          sp.simplify(Lshift_wrong - target_shift) == 0, False)

    # boundary limit chi -> 0: GO + (SMvertex)
    L0 = sp.simplify(Lshift.subs(chi, 0))
    GO_plus = dy * by + beta * bxp + betab * dxm + (W / 2) * dxp * bxm
    expected_go = False if selftest_fail else True
    check("chi->0: L = GO + (W/2) dx+ dbx-  [= (SMGO)+(SMvertex)]",
          sp.simplify(L0 - GO_plus) == 0, expected_go)
    GO_minus = dy * by + beta * bxp + betab * dxm - (W / 2) * dxp * bxm
    check("negative control: chi->0 with -(W/2) dx+ dbx-",
          sp.simplify(L0 - GO_minus) == 0, False)

    # consistency with the actual coordinate B-pullback (antisymmetric) route
    sym_c = (W / 2) * (tauP_dx * tauM_bx + tauM_dx * tauP_bx)
    pullback_c = (W / 2) * coordinate_wedge
    check("clock antisym = coordinate B-pullback at finite chi",
          sp.simplify((W / 2) * clock_wedge - pullback_c) == 0, True)
    check("direct(sym) - coordinate B-pullback = W (tau-.dx)(tau+.dbx)",
          sp.simplify((sym_c - pullback_c) - W * tauM_dx * tauP_bx) == 0,
          True)

    tex_contracts()

    if FAILURES:
        print("FAILED CHECKS: " + "; ".join(FAILURES))
        sys.exit(1)
    print("ALL CHECKS PASSED")
    sys.exit(0)


if __name__ == "__main__":
    main(selftest_fail="--selftest-fail" in sys.argv[1:],
         strict_pin="--strict-pin" in sys.argv[1:])
