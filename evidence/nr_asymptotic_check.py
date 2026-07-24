# -*- coding: utf-8 -*-
# Dated amendment (2026-07-24) to claude_nr_asympt.py — original kept
# unmodified per REVIEW_PROTOCOL.md append-only rule.
#
# Codex cross-review finding (CODEX_RESPONSE.md, CLAUDE-C7 item 2):
# the original order classifier substituted exp-factors by u where
# u := exp(-2y/l) itself (a no-op), so atomic exp(-4y/l) residual
# factors were misclassified as O(u^0) and the summary Boolean printed
# False even though every residual is O(e^{-4y/l}).
#
# This amendment maps exponentials to a FORMAL symbol t after
# powsimp(force=True), classifies by Poly degree in t, HARD-FAILS on a
# genuine O(t^0) or O(t^1) residual, and carries a negative control
# that injects a fake O(1) residual and requires the classifier to
# catch it.
import sys
import sympy as sp

l = sp.symbols('l', positive=True)
xp, xm, y = sp.symbols('x_p x_m y', real=True)
u = sp.exp(-2*y/l)
Lp = sp.Function('L_p')(xp); Lm = sp.Function('L_m')(xm)
W1 = sp.Function('W_1')(xp, xm)
ep = sp.Function('epsilon_p')(xp); em = sp.Function('epsilon_m')(xm)
Hinf = sp.Matrix([[0,0,0,1,0,0],[0,0,0,0,-1,0],[0,0,1,0,0,0],
                  [1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,0,0,0,1]])
h1 = sp.zeros(6,6)
h1[0,4] = -2*Lm; h1[4,0] = -2*Lm
h1[1,3] =  2*Lp; h1[3,1] =  2*Lp
h1[3,4] = W1;    h1[4,3] = W1
H = Hinf + u*h1
J = sp.Matrix(sp.BlockMatrix([[sp.zeros(3), sp.eye(3)],
                              [sp.eye(3), sp.zeros(3)]]))
coords = [xp, xm, y]
xi_p  = ep
xi_m  = em
xi_y  = -l/2*(sp.diff(ep,xp) + sp.diff(em,xm))
txi_p = -l**2/2*u*Lp*sp.diff(em, xm, 2)
txi_m = +l**2/2*u*Lm*sp.diff(ep, xp, 2)
txi_y = -l/2*(sp.diff(ep,xp) - sp.diff(em,xm))
xiU = sp.Matrix([txi_p, txi_m, txi_y, xi_p, xi_m, xi_y])
xiD = J*xiU

def dA(f, A):
    if A < 3: return sp.S.Zero
    return sp.diff(f, coords[A-3])

dxi = sp.zeros(6,6)
for A in range(6):
    for C in range(6):
        dxi[A,C] = dA(xiD[C], A) - dA(xiD[A], C)
transport = sp.zeros(6,6)
for A in range(6):
    for Bb in range(6):
        transport[A,Bb] = sum(xiU[C]*dA(H[A,Bb], C) for C in range(6))
LieH = transport + dxi*J*H + H*J*dxi.T
dLp = ep*sp.diff(Lp,xp) + 2*Lp*sp.diff(ep,xp)
dLm = em*sp.diff(Lm,xm) + 2*Lm*sp.diff(em,xm)
dW1 = (ep*sp.diff(W1,xp) + em*sp.diff(W1,xm)
       + 2*W1*(sp.diff(ep,xp)+sp.diff(em,xm))
       - l**2*(Lm*sp.diff(ep,xp,3) + Lp*sp.diff(em,xm,3)))
s = sp.symbols('s')
Hs = H.subs([(Lp, Lp+s*dLp), (Lm, Lm+s*dLm), (W1, W1+s*dW1)])
deltaH = sp.diff(Hs, s).subs(s, 0)
Delta = LieH - deltaH

t = sp.Symbol('t_order', positive=True)
SUB = {sp.exp(-2*y/l): t, sp.exp(-4*y/l): t**2, sp.exp(-6*y/l): t**3,
       sp.exp(2*y/l): 1/t, sp.exp(4*y/l): 1/t**2}

def low_orders(expr):
    """coefficients of t^0 and t^1 after proper formal-symbol mapping;
    also flag any negative power (would signal a real inconsistency)."""
    e = sp.powsimp(sp.expand(expr), force=True).subs(SUB)
    e = sp.expand(e)
    if e == 0:
        return []
    num, den = sp.fraction(sp.cancel(sp.together(e)))
    shift = sp.degree(sp.Poly(den, t)) if den.has(t) else 0
    poly = sp.Poly(sp.expand(num), t)
    bad = []
    for (k,), c in zip(poly.monoms(), poly.coeffs()):
        if k - shift < 2 and sp.simplify(c) != 0:
            bad.append((k - shift, sp.simplify(c)))
    return bad

bad = []
for i in range(6):
    for j in range(6):
        for entry in low_orders(Delta[i, j]):
            bad.append(((i, j),) + entry)
print('falloff preserved to O(e^{-2y/l}) inclusive (residuals only at '
      'O(e^{-4y/l}) or higher):', not bad)
for b in bad[:6]:
    print(' genuine low-order residual', b)

# negative control: inject a fake O(1) residual; classifier must catch it
fake = Delta[3, 3] + Lp*sp.diff(ep, xp)
neg = low_orders(fake)
print('negative control (injected O(1) residual caught):', bool(neg))

if bad or not neg:
    print('FAILED')
    sys.exit(1)
print('ALL CHECKS PASSED')

