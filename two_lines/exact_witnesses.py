"""Exact check of the fourteen witnesses of Table 3 of the paper (G < 0 for every allowed pair).

Setting (Section 3): rays A = R_{>0}(1,0), B = R_{>0}(cos t, sin t), C = -A, D = -B; a point at distance r from O on
ray X is r*X. For two cut edges e = io (type X->Y, |i| = x, |o| = y) and e' = i'o' (type X'->Y', |i'| = u, |o'| = v)
  G(e,e') = |io| + |i'o'| - |ii'| - |oo'|.
For t = 60 or 90 degrees, cos t is rational (1/2 or 0) and every squared distance r^2 + s^2 - 2rs<X,Y> is rational,
so G is a signed sum of four square roots of rationals. With x, y, u, v in (1/10)Z, 100 * squared distance is an
integer m, and G = (sqrt m1 + sqrt m2 - sqrt m3 - sqrt m4) / 10.

For each witness the script
  1. checks max(x, u) < min(y, v) in exact rational arithmetic;
  2. builds G as an exact sympy expression (radicals) and prints its closed form;
  3. proves G < 0 in two independent ways:
     (a) rigorous rational interval arithmetic: floor/ceil integer square roots (math.isqrt) at scale 10^30
         enclose each sqrt m_j, giving rationals lo <= 10*G <= hi; certified iff hi < 0;
     (b) sympy: the minimal polynomial of G over Q and exact rational isolating intervals of its real roots
         (Poly.intervals, exact arithmetic); the interval from (a) must meet exactly one isolating interval, and
         that interval must lie strictly left of 0;
  4. prints a certified upper bound G <= -c with c rounded down to 3 decimals.

usage: python exact_witnesses.py          (exit status 0 iff all 14 witnesses are certified)
"""
import math, sys
from fractions import Fraction as F
import sympy as sp

DOT = {}   # <X,Y> for unit vectors of the rays, as a function of c = cos t
for X, sx in (('A', 1), ('C', -1)):
    for Y, sy in (('B', 1), ('D', -1)):
        DOT[(X, Y)] = DOT[(Y, X)] = ('c', sx * sy)
DOT[('A', 'C')] = DOT[('C', 'A')] = DOT[('B', 'D')] = DOT[('D', 'B')] = ('1', -1)
for X in 'ABCD':
    DOT[(X, X)] = ('1', 1)

def dot(X, Y, c):
    kind, sgn = DOT[(X, Y)]
    return sgn * (c if kind == 'c' else 1)

COS = {60: F(1, 2), 90: F(0)}
# (first type, second type, theta in degrees, x, y, u, v) as in Table 3
W = [("AB", "AD", 90, '9/10', 10, '1/10', 10), ("AB", "BC", 60, '9/10', 10, '1/10', 1),
     ("AB", "CD", 90, '9/10', 10, '9/10', 10), ("AB", "DC", 60, '9/10', 1, '9/10', 10),
     ("BA", "AD", 60, '9/10', 10, '1/10', 1), ("CB", "AD", 90, '9/10', 10, '9/10', 10),
     ("BA", "BC", 90, '1/10', 10, '9/10', 10), ("BA", "CD", 60, '9/10', 1, '9/10', 10),
     ("BA", "DC", 90, '9/10', 10, '9/10', 10), ("DA", "BC", 90, '9/10', 10, '9/10', 10),
     ("CB", "CD", 90, '1/10', 10, '9/10', 10), ("CB", "DC", 60, '1/10', 1, '9/10', 10),
     ("CD", "DA", 60, '9/10', 10, '1/10', 1), ("DA", "DC", 90, '1/10', 10, '9/10', 10)]

def sqdist(r, X, s, Y, c):
    return r * r + s * s - 2 * r * s * dot(X, Y, c)

SCALE = 10 ** 30

def sqrt_bounds(m):
    """Rationals lo <= sqrt(m) <= hi for a nonnegative integer m."""
    a = math.isqrt(m * SCALE * SCALE)
    return F(a, SCALE), F(a if a * a == m * SCALE * SCALE else a + 1, SCALE)

ok_all = True
for a, b, td, x, y, u, v in W:
    x, y, u, v = F(x), F(y), F(u), F(v)
    c = COS[td]
    cond = max(x, u) < min(y, v)
    terms = [(+1, sqdist(x, a[0], y, a[1], c)), (+1, sqdist(u, b[0], v, b[1], c)),   # |io|, |i'o'|
             (-1, sqdist(x, a[0], u, b[0], c)), (-1, sqdist(y, a[1], v, b[1], c))]   # |ii'|, |oo'|
    ms = []
    for sgn, q in terms:
        m = 100 * q
        assert m.denominator == 1 and m >= 0
        ms.append((sgn, int(m)))
    # (a) rational interval enclosure of 10*G
    lo = hi = F(0)
    for sgn, m in ms:
        l, h = sqrt_bounds(m)
        if sgn > 0:
            lo, hi = lo + l, hi + h
        else:
            lo, hi = lo - h, hi - l
    cert_a = hi < 0
    # (b) sympy: exact expression, minimal polynomial, isolating intervals
    G = sum(sgn * sp.sqrt(sp.Integer(m)) for sgn, m in ms) / 10
    z = sp.Symbol('z')
    mp = sp.Poly(sp.minimal_polynomial(G, z), z)
    glo, ghi = lo / 10, hi / 10
    hits = [(sp.Rational(p), sp.Rational(q)) for (p, q), _ in mp.intervals()
            if not (q < sp.Rational(glo.numerator, glo.denominator) or p > sp.Rational(ghi.numerator, ghi.denominator))]
    # refine until exactly one isolating interval meets the enclosure and it lies left of 0
    # (an irreducible minimal polynomial of degree >= 2 does not vanish at 0, so this terminates)
    eps = sp.Rational(1, 10 ** 6)
    while (len(hits) != 1 or hits[0][1] >= 0) and eps > sp.Rational(1, 10 ** 40):
        hits = [(sp.Rational(p), sp.Rational(q)) for (p, q), _ in mp.intervals(eps=eps)
                if not (q < sp.Rational(glo.numerator, glo.denominator) or p > sp.Rational(ghi.numerator, ghi.denominator))]
        eps /= 10 ** 6
    cert_b = len(hits) == 1 and hits[0][1] < 0
    bound = math.floor(-hi / 10 * 1000) / 1000          # G <= hi/10 <= -bound
    ok = cond and cert_a and cert_b and bound > 0
    ok_all &= ok
    closed = '(%s)/10' % ' '.join(('+' if s > 0 else '-') + 'sqrt(%d)' % m for s, m in ms)
    print('%s->%s  %s->%s  %2d deg  x=%s y=%s u=%s v=%s  max(x,u)<min(y,v): %s' %
          (a[0], a[1], b[0], b[1], td, x, y, u, v, cond))
    print('    G = %s = %s' % (closed, sp.simplify(G)))
    print('    10G in [%.12f, %.12f]  (a) interval: %s   (b) minpoly degree %d, isolating interval [%s, %s]: %s'
          % (float(lo), float(hi), 'G<0' if cert_a else 'FAIL', mp.degree(),
             float(hits[0][0]) if hits else float('nan'), float(hits[0][1]) if hits else float('nan'),
             'G<0' if cert_b else 'FAIL'))
    print('    certified: G <= -%.3f   %s' % (bound, 'OK' if ok else 'FAIL'))

print('\nall 14 witnesses certified exactly:', ok_all and len(W) == 14)
sys.exit(0 if ok_all else 1)
