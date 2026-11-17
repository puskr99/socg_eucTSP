"""Table 1 and Table 2 of the paper.

1. Enumerates the 24 pairs of wedge cut edge types in different wedges, marks the ten excluded by L4 (same outer ray),
   L5 (chains) and L6 (obtuse pair), checks that the fourteen allowed pairs form a triangle-free graph and that no
   three types satisfy (F1)-(F4) of Proposition 12.
2. Evaluates G(e,e') = |io| + |i'o'| - |ii'| - |oo'| on the fourteen witnesses of Table 2 (all must be < 0).
3. Samples random configurations (theta in (0, pi/2]) for every pair: excluded pairs must give min G > 0.

usage: python pair_witnesses.py [SAMPLES=200000] [SEED=0]
"""
import cmath, itertools, math, sys
import numpy as np

R = 'ABCD'
types = [(X, Y) for X in R for Y in R if (R.index(X) - R.index(Y)) % 2 == 1]   # adjacent rays: inner -> outer
wedge = lambda t: frozenset(t)
obtuse = {frozenset('BC'), frozenset('DA')}

def L4(a, b): return a[1] == b[1]
def chain(a, b):   # X->Y in obtuse wedge, Y->Z in acute wedge
    return a[1] == b[0] and wedge(a) in obtuse and wedge(b) not in obtuse
def L5(a, b): return chain(a, b) or chain(b, a)
def L6(a, b): return {a, b} in ({('C', 'B'), ('D', 'A')}, {('A', 'D'), ('B', 'C')})

pairs = [(a, b) for a, b in itertools.combinations(types, 2) if wedge(a) != wedge(b)]
rule = {}
for a, b in pairs:
    rule[(a, b)] = 'L4' if L4(a, b) else 'L5' if L5(a, b) else 'L6' if L6(a, b) else None
name = lambda t: t[0] + '->' + t[1]
print('pairs in different wedges:', len(pairs))
for r in ('L4', 'L5', 'L6'):
    print(r, 'excludes', sorted('{%s, %s}' % (name(a), name(b)) for (a, b), x in rule.items() if x == r))
allowed = [p for p, x in rule.items() if x is None]
print('allowed:', len(allowed))
adj = {frozenset(p) for p in allowed}
tri = [T for T in itertools.combinations(types, 3) if all(frozenset(q) in adj for q in itertools.combinations(T, 2))]
print('triangles in the allowed graph:', tri)
# Proposition 12: largest set of types, at most one per wedge, pairwise allowed
best = max(len(S) for m in range(1, 5) for S in itertools.combinations(types, m)
           if len({wedge(t) for t in S}) == len(S) and all(frozenset(q) in adj for q in itertools.combinations(S, 2)))
print('largest same-orientation set satisfying (F1)-(F4):', best)

def G(a, b, th, x, y, u, v):
    d = {'A': 1, 'B': cmath.exp(1j * th), 'C': -1, 'D': -cmath.exp(1j * th)}
    i, o, i2, o2 = x * d[a[0]], y * d[a[1]], u * d[b[0]], v * d[b[1]]
    return abs(i - o) + abs(i2 - o2) - abs(i - i2) - abs(o - o2)

W = [("AB", "AD", 90, .9, 10, .1, 10), ("AB", "BC", 60, .9, 10, .1, 1), ("AB", "CD", 90, .9, 10, .9, 10),
     ("AB", "DC", 60, .9, 1, .9, 10), ("BA", "AD", 60, .9, 10, .1, 1), ("CB", "AD", 90, .9, 10, .9, 10),
     ("BA", "BC", 90, .1, 10, .9, 10), ("BA", "CD", 60, .9, 1, .9, 10), ("BA", "DC", 90, .9, 10, .9, 10),
     ("DA", "BC", 90, .9, 10, .9, 10), ("CB", "CD", 90, .1, 10, .9, 10), ("CB", "DC", 60, .1, 1, .9, 10),
     ("CD", "DA", 60, .9, 10, .1, 1), ("DA", "DC", 90, .1, 10, .9, 10)]
print('\nTable 2 witnesses (first type, second type, theta, x=|i|, y=|o|, u=|i\'|, v=|o\'|, G):')
covered = set()
for a, b, td, x, y, u, v in W:
    assert max(x, u) < min(y, v)
    g = G(a, b, math.radians(td), x, y, u, v)
    covered.add(frozenset((tuple(a), tuple(b))))
    print('  %s  %s  %3d  %4g %4g %4g %4g  %+.4f  %s' % (name(a), name(b), td, x, y, u, v, g, 'OK' if g < 0 else 'FAIL'))
print('witnesses cover exactly the allowed pairs:', covered == adj)

M = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
th = rng.uniform(1e-3, np.pi / 2, M); t = np.exp(rng.uniform(-3, 3, M))
x, u = t * rng.uniform(0, 1, M), t * rng.uniform(0, 1, M)
y, v = t / rng.uniform(0, 1, M), t / rng.uniform(0, 1, M)
dirs = lambda r: np.exp(1j * ({'A': 0, 'B': th, 'C': np.pi, 'D': np.pi + th}[r]))
print('\nsampled min G over %d configurations with theta in (0, pi/2]:' % M)
for (a, b), r in sorted(rule.items(), key=lambda z: (z[1] is None, str(z[1]))):
    i, o, i2, o2 = x * dirs(a[0]), y * dirs(a[1]), u * dirs(b[0]), v * dirs(b[1])
    g = np.abs(i - o) + np.abs(i2 - o2) - np.abs(i - i2) - np.abs(o - o2)
    print('  %-7s %s, %s: min G = %+.5f' % (r or 'allowed', name(a), name(b), g.min()))
