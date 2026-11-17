"""Adversarial hill climbing on point positions maximising max_t cut/(4V-8) (and E/(2V-4)) for Held-Karp optima;
prints any violation of Lemmas 22-24 or Theorem 21.

usage: python halfslot_hill.py SEED ROUNDS N
"""
import numpy as np, sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halfslot_check import *
rng = np.random.default_rng(int(sys.argv[1])); rounds = int(sys.argv[2]); n = int(sys.argv[3])
best_all = 0
for r in range(rounds):
    kind, A = make(rng)
    li = rng.integers(len(A.Ls), size=n); sp = rng.uniform(-2, 2, n)
    def pts(li, sp): return np.array([A.Ls[l][0] + s * A.Ls[l][1] for l, s in zip(li, sp)])
    def score(li, sp):
        P = pts(li, sp)
        if len(np.unique(np.round(P, 7), axis=0)) < n: return -1, None
        c, t = hk(P); best = (0, None); viol = None
        for tt in thresholds(A, P, rng):
            res = analyse(A, P, t, tt)
            if res['cut'] > res['bound'] or res['E'] > 2 * res['V'] - 4 or not res['planar'] or res['mult'] or res['twocell'] or res['samel']:
                viol = res
            s = res['cut'] / res['bound'] + 0.01 * res['E'] / (2 * res['V'] - 4)
            if s > best[0]: best = (s, res)
        if viol: print('VIOLATION', kind, viol, li.tolist(), sp.tolist(), flush=True)
        return best
    cur, res = score(li, sp)
    for step in range(300):
        l2, s2 = li.copy(), sp.copy(); j = rng.integers(n)
        if rng.random() < 0.2: l2[j] = rng.integers(len(A.Ls))
        s2[j] += rng.normal() * rng.choice([0.01, 0.1, 0.5])
        v, r2 = score(l2, s2)
        if v >= cur: cur, res, li, sp = v, r2, l2, s2
    best_all = max(best_all, cur)
    print(r, kind, round(cur, 3), res, flush=True)
print('best', best_all)
