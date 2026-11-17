"""Adversarial hill climbing over the angle and the points (exact Held-Karp optima) maximising the largest prefix
cut of the radial order.  Writes the best instance found to cut_search_SEED_N.json.

usage: python cut_search.py SEED N ITERS RESTARTS     (paper: N = 12..14; the maximum found was 6)
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tsp_exact import hk
from radial import mk, coords, maxcut

seed, n, iters, restarts = (int(x) for x in sys.argv[1:5])
rng = np.random.default_rng(seed)

def score(th, r, lab):
    if lab.sum() in (0, n):
        return -1
    P = coords(th, mk(th, r[lab == 0], r[lab == 1])); c, t = hk(P)
    return maxcut(P, t)

best = (-1,)
for rs in range(restarts):
    th = rng.uniform(0.05, np.pi / 2); r = rng.uniform(-1, 1, n); lab = rng.integers(0, 2, n)
    s = score(th, r, lab)
    for it in range(iters):
        th2, r2, lab2 = th, r.copy(), lab.copy(); k = rng.integers(4)
        if k == 0: th2 = float(np.clip(th + rng.normal(0, 0.1), 0.01, np.pi / 2))
        elif k == 1: i = rng.integers(n); r2[i] += rng.normal(0, 0.1)
        elif k == 2: i = rng.integers(n); lab2[i] ^= 1
        else: r2 *= np.exp(rng.normal(0, 0.3, n))
        s2 = score(th2, r2, lab2)
        if s2 >= s:
            th, r, lab, s = th2, r2, lab2, s2
    if s > best[0]:
        best = (s, th, r, lab); print('restart', rs, 'max cut', s, 'theta %.4f' % th, flush=True)
        json.dump(dict(maxcut=int(s), theta=float(th), r=r.tolist(), line=lab.tolist()), open('cut_search_%d_%d.json' % (seed, n), 'w'))
print('best max cut', best[0])
