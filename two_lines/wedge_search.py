"""Adversarial hill climbing (exact Held-Karp optima) on the wedge cut edges at one level.
Objectives: sumw = total wedge cut edges at a level (the paper's 4-edge instance at theta = pi/2 came from this),
two = number of wedges with two wedge cut edges, cut = largest prefix cut, nz = wedges met.
Writes the best instance found to wedge_search_OBJ_SEED_N.json.

usage: python wedge_search.py SEED N ITERS RESTARTS OBJ [TIME_LIMIT_S]
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tsp_exact import hk
from radial import mk, coords, profile

seed, n, iters, restarts = (int(x) for x in sys.argv[1:5]); obj = sys.argv[5]
tlim = float(sys.argv[6]) if len(sys.argv) > 6 else 1e18
rng = np.random.default_rng(seed)

def score(th, r, lab):
    if lab.sum() in (0, n):
        return -1
    pts = mk(th, r[lab == 0], r[lab == 1]); c, t = hk(coords(th, pts)); res = profile(th, pts, t)
    if obj == 'sumw': v = max(sum(om) for _, om, _ in res); v2 = sum(sum(om) >= v for _, om, _ in res)
    elif obj == 'two': v = max(sum(o == 2 for o in om) for _, om, _ in res); v2 = max(sum(om) for _, om, _ in res)
    elif obj == 'cut': v = max(c for c, _, _ in res); v2 = max(sum(om) for _, om, _ in res)
    else: v = max(sum(o > 0 for o in om) + 0.5 * sum(o == 2 for o in om) for _, om, _ in res); v2 = 0
    return v + 0.01 * v2

best = (-1,); t0 = time.time()
for rs in range(restarts):
    th = rng.uniform(0.05, np.pi / 2); r = rng.uniform(-1, 1, n) * np.exp(rng.uniform(-3, 0, n)); lab = rng.integers(0, 2, n)
    s = score(th, r, lab)
    for it in range(iters):
        th2, r2, lab2 = th, r.copy(), lab.copy(); k = rng.integers(5)
        if k == 0: th2 = float(np.clip(th + rng.normal(0, 0.1), 0.01, np.pi / 2))
        elif k == 1: i = rng.integers(n); r2[i] += rng.normal(0, 0.1) * max(abs(r2[i]), 0.01) * 3
        elif k == 2: i = rng.integers(n); lab2[i] ^= 1
        elif k == 3: r2 *= np.exp(rng.normal(0, 0.3, n))
        else: i = rng.integers(n); r2[i] = -r2[i]
        s2 = score(th2, r2, lab2)
        if s2 >= s:
            th, r, lab, s = th2, r2, lab2, s2
    if s > best[0]:
        best = (s, th, r, lab); print('restart', rs, 'score', s, 'theta %.4f' % th, flush=True)
        json.dump(dict(score=float(s), theta=float(th), r=r.tolist(), line=lab.tolist()),
                  open('wedge_search_%s_%d_%d.json' % (obj, seed, n), 'w'))
    if time.time() - t0 > tlim:
        break
print('best score', best[0])
