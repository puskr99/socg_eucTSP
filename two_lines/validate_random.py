"""Compare the radial-sweep DP (dp.py) with Held-Karp on random instances on two crossing lines.

usage: python validate_random.py SEED COUNT NMIN NMAX [TOTAL_WEDGE_CAP=4]
Modes cycle through uniform, log-radial, clustered, one-sided and integer (ties across rays) positions;
every 7th instance uses the right angle.  Prints any mismatch and a summary line.
"""
import math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tsp_exact import hk
from radial import mk, coords
from dp import solve

seed, N, nlo, nhi = (int(x) for x in sys.argv[1:5])
ct = int(sys.argv[5]) if len(sys.argv) > 5 else 4
rng = np.random.default_rng(seed)
bad = 0; t0 = time.time(); mxs = 0
for k in range(N):
    n = int(rng.integers(nlo, nhi + 1)); th = rng.uniform(0.02, np.pi / 2)
    if k % 7 == 6:
        th = np.pi / 2
    mode = k % 5
    if mode == 0: r = rng.uniform(-1, 1, n)
    elif mode == 1: r = rng.choice([-1, 1], n) * np.exp(rng.uniform(-4, 0, n))
    elif mode == 2: r = rng.choice([-1, 1], n) * (rng.integers(1, 4, n) + rng.uniform(0, 0.05, n))
    elif mode == 3: r = rng.uniform(0, 1, n)
    else: r = rng.choice([-1, 1], n) * rng.integers(1, 6, n).astype(float)
    lab = rng.integers(0, 2, n)
    if lab.sum() in (0, n):
        lab[0] = 1 - lab[0]
    pts = mk(th, r[lab == 0], r[lab == 1]); P = coords(th, pts); c, t = hk(P)
    st = []; v = solve(th, pts, ct=ct, stats=st); mxs = max(mxs, max(st))
    if abs(v - c) > 1e-7 * max(1, c):
        bad += 1; print('MISMATCH', k, n, th, c, v, flush=True)
print('done', N, 'mismatches', bad, 'max states/level', mxs, 'time %.1f s' % (time.time() - t0))
