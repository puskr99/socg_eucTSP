"""Compare the DP with Held-Karp on the stored adversarial instances (data/adversarial_instances.json) and on
random perturbations of them (angle +- N(0, 0.02), radii * exp(N(0, 0.05)), every 5th: one point to the opposite ray).

usage: python validate_adversarial.py SEED PERTURBATIONS_PER_INSTANCE [INDICES=all, e.g. 0,3,6]
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from tsp_exact import hk
from radial import mk, coords, profile
from dp import solve

rng = np.random.default_rng(int(sys.argv[1])); per = int(sys.argv[2])
base = json.load(open(os.path.join(HERE, 'data', 'adversarial_instances.json')))
sel = range(len(base)) if len(sys.argv) < 4 else [int(x) for x in sys.argv[3].split(',')]
bad = tot = 0; t0 = time.time(); mxw = 0
for bi in sel:
    b = base[bi]; th = b['theta']; r = np.array(b['r']); lab = np.array(b['line'])
    for k in range(per + 1):
        if k == 0:
            th2, r2 = th, r.copy()
        else:
            th2 = float(np.clip(th + rng.normal(0, 0.02), 0.01, np.pi / 2)); r2 = r * np.exp(rng.normal(0, 0.05, len(r)))
            if k % 5 == 0:
                i = rng.integers(len(r)); r2[i] = -r2[i]
        pts = mk(th2, r2[lab == 0], r2[lab == 1]); P = coords(th2, pts); c, t = hk(P)
        mxw = max(mxw, max(sum(om) for _, om, _ in profile(th2, pts, t)))
        v = solve(th2, pts, ct=4); tot += 1
        ok = abs(v - c) <= 1e-7 * max(1, c); bad += not ok
        print(b['name'], k, 'ok' if ok else 'MISMATCH %r %r' % (c, v), flush=True)
print('TOTAL', tot, 'mismatches', bad, 'max wedge cut edges at one level', mxw, 'time %.0f s' % (time.time() - t0))
