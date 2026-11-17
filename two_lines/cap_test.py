"""The 12-point instance at theta = pi/2 with four wedge cut edges at one level (found by wedge_search.py with
objective 'sumw').  Shows that the DP with the total wedge-end cap lowered from 4 to 3 is no longer exact.

usage: python cap_test.py
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from tsp_exact import hk
from radial import mk, coords, profile
from dp import solve

b = next(x for x in json.load(open(os.path.join(HERE, 'data', 'adversarial_instances.json'))) if x['name'] == 'hw_sumw_11_12')
th = b['theta']; r = np.array(b['r']); lab = np.array(b['line'])
pts = mk(th, r[lab == 0], r[lab == 1]); P = coords(th, pts)
c, t = hk(P)
print('theta = %.6f (pi/2 = %.6f), n = %d' % (th, np.pi / 2, len(pts)))
print('Held-Karp optimum           %.9f' % c)
print('max wedge cut edges / level ', max(sum(om) for _, om, _ in profile(th, pts, t)))
v4 = solve(th, pts, ct=4); v3 = solve(th, pts, ct=3)
print('DP, total wedge cap 4       %.9f' % v4)
print('DP, total wedge cap 3       %.9f  (+%.2f%%)' % (v3, 100 * (v3 / c - 1)))
