"""Optional (needs Concorde, see ../concorde_tsp.py): radial cut statistics of Concorde optima on two crossing lines,
n points, angles 3, 15, 40, 70, 90 degrees, four radial distributions.  For every instance: the largest prefix cut
(cut), the largest number of wedge cut edges of one orientation (w_or) and in total (sumw) at a level, and the largest
number of along cut edges of one orientation on one ray (al_or).  Stored outputs of earlier runs:
results/concorde_maxcut.json (n = 100..1000, cut only) and results/concorde_wedge_stats.json.

usage: python concorde_two_lines.py SEED OUT.json N1,N2,...
"""
import json, math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from concorde_tsp import solve


def stats(P, ray, t):
    n = len(P); f = [math.hypot(*p) for p in P]; order = sorted(range(n), key=lambda i: (f[i], i))
    rank = [0] * n
    for k, i in enumerate(order):
        rank[i] = k
    W = [[0] * (n + 2) for _ in range(2)]; A = [[[0] * (n + 2) for _ in range(2)] for _ in range(4)]; C = [0] * (n + 2)
    for i in range(n):
        a, b = t[i], t[(i + 1) % n]
        lo, hi = sorted((rank[a], rank[b])); ori = 0 if rank[a] < rank[b] else 1
        C[lo + 1] += 1; C[hi + 1] -= 1
        if ray[a] % 2 != ray[b] % 2:
            W[ori][lo + 1] += 1; W[ori][hi + 1] -= 1
        else:
            o = b if rank[a] < rank[b] else a
            A[ray[o]][ori][lo + 1] += 1; A[ray[o]][ori][hi + 1] -= 1
    best = dict(w_or=0, sumw=0, al_or=0, cut=0); cw = [0, 0]; ca = [[0, 0] for _ in range(4)]; cc = 0
    for m in range(1, n):
        for o in range(2):
            cw[o] += W[o][m]
            for r in range(4):
                ca[r][o] += A[r][o][m]
        cc += C[m]
        if abs(f[order[m]] - f[order[m - 1]]) < 1e-15:
            continue
        best['w_or'] = max(best['w_or'], max(cw)); best['sumw'] = max(best['sumw'], sum(cw))
        best['al_or'] = max(best['al_or'], max(max(x) for x in ca)); best['cut'] = max(best['cut'], cc)
    return best


rng = random.Random(int(sys.argv[1])); out = sys.argv[2]; res = []
for n in [int(x) for x in sys.argv[3].split(',')]:
    for thd in [3, 15, 40, 70, 90]:
        for dist in ['unif', 'log', 'rings', 'near']:
            th = math.radians(thd); u = (1, 0); v = (math.cos(th), math.sin(th)); P = []; ray = []
            for i in range(n):
                line = rng.random() < 0.5; s = rng.choice([-1, 1])
                if dist == 'unif': r = rng.random()
                elif dist == 'log': r = math.exp(rng.uniform(-7, 0))
                elif dist == 'rings': r = (1 + rng.randrange(8)) / 8 + rng.gauss(0, 1e-3)
                else: r = 1 + rng.gauss(0, 0.02) if rng.random() < 0.7 else rng.random()
                d = v if line else u; P.append((s * r * d[0], s * r * d[1])); ray.append((1 if line else 0) + (0 if s > 0 else 2))
            c, t = solve(P, timeout=1800, tries=4)
            st = stats(P, ray, t); st.update(n=n, theta=thd, dist=dist); res.append(st); print(st, flush=True)
            json.dump(res, open(out, 'w'))
