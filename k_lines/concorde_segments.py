"""Prefix cuts of stored Concorde optima (data/segments_concorde_optima.json: n = 40..1000 points on k = 2, 3, 4
random segments, three instances per (k, n); Concorde optima of the instances with coordinates scaled by 1e7).
Net: EXACT rational faithful-net closure (faithful_q.closure) of the segments' supporting lines (rationalised).
phi(x) = distance to the net along x's own line; points sorted by phi.  For every prefix cut: OPT edges crossing it;
per group (piece of inner end, piece of outer end) the size (lemma <= 4) and same-direction vertex-disjoint pairs
whose 2-opt exchange would shorten the tour (none expected).  Also checks faithfulness at every input point.
Output: one JSON line per instance; 'pieces' is N and 'maxcross' the largest prefix cut (paper: maxcross <= N).
Stored output: results/concorde_segments_cuts.json.

usage: python concorde_segments.py NMAX [KMAX=4] [OUT.json]
"""
import json, math, sys, collections, signal
from fractions import Fraction as F
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import faithful_q as fq

def rat(x): return F(x).limit_denominator(10**7)

def net_for(P, curves, tlimit=600):
    L = []
    for c in curves:
        a, b = P[c[0]], P[c[-1]]
        L.append(((rat(a[0]), rat(a[1])), (rat(b[0]) - rat(a[0]), rat(b[1]) - rat(a[1]))))
    def handler(s, f): raise TimeoutError
    signal.signal(signal.SIGALRM, handler); signal.alarm(tlimit)
    try:
        k = len(L); pts = []; lab = []; seen = set()
        r = fq.closure(L, max_pts=3000, max_rounds=80)
    except TimeoutError:
        return None, L
    finally:
        signal.alarm(0)
    return r, L

def closure_points(L):
    # rerun closure keeping the points (faithful_q.closure keeps them internally; re-implement minimal wrapper)
    import itertools
    k = len(L); pts = []; seen = set()
    def add(z):
        if z in seen: return False
        seen.add(z); pts.append(z); return True
    for a, b in itertools.combinations(range(k), 2):
        c = fq.crossing(L[a], L[b])
        if c is not None: add(c)
    for i in range(k):
        if not any(fq.crs(fq.sub(z, L[i][0]), L[i][1]) == 0 for z in pts): add(L[i][0])
    while True:
        lab = [{i for i in range(k) if fq.crs(fq.sub(z, L[i][0]), L[i][1]) == 0} for z in pts]
        newf = []
        for i in range(k):
            p, d = L[i]; dd = fq.dot(d, d)
            par = [fq.dot(fq.sub(z, p), d) / dd for z in pts]
            ton = sorted(par[j] for j in range(len(pts)) if i in lab[j])
            for j in range(len(pts)):
                if i in lab[j]: continue
                a = par[j]; w = fq.sub(pts[j], p); H = fq.dot(w, w)/dd - a*a
                if a in ton: continue
                left = [t for t in ton if t < a]; right = [t for t in ton if t > a]
                lo = max(((t + a)/2 + H/(2*(a - t)) for t in left), default=None)
                hi = min(((t + a)/2 - H/(2*(t - a)) for t in right), default=None)
                if lo is None or hi is None or lo < hi: newf.append(fq.foot(pts[j], L[i]))
        if not newf: return pts, lab
        if sum(add(z) for z in newf) == 0: raise RuntimeError('stuck')

def check(r, tlimit=600):
    P = r['P']; n = len(P)
    stat, L = net_for(P, r['curves'], tlimit)
    if stat is None or not stat['ok']: return dict(ok=False, why='net', stat=stat)
    pts, lab = closure_points(L)
    Z = np.array([[float(z[0]), float(z[1])] for z in pts])
    phi = np.zeros(n); piece = [None]*n; ferr = 0.0
    for i, c in enumerate(r['curves']):
        p = np.array([float(L[i][0][0]), float(L[i][0][1])]); u = np.array([float(L[i][1][0]), float(L[i][1][1])]); u /= np.linalg.norm(u)
        on = [j for j, s in enumerate(lab) if i in s]; tz = (Z[on] - p) @ u
        for v in c:
            x = np.array(P[v], float); t = (x - p) @ u
            j = int(np.argmin(np.abs(tz - t))); phi[v] = abs(tz[j] - t)
            ferr = max(ferr, phi[v] - np.min(np.hypot(*(Z - x).T)))
            piece[v] = (i, on[j], int(np.sign(t - tz[j])))
    tour = r['opt']; order = np.argsort(phi, kind='stable'); rank = np.empty(n, int); rank[order] = np.arange(n)
    E = [(tour[i], tour[(i+1) % n]) for i in range(n)]
    maxcross = 0; maxgroup = 0; bad = 0; D = lambda s, t: math.dist(P[s], P[t])
    for m in range(1, n):
        inside = rank < m; G = collections.defaultdict(list); cnt = 0
        for a, b in E:
            if inside[a] != inside[b]:
                cnt += 1; ai, bo = (a, b) if inside[a] else (b, a)
                G[(piece[ai], piece[bo])].append((ai, bo, bool(inside[a])))
        maxcross = max(maxcross, cnt)
        for g in G.values():
            maxgroup = max(maxgroup, len(g))
            for x in range(len(g)):
                for y in range(x+1, len(g)):
                    (a1, b1, d1), (a2, b2, d2) = g[x], g[y]
                    if d1 == d2 and len({a1, a2, b1, b2}) == 4 and D(a1, b1) + D(a2, b2) - D(a1, a2) - D(b1, b2) > 1e-9: bad += 1
    return dict(ok=True, netsize=len(pts), pieces=sum(2*len([1 for s in lab if i in s]) for i in range(len(L))),
                maxcross=maxcross, maxgroup=maxgroup, improving_pairs=bad, faithful_err=float(ferr))

if __name__ == '__main__':
    nmax = int(sys.argv[1]); kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    here = os.path.dirname(os.path.abspath(__file__))
    R = [r for r in json.load(open(os.path.join(here, 'data', 'segments_concorde_optima.json')))
         if int(r['n']) <= nmax and int(r['k']) <= kmax]
    out = []
    for r in sorted(R, key=lambda r: (r['k'], r['n'], r['trial'])):
        x = check(r); x.update(n=r['n'], k=r['k'], trial=r['trial']); out.append(x)
        print(json.dumps({q: x[q] for q in x if q != 'stat'}), flush=True)
    if len(sys.argv) > 3:
        json.dump(out, open(sys.argv[3], 'w'), default=str)
