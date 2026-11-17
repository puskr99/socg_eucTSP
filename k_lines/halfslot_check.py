"""Independent check of the linear cut lemma (Theorem 21) through the half-slot graph H of Section 5:
multiplicity (Lemma 22), simplicity, planarity (Lemma 24), |E(H)| <= 2V-4 and cut <= 4(N + Lambda_t) - 8, with exact
level points, on Held-Karp optima of random instances on pencils (2-4 concurrent lines), triangles with the net of
Lemma 27 (vertices + altitude feet), nonobtuse triangles with the vertex net, two parallel lines with a transversal,
and 2-4 parallel lines.  Also reports the largest cut / N.

usage: python halfslot_check.py SEED INSTANCES NMAX      (paper: 3,200 instances, n = 7..13)
"""
import itertools, os, sys
import numpy as np
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tsp_exact import hk
EPS = 1e-9
def mkline(P, Q):
    d = (Q - P) / np.linalg.norm(Q - P); return (np.array(P, float), d)
def crs(a, b): return a[0] * b[1] - a[1] * b[0]
def foot(z, L): p, d = L; return p + np.dot(z - p, d) * d
def onl(z, L): p, d = L; return abs(crs(z - p, d)) < 1e-9
class Arr:
    def __init__(s, Ls, Z):
        s.Ls = Ls; s.Z = [np.array(z, float) for z in Z]; s.Q = []; s.B = []
        for L in Ls:
            p, d = L
            q = sorted(set(round(float(np.dot(z - p, d)), 12) for z in s.Z if onl(z, L)))
            assert q, 'net misses a line'
            s.Q.append(np.array(q))
            s.B.append(np.array(sorted(set(q + [(x + y) / 2 for x, y in zip(q, q[1:])]))))
        s.N = sum(len(b) + 1 for b in s.B)
    def f(s, x): return min(np.linalg.norm(x - z) for z in s.Z)
    def fl(s, li, sp):  # potential on line li at param sp via faithful formula
        return float(np.min(np.abs(s.Q[li] - sp)))
    def piece(s, li, sp): return int(np.searchsorted(s.B[li], sp))
    def lines_of(s, x): return [i for i in range(len(s.Ls)) if onl(x, s.Ls[i])]
    def param(s, li, x): p, d = s.Ls[li]; return float(np.dot(x - p, d))
    def levels(s, t):
        out = []  # (li, piece, param)
        for li in range(len(s.Ls)):
            B = s.B[li]; Q = s.Q[li]
            edges = [-np.inf] + list(B) + [np.inf]
            for pc in range(len(B) + 1):
                lo, hi = edges[pc], edges[pc + 1]
                # f on piece is |s - q| for q = nearest Q; candidates q +- t
                for q in Q:
                    for cand in (q - t, q + t):
                        if lo + EPS < cand < hi - EPS and abs(s.fl(li, cand) - t) < 1e-9:
                            out.append((li, pc, cand))
        out = sorted(set((a, b, round(c, 10)) for a, b, c in out))
        return out
def analyse(A, P, tour, tt):
    """returns dict for threshold tt (assumed generic)"""
    n = len(P); fv = np.array([A.f(x) for x in P])
    lev = A.levels(tt); Lam = len(lev); V = A.N + Lam
    H = {}; cut = 0; bad = 0
    for i in range(n):
        a, b = tour[i], tour[(i + 1) % n]
        if (fv[a] < tt) == (fv[b] < tt): continue
        cut += 1; o = 1 if fv[a] < tt else -1
        if o < 0: a, b = b, a
        xa, xb = P[a], P[b]
        common = [li for li in range(len(A.Ls)) if onl(xa, A.Ls[li]) and onl(xb, A.Ls[li])]
        if common:
            li = common[0]; sa, sb = A.param(li, xa), A.param(li, xb)
            cands = [c for (l2, pc, c) in lev if l2 == li and min(sa, sb) < c < max(sa, sb)]
            c = min(cands, key=lambda c: abs(c - sa)); pc = A.piece(li, c)
            key = ((li, pc, 'L'), (li, pc, 'H')); cell = 'along'
        else:
            e = xb - xa; ws = [(0.0, xa)]
            for li, (p, d) in enumerate(A.Ls):
                den = crs(e, d)
                if abs(den) < 1e-14: continue
                s_ = crs(p - xa, d) / den
                if 1e-9 < s_ < 1 - 1e-9: ws.append((s_, xa + s_ * e))
            ws.append((1.0, xb)); ws.sort(key=lambda w: w[0])
            for j in range(len(ws) - 1):
                if A.f(ws[j][1]) < tt < A.f(ws[j + 1][1]): break
            u, v = ws[j][1], ws[j + 1][1]
            lu, lv = A.lines_of(u), A.lines_of(v)
            if len(lu) != 1 or len(lv) != 1: bad += 1; continue
            pu, pv = A.piece(lu[0], A.param(lu[0], u)), A.piece(lv[0], A.param(lv[0], v))
            # breakpoint check
            if min(abs(A.B[lu[0]] - A.param(lu[0], u)).min(initial=9), 9) < 1e-9 or \
               min(abs(A.B[lv[0]] - A.param(lv[0], v)).min(initial=9), 9) < 1e-9: bad += 1
            key = ((lu[0], pu, 'L'), (lv[0], pv, 'H'))
            mid = (u + v) / 2; cell = tuple(int(np.sign(crs(mid - p, d))) for p, d in A.Ls)
        H.setdefault(key, []).append((cell, o))
    E = len(H); G = nx.Graph(); G.add_edges_from(H.keys())
    planar = nx.check_planarity(G)[0] if E else True
    mult = sum(1 for k, l in H.items() if len(l) > 2 or len(set(o for c, o in l)) < len(l))
    twocell = sum(1 for k, l in H.items() if len(set(c for c, o in l)) > 1)
    samel = sum(1 for k in H if k[0][0] == k[1][0] and k[0][1] != k[1][1])
    return dict(cut=cut, V=V, E=E, N=A.N, Lam=Lam, planar=planar, mult=mult, twocell=twocell,
                samel=samel, bad=bad, bound=4 * V - 8)
def thresholds(A, P, rng):
    fv = np.sort(np.array([A.f(x) for x in P])); u = np.unique(np.round(fv, 11))
    return [(x + y) / 2 for x, y in zip(u, u[1:])]
# ---------- arrangements ----------
def arr_pencil(rng, k):
    ang = np.sort(rng.uniform(0, np.pi, k)); O = rng.normal(size=2) * 0.2
    Ls = [(O, np.array([np.cos(a), np.sin(a)])) for a in ang]; return Arr(Ls, [O])
def arr_triangle(rng, vertex_only=False):
    while True:
        A_, B_, C_ = rng.normal(size=(3, 2))
        if abs(crs(B_ - A_, C_ - A_)) > 0.3: break
    Ls = [mkline(B_, C_), mkline(C_, A_), mkline(A_, B_)]
    Z = [A_, B_, C_]
    if vertex_only:
        # need nonobtuse: test
        def ang(X, Y, W): u, v = Y - X, W - X; return np.dot(u, v)
        if min(ang(A_, B_, C_), ang(B_, C_, A_), ang(C_, A_, B_)) < 0: return None
    else:
        Z += [foot(A_, Ls[0]), foot(B_, Ls[1]), foot(C_, Ls[2])]
    return Arr(Ls, Z)
def arr_par_trans(rng):
    a = rng.uniform(0.3, np.pi - 0.3); h = rng.uniform(0.5, 2)
    lb = (np.array([0., 0.]), np.array([1., 0.])); lc = (np.array([0., h]), np.array([1., 0.]))
    la = (np.array([0., 0.]), np.array([np.cos(a), np.sin(a)]))
    Xb = np.array([0., 0.]); Xc = np.array([h / np.tan(a), h])
    return Arr([la, lb, lc], [Xb, Xc, foot(Xc, lb), foot(Xb, lc)])
def arr_par(rng, k):
    hs = np.sort(rng.uniform(0, 2, k)); Ls = [(np.array([0., y]), np.array([1., 0.])) for y in hs]
    return Arr(Ls, [np.array([0., y]) for y in hs])
def sample_points(A, rng, n, mode):
    pts = []
    for i in range(n):
        li = int(rng.integers(len(A.Ls))); p, d = A.Ls[li]
        if mode == 0: s_ = rng.uniform(-3, 3)
        elif mode == 1:  # near breakpoints
            s_ = rng.choice(A.B[li]) + rng.normal() * 0.05
        else: s_ = np.sign(rng.normal()) * np.exp(rng.uniform(-4, 1.5))
        pts.append(p + s_ * d)
    return np.array(pts)
def make(rng):
    kind = int(rng.integers(6))
    if kind == 0: A = arr_pencil(rng, int(rng.integers(2, 5)))
    elif kind == 1: A = arr_triangle(rng)
    elif kind == 2:
        A = None
        while A is None: A = arr_triangle(rng, True)
    elif kind == 3: A = arr_par_trans(rng)
    elif kind == 4: A = arr_par(rng, int(rng.integers(2, 5)))
    else: A = arr_triangle(rng)
    return kind, A
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1])); its = int(sys.argv[2]); nmax = int(sys.argv[3])
    agg = dict(inst=0, cuts=0, viol_bound=0, viol_E=0, nonplanar=0, mult=0, twocell=0, samel=0, bad=0,
               maxratio=0.0, maxEratio=0.0, maxcutN=0.0)
    for it in range(its):
        kind, A = make(rng); n = int(rng.integers(7, nmax + 1))
        P = sample_points(A, rng, n, int(rng.integers(3)))
        c, t = hk(P); agg['inst'] += 1
        for tt in thresholds(A, P, rng):
            r = analyse(A, P, t, tt); agg['cuts'] += 1
            agg['viol_bound'] += r['cut'] > r['bound']; agg['viol_E'] += r['E'] > 2 * r['V'] - 4
            agg['nonplanar'] += not r['planar']
            for k_ in ('mult', 'twocell', 'samel', 'bad'): agg[k_] += r[k_]
            agg['maxratio'] = max(agg['maxratio'], r['cut'] / r['bound'])
            agg['maxEratio'] = max(agg['maxEratio'], r['E'] / (2 * r['V'] - 4))
            agg['maxcutN'] = max(agg['maxcutN'], r['cut'] / r['N'])
    print(agg)
