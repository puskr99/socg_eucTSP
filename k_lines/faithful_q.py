"""Exact (rational) faithful-net closure for line arrangements; same rule as faithful_exact.py but in Fractions,
so ties (Thales / altitude-theorem equalities) are decided exactly.  Lines: p + t d with rational p, d."""
import itertools, argparse, json, math, random
from fractions import Fraction as F

def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def dot(a, b): return a[0]*b[0] + a[1]*b[1]
def crs(a, b): return a[0]*b[1] - a[1]*b[0]

def crossing(L1, L2):
    (p1, d1), (p2, d2) = L1, L2
    den = crs(d1, d2)
    if den == 0: return None
    s = crs(sub(p2, p1), d2) / den
    return (p1[0] + s*d1[0], p1[1] + s*d1[1])

def foot(z, L):
    p, d = L
    a = dot(sub(z, p), d) / dot(d, d)
    return (p[0] + a*d[0], p[1] + a*d[1])

def closure(L, max_pts=3000, max_rounds=60, verbose=False):
    k = len(L); pts = []; lab = []
    def add(z):
        if z in seen: return False
        seen.add(z); pts.append(z)
        lab.append({i for i in range(k) if crs(sub(z, L[i][0]), L[i][1]) == 0}); return True
    seen = set()
    for a, b in itertools.combinations(range(k), 2):
        c = crossing(L[a], L[b])
        if c is not None: add(c)
    for i in range(k):
        if not any(i in s for s in lab): add(L[i][0])
    for r in range(max_rounds):
        newf = []
        for i in range(k):
            p, d = L[i]; dd = dot(d, d)
            par = [dot(sub(z, p), d) / dd for z in pts]
            ton = sorted(par[j] for j in range(len(pts)) if i in lab[j])
            for j in range(len(pts)):
                if i in lab[j]: continue
                a = par[j]; w = sub(pts[j], p); H = dot(w, w)/dd - a*a
                if a in ton: continue
                left = [t for t in ton if t < a]; right = [t for t in ton if t > a]
                lo = max(((t + a)/2 + H/(2*(a - t)) for t in left), default=None)
                hi = min(((t + a)/2 - H/(2*(t - a)) for t in right), default=None)
                if lo is None or hi is None or lo < hi:
                    newf.append((foot(pts[j], L[i]), i, j))
        if not newf:
            return dict(ok=True, size=len(pts), rounds=r, per_line=[sum(i in s for s in lab) for i in range(k)])
        added = sum(add(z) for z, _, _ in newf)
        if verbose: print('round', r, 'violations', len(newf), 'added', added, 'size', len(pts))
        if added == 0: return dict(ok=False, size=len(pts), rounds=r, why='stuck')
        if len(pts) > max_pts: return dict(ok=False, size=len(pts), rounds=r, why='too big')
    return dict(ok=False, size=len(pts), rounds=max_rounds, why='rounds')

def q(x, den=10**6): return F(round(x*den), den)

def arrangement(k, rnd, kind, eps):
    L = []
    for i in range(k):
        if kind == 'random':
            th = rnd.uniform(0, math.pi); p = (rnd.uniform(-1, 1), rnd.uniform(-1, 1))
        elif kind == 'flat':
            th = rnd.uniform(-eps, eps); p = (rnd.uniform(-1, 1), rnd.uniform(-eps, eps))
        elif kind == 'nearconc':
            th = rnd.uniform(0, math.pi); p = (rnd.uniform(-eps, eps), rnd.uniform(-eps, eps))
        L.append(((q(p[0]), q(p[1])), (q(math.cos(th)), q(math.sin(th)))))
    return L

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', default='3,4'); ap.add_argument('--kinds', default='random,flat,nearconc')
    ap.add_argument('--eps', default='0.1,0.01,0.001'); ap.add_argument('--trials', type=int, default=15)
    ap.add_argument('--seed', type=int, default=5); ap.add_argument('--out', default='')
    a = ap.parse_args(); rnd = random.Random(a.seed); res = []
    for kind in a.kinds.split(','):
        for eps in map(float, a.eps.split(',')):
            if kind == 'random' and eps != float(a.eps.split(',')[0]): continue
            for k in map(int, a.k.split(',')):
                rows = []
                for t in range(a.trials):
                    r = closure(arrangement(k, rnd, kind, eps)); r.update(kind=kind, k=k, eps=eps, trial=t); rows.append(r)
                ok = [r for r in rows if r['ok']]; sz = sorted(r['size'] for r in ok)
                print(f"{kind:8s} eps={eps:<6g} k={k}: faithful {len(ok)}/{len(rows)} |Z| max {sz[-1] if sz else '-'} sizes {sz}"
                      f" fails {[(r['why'], r['size']) for r in rows if not r['ok']]}", flush=True)
                res += rows
    if a.out: json.dump(res, open(a.out, 'w'))
