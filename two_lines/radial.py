"""Instances on two crossing lines and radial cut statistics.

A point is (ray, dist) with rays 0=A (angle 0), 1=B (angle theta), 2=C (angle pi), 3=D (angle pi+theta).
mk(theta, ra, rb): signed positions ra on line 1 (A/C) and rb on line 2 (B/D) -> list of (ray, dist).
profile(theta, pts, tour): for every generic level t between consecutive distinct radii, returns
  (cut size, wedge cut edges per wedge [AB, BC, CD, DA], along cut edges per ray of the outer end).
maxcut(P, tour): largest prefix cut of the radial order (ties by index).
"""
import numpy as np


def mk(th, ra, rb):
    pts = []
    for x in ra:
        pts.append((0 if x > 0 else 2, abs(x)))
    for x in rb:
        pts.append((1 if x > 0 else 3, abs(x)))
    return pts


def coords(th, pts):
    ang = [0, th, np.pi, np.pi + th]
    return np.array([[d * np.cos(ang[r]), d * np.sin(ang[r])] for r, d in pts])


def profile(th, pts, tour):
    n = len(pts)
    d = np.array([p[1] for p in pts])
    ray = [p[0] for p in pts]
    lev = np.sort(d)
    res = []
    E = [(tour[i], tour[(i + 1) % n]) for i in range(n)]
    wedge_of = {(0, 1): 0, (1, 0): 0, (1, 2): 1, (2, 1): 1, (2, 3): 2, (3, 2): 2, (3, 0): 3, (0, 3): 3}
    for k in range(n - 1):
        if lev[k + 1] - lev[k] < 1e-12:
            continue
        t = (lev[k] + lev[k + 1]) / 2
        om = [0] * 4
        al = [0] * 4
        for a, b in E:
            ia, ib = d[a] < t, d[b] < t
            if ia == ib:
                continue
            o = b if ia else a
            if ray[a] % 2 == ray[b] % 2:
                al[ray[o]] += 1
            else:
                om[wedge_of[(ray[a], ray[b])]] += 1
        res.append((sum(om) + sum(al), om, al))
    return res


def maxcut(P, tour):
    n = len(P)
    f = np.sqrt((np.asarray(P) ** 2).sum(1))
    order = np.argsort(f, kind='stable')
    rank = np.empty(n, int)
    rank[order] = np.arange(n)
    diff = np.zeros(n + 1, int)
    for i in range(n):
        a, b = tour[i], tour[(i + 1) % n]
        lo, hi = sorted((rank[a], rank[b]))
        diff[lo + 1] += 1
        diff[hi + 1] -= 1
    return int(np.cumsum(diff)[1:n].max())
