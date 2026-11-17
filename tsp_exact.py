"""Exact and local TSP helpers used by all experiments.

hk(P)          Held-Karp optimum of the points P (n x 2 array), returns (cost, tour as list of indices).
two_opt(P, t)  2-opt local search from tour t (strict improvements only), returns a 2-optimal tour.

If the shared libraries libhk.so / libtwo.so (built from hk.c / two_opt.c, see README) are present next to this
file they are used; otherwise pure-Python fallbacks are used (fine for n <= 12, slow beyond).
"""
import ctypes
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    try:
        return ctypes.CDLL(os.path.join(_HERE, name))
    except OSError:
        return None


_hk = _load('libhk.so')
if _hk is not None:
    _hk.hk.restype = ctypes.c_double
_two = _load('libtwo.so')


def dist_matrix(P):
    P = np.asarray(P, dtype=np.float64)
    return np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))


def hk_d(D):
    """Held-Karp on a distance matrix D; returns (cost, tour)."""
    n = len(D)
    D = np.ascontiguousarray(D, dtype=np.float64)
    if _hk is not None:
        t = (ctypes.c_int * n)()
        c = _hk.hk(n, D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), t)
        return c, list(t)
    return _hk_py(D)


def hk(P):
    return hk_d(dist_matrix(P))


def _hk_py(D):
    n = len(D)
    m = n - 1
    full = 1 << m
    INF = float('inf')
    dp = np.full((full, m), INF)
    par = np.full((full, m), -1, dtype=np.int64)
    for j in range(m):
        dp[1 << j, j] = D[0, j + 1]
    Dm = D[1:, 1:]
    for S in range(1, full):
        row = dp[S]
        for j in range(m):
            v = row[j]
            if v == INF or not (S >> j) & 1:
                continue
            for k in range(m):
                if (S >> k) & 1:
                    continue
                T = S | (1 << k)
                w = v + Dm[j, k]
                if w < dp[T, k]:
                    dp[T, k] = w
                    par[T, k] = j
    S = full - 1
    best, bj = INF, -1
    for j in range(m):
        w = dp[S, j] + D[j + 1, 0]
        if w < best:
            best, bj = w, j
    tour = [0] * n
    pos = n - 1
    j = bj
    while j >= 0:
        tour[pos] = j + 1
        pos -= 1
        pj = int(par[S, j])
        S &= ~(1 << j)
        j = pj
    return float(best), tour


def two_opt(P, tour):
    D = np.ascontiguousarray(dist_matrix(P))
    n = len(D)
    if _two is not None:
        arr = (ctypes.c_int * n)(*tour)
        _two.twoopt(n, D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), arr)
        return list(arr)
    t = list(tour)
    improved = True
    while improved:
        improved = False
        for i in range(n - 1):
            for j in range(i + 2, n):
                a, b, c, e = t[i], t[i + 1], t[j], t[(j + 1) % n]
                if e == a:
                    continue
                if D[a, c] + D[b, e] - D[a, b] - D[c, e] < -1e-12:
                    t[i + 1:j + 1] = t[i + 1:j + 1][::-1]
                    improved = True
    return t


def tour_length(P, tour):
    P = np.asarray(P, float)
    return float(sum(np.linalg.norm(P[tour[i]] - P[tour[(i + 1) % len(tour)]]) for i in range(len(tour))))
