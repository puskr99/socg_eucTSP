"""Optional: exact TSP via an external Concorde binary (not included; set CONCORDE_BIN to its path).

Points are written as TSPLIB EUC_2D with coordinates multiplied by `scale` (default 1e7); Concorde rounds each
distance to an integer, so the tour is optimal for the rounded instance.  On a crash, a non-zero exit or a timeout the
solve is retried with the next random seed.  Returns (Euclidean length, tour starting at 0).
"""
import math, os, shutil, subprocess, tempfile

BIN = os.environ.get('CONCORDE_BIN', 'concorde')


def _attempt(P, scale, seed, timeout):
    n = len(P); wd = tempfile.mkdtemp(prefix='concorde_')
    try:
        with open(os.path.join(wd, 'p.tsp'), 'w') as f:
            f.write('NAME: p\nTYPE: TSP\nDIMENSION: %d\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION\n' % n)
            for i, (x, y) in enumerate(P):
                f.write('%d %.3f %.3f\n' % (i + 1, x * scale, y * scale))
            f.write('EOF\n')
        try:
            out = subprocess.run([BIN, '-s', str(seed), '-o', 'p.sol', 'p.tsp'], cwd=wd, capture_output=True,
                                 text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return None
        if out.returncode != 0 or 'FATAL ERROR' in out.stdout:
            return None
        tok = open(os.path.join(wd, 'p.sol')).read().split()
        assert int(tok[0]) == n
        return [int(v) for v in tok[1:n + 1]]
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def solve(P, scale=1e7, seed=1, timeout=None, tries=4):
    n = len(P)
    for a in range(tries):
        t = _attempt(P, scale, seed + a, timeout)
        if t is not None:
            break
    else:
        raise RuntimeError('concorde failed %d times' % tries)
    i = t.index(0); t = t[i:] + t[:i]
    return sum(math.dist(P[t[i]], P[t[(i + 1) % n]]) for i in range(n)), t
