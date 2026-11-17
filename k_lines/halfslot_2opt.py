"""The checks of halfslot_check.py on 2-optimal tours (2-opt from a nearest-neighbour tour) with n = 100, 200, 400;
Theorem 21 holds for every 2-optimal tour.  Uses at most 80 random thresholds per instance.

usage: python halfslot_2opt.py SEED INSTANCES     (paper: 320 instances)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halfslot_check import *
from tsp_exact import two_opt


def two(P, rng):
    n = len(P); D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
    t = [0]; left = set(range(1, n))
    while left:
        j = min(left, key=lambda j: D[t[-1], j]); t.append(j); left.remove(j)
    return two_opt(P, t)


rng=np.random.default_rng(int(sys.argv[1])); its=int(sys.argv[2])
agg=dict(inst=0,cuts=0,viol_bound=0,viol_E=0,nonplanar=0,mult=0,twocell=0,samel=0,bad=0,maxratio=0.0,maxEratio=0.0,maxcutN=0.0)
for it in range(its):
    kind,A=make(rng); n=int(rng.choice([100,200,400]))
    P=sample_points(A,rng,n,int(rng.integers(3))); t=two(P,rng); agg['inst']+=1
    ths=thresholds(A,P,rng); ths=[ths[i] for i in sorted(rng.choice(len(ths),min(80,len(ths)),replace=False))]
    for tt in ths:
        r=analyse(A,P,t,tt); agg['cuts']+=1
        agg['viol_bound']+=r['cut']>r['bound']; agg['viol_E']+=r['E']>2*r['V']-4; agg['nonplanar']+=not r['planar']
        for k_ in ('mult','twocell','samel','bad'): agg[k_]+=r[k_]
        agg['maxratio']=max(agg['maxratio'],r['cut']/r['bound']); agg['maxEratio']=max(agg['maxEratio'],r['E']/(2*r['V']-4)); agg['maxcutN']=max(agg['maxcutN'],r['cut']/r['N'])
print(agg)
