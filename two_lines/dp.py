"""Radial sweep DP for Euclidean TSP on two crossing lines.
Points: list of (ray, dist) with rays 0=A(angle 0),1=B(theta),2=C(pi),3=D(pi+theta); dist>=0.
Along edges are paid incrementally (cost additive along a ray); only wedge cut edges carry an
identity (their inner endpoint). caps: max pending per wedge (cw), max along stubs per ray (ca),
max through-O along edges per line (co)."""
import math, itertools
WEDGES=[(0,1),(1,2),(2,3),(3,0)]
def wedges_of(r): return [w for w,(a,b) in enumerate(WEDGES) if r in (a,b)]
def other(w,r): a,b=WEDGES[w]; return b if r==a else a

def solve(theta, pts, cw=2, ca=2, co=2, ct=4, stats=None):
    n=len(pts); ang=[0,theta,math.pi,math.pi+theta]
    X=[(d*math.cos(ang[r]), d*math.sin(ang[r])) for r,d in pts]
    def dist(i,j): return math.hypot(X[i][0]-X[j][0],X[i][1]-X[j][1])
    order=sorted(range(n),key=lambda i:pts[i][1])
    # state: (W, R, lastd) where W = tuple over wedges of sorted tuple of (pid,label); R = tuple over rays of sorted tuple of labels
    # lastd per ray stored separately as it is determined by processed prefix. labels: each open path has a label, appears exactly twice.
    def canon(W,R):
        mp={}; c=[0]
        def L(x):
            if x not in mp: mp[x]=c[0]; c[0]+=1
            return mp[x]
        W2=tuple(tuple(sorted((p,L(l)) for p,l in sorted(w,key=lambda z:(z[0],z[1])))) for w in W)
        R2=tuple(tuple(sorted(L(l) for l in rr)) for rr in R)
        # relabel needs deterministic order: relabel by traversal order of W then R
        return W2,R2
    def canon2(W,R):
        # deterministic canonical labelling: order labels by first appearance scanning W (sorted by pid) then R
        seq=[]
        for w in W:
            for p,l in sorted(w): seq.append(l)
        for rr in R:
            for l in sorted(rr): seq.append(l)
        mp={}
        for l in seq:
            if l not in mp: mp[l]=len(mp)
        W2=tuple(tuple(sorted((p,mp[l]) for p,l in w)) for w in W)
        R2=tuple(tuple(sorted(mp[l] for l in rr)) for rr in R)
        return (W2,R2)
    # initial states: through-O along edges on line 1 (rays 0,2) and line 2 (rays 1,3)
    init={}
    for k1 in range(co+1):
        for k2 in range(co+1):
            R=[[] for _ in range(4)]; lab=0
            for _ in range(k1): R[0].append(lab); R[2].append(lab); lab+=1
            for _ in range(k2): R[1].append(lab); R[3].append(lab); lab+=1
            if any(len(r)>ca for r in R): continue
            init[canon2(((),(),(),()),R)]=0.0
    states=init; lastd=[0.0]*4; best_final=math.inf; nproc=0
    # remaining points per ray after processing index idx
    rem=[[0]*4 for _ in range(n+1)]
    for k in range(n-1,-1,-1):
        rem[k]=list(rem[k+1]); rem[k][pts[order[k]][0]]+=1
    for idx,v in enumerate(order):
        r,d=pts[v]; last=(idx==n-1)
        new={}
        for (W,R),cost in states.items():
            # pay along stubs on ray r
            c0=cost+len(R[r])*(d-lastd[r])
            # enumerate the two edge ends of v
            # options: ('close',l) along stub on r; ('open',) new along stub; ('match',w,p,l); ('new',w)
            opts=[]
            for l in set(R[r]): opts.append(('close',l))
            opts.append(('open',)); 
            for w in wedges_of(r):
                o=other(w,r)
                seen=set()
                for p,l in W[w]:
                    if pts[p][0]==o and (p,l) not in seen: seen.add((p,l)); opts.append(('match',w,p,l))
                opts.append(('new',w))
            for a,b in itertools.combinations_with_replacement(range(len(opts)),2):
                oa,ob=opts[a],opts[b]
                if a==b and oa[0] in ('close','match'):
                    # same stub twice only possible if multiplicity>=2 - handle: close with same label twice means two stubs of same path (cycle)
                    pass
                Wl=[list(w) for w in W]; Rl=[list(x) for x in R]; c=c0; ok=True
                ends=[]  # labels joined at v
                newlabs=[]
                for o in (oa,ob):
                    if o[0]=='close':
                        if o[1] not in Rl[r]: ok=False;break
                        Rl[r].remove(o[1]); ends.append(o[1])
                    elif o[0]=='match':
                        _,w,p,l=o
                        if (p,l) not in Wl[w]: ok=False;break
                        Wl[w].remove((p,l)); ends.append(l); c+=dist(p,v)
                    else:
                        newlabs.append(o)
                if not ok: continue
                if oa[0]=='match' and ob[0]=='match' and oa[2]==ob[2]: continue  # double edge
                # combine: v joins up to two paths
                if len(ends)==2:
                    l1,l2=ends
                    if l1==l2:
                        # closing a cycle: only allowed at the very end with nothing else open
                        if last and all(len(w)==0 for w in Wl) and all(len(x)==0 for x in Rl):
                            best_final=min(best_final,c)
                        continue
                    # merge l2 into l1
                    Wl=[[(p,(l1 if l==l2 else l)) for p,l in w] for w in Wl]
                    Rl=[[(l1 if l==l2 else l) for l in x] for x in Rl]
                elif len(ends)==1:
                    l1=ends[0]
                    o=newlabs[0]
                    if o[0]=='open': Rl[r].append(l1)
                    else: Wl[o[1]].append((v,l1))
                else:
                    l1=10**6+idx
                    for o in newlabs:
                        if o[0]=='open': Rl[r].append(l1)
                        else: Wl[o[1]].append((v,l1))
                if any(len(w)>cw for w in Wl) or any(len(x)>ca for x in Rl) or sum(len(w) for w in Wl)>ct: continue
                if last: continue
                rm=rem[idx+1]
                if any(len(Rl[q])>0 and rm[q]==0 for q in range(4)): continue
                if any(rm[other(w,p_r)]==0 for w in range(4) for p_r in [pts[p][0] for p,l in Wl[w]]): continue
                tot=sum(len(w) for w in Wl)+sum(len(x) for x in Rl)
                if tot>2*sum(rm): continue
                key=canon2(Wl,Rl)
                if c<new.get(key,math.inf): new[key]=c
        states=new; lastd[r]=d
        if stats is not None: stats.append(len(states))
    return best_final
