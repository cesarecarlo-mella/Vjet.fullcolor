import numpy as np, random
p=1048573
class Sketch:
    def __init__(s,K,seed=1): s.K=K; s.cache={}; s.seed=seed
    def __call__(s,v):
        out=np.zeros(s.K,dtype=np.int64); ca=s.cache
        for k,c in v.items():
            h=ca.get(k)
            if h is None:
                r=random.Random(hash((s.seed,k))); h=[(r.randrange(s.K),r.randrange(1,p)) for _ in range(3)]; ca[k]=h
            for j,a in h: out[j]+=c*a
        return out%p
class Basis:
    """incremental reduced row echelon basis mod p (p<2^20 so int64 matmul is exact for r<2^20)"""
    def __init__(s,K): s.K=K; s.R=np.zeros((0,K),dtype=np.int64); s.piv=[]
    def reduce(s,B):
        B=np.atleast_2d(B)%p
        if len(s.piv): B=(B-(B[:,s.piv]@s.R)%p)%p
        return B
    def add(s,B):
        B=s.reduce(B); new=0
        # rref of B
        rows=[]; pv=[]
        B=B[np.any(B,axis=1)]
        for c in range(s.K):
            if len(B)==0: break
            nz=np.nonzero(B[:,c])[0]
            if len(nz)==0: continue
            i=nz[0]; row=B[i]*pow(int(B[i,c]),p-2,p)%p
            B=np.delete(B,i,axis=0)
            if len(B):
                B=(B-(B[:,c][:,None]*row[None,:])%p)%p; B=B[np.any(B,axis=1)]
            # back-substitute into previous new rows
            for j in range(len(rows)):
                if rows[j][c]: rows[j]=(rows[j]-rows[j][c]*row)%p
            rows.append(row); pv.append(c)
        if rows:
            Nw=np.array(rows)
            if len(s.piv): s.R=(s.R-(s.R[:,pv]@Nw)%p)%p
            s.R=np.vstack([s.R,Nw]); s.piv+=pv
        return len(rows)
    def contains(s,B): return ~np.any(s.reduce(B),axis=1)
    def __len__(s): return len(s.piv)
