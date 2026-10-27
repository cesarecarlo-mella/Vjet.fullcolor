import numpy as np
p=1048573
class BasisF:
    """incremental rref mod p; reductions via float64 BLAS (exact: entries<2^20, rank<2^13)"""
    def __init__(s,K): s.K=K; s.R=np.zeros((0,K),dtype=np.int64); s.Rf=np.zeros((0,K)); s.piv=[]
    def reduce(s,B):
        B=np.atleast_2d(B)%p
        if len(s.piv):
            Bp=B[:,s.piv].astype(np.float64)
            # split Rf product to keep exactness: sum over r terms of (<2^20 * <2^20) < 2^53 needs r<2^13
            acc=np.zeros(B.shape,dtype=np.int64)
            for i in range(0,len(s.piv),8000):
                acc=(acc+np.rint(Bp[:,i:i+8000]@s.Rf[i:i+8000]).astype(np.int64))%p
            B=(B-acc)%p
        return B
    def add(s,B):
        B=s.reduce(B); B=B[np.any(B,axis=1)]; rows=[]; pv=[]
        while len(B):
            nzc=np.nonzero(B[0])[0]
            if len(nzc)==0: B=B[1:]; continue
            c=nzc[0]; row=B[0]*pow(int(B[0,c]),p-2,p)%p; B=B[1:]
            if len(B):
                f=B[:,c]; B=(B-(f[:,None]*row[None,:])%p)%p; B=B[np.any(B,axis=1)]
            for j in range(len(rows)):
                if rows[j][c]: rows[j]=(rows[j]-rows[j][c]*row)%p
            rows.append(row); pv.append(c)
        if rows:
            Nw=np.array(rows)
            if len(s.piv):
                s.R=(s.R-((s.R[:,pv].astype(np.float64)@Nw.astype(np.float64)).round().astype(np.int64))%p)%p
            s.R=np.vstack([s.R,Nw]); s.piv+=pv; s.Rf=s.R.astype(np.float64)
        return len(rows)
    def copy(s):
        n=BasisF(s.K); n.R=s.R.copy(); n.Rf=s.Rf.copy(); n.piv=list(s.piv); return n
    def __len__(s): return len(s.piv)
