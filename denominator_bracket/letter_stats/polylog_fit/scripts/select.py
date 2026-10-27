# greedy selection of non-classical functions, working modulo the classical span
import numpy as np
def quotient_coords(Bcl, cands, data, K):
    """reduce everything mod classical, then express in a common basis of the quotient"""
    Cr=Bcl.reduce(np.array(cands)); Dr=Bcl.reduce(np.array(data))
    Q=Basis(K); Q.add(np.vstack([Cr,Dr]))
    # coordinates: entries at pivot columns of rref basis Q
    return Cr[:,Q.piv]%p, Dr[:,Q.piv]%p, len(Q.piv)
def greedy(Cq, Dq, order):
    n=Cq.shape[1]; S=Basis(n); SD=Basis(n); SD.add(Dq)
    chosen=[]; unc=len(SD)
    while unc>0:
        found=None
        for i in order:
            if i in chosen: continue
            v=Cq[i:i+1]
            if SD.contains(v)[0] and not S.contains(v)[0]: found=i; break
        if found is None: break
        chosen.append(found); S.add(Cq[found:found+1]); SD.add(Cq[found:found+1])
        T=Basis(n); T.R=S.R.copy(); T.piv=list(S.piv); unc=T.add(Dq)
    return chosen, unc
def support_basis(Cq, Dq, order):
    """independent candidates in 'order' spanning everything; coordinates of data -> support"""
    n=Cq.shape[1]; S=Basis(n); basis=[]
    for i in order:
        if S.add(Cq[i:i+1]): basis.append(i)
    M=Cq[basis]                     # (r x n), invertible on its span
    # solve c M = d for each data row: use rref of [M^T | D^T]
    A=np.concatenate([M.T, Dq.T],axis=1)%p
    R,piv=rref(A) if False else (None,None)
    # Gaussian elimination: express D rows in basis rows
    Bm=Basis(n+len(basis)); aug=np.concatenate([M,np.eye(len(basis),dtype=np.int64)],axis=1); Bm.add(aug)
    red=Bm.reduce(np.concatenate([Dq,np.zeros((len(Dq),len(basis)),dtype=np.int64)],axis=1))
    coef=red[:,n:]                  # = -coefficients
    used=sorted({basis[j] for j in np.nonzero(np.any(coef,axis=0))[0]})
    return used, basis
