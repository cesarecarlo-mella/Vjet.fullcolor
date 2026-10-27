exec(open('fast.py').read()); exec(open('sk3.py').read())
import pickle
t0=time.time(); K=8000; sk=R1(K,6,seed=21)
ALG={a for a in by if a[:6] in ('C1F0Z0','C7F0Z0','C9F0Z0')}
B=BasisF(K)
W,Cf=abstract_rho(C['symbols']['Li6'])
B.add(np.array([sk.abstract(W,Cf,{0:[m(x) for x in z],2:[m(x) for x in o]}) for z,o in E['args']])); r6=len(B)
W,Cf=abstract_rho(C['symbols']['Li4_2']); rows=[]
for x1,x2,o1,o2,o12 in E['pairs']:
    rows.append(sk.abstract(W,Cf,{0:[m(x) for x in x1],1:[m(x) for x in x2],2:[m(x) for x in o1],3:[m(x) for x in o2],4:[m(x) for x in o12]}))
    if len(rows)==500: B.add(np.array(rows)); rows=[]
if rows: B.add(np.array(rows))
r62=len(B); print("Li6",r6,"+Li42",r62-r6,f"{time.time()-t0:.0f}s",flush=True)
T=json.load(open('triples.json')); used=json.load(open('li222_used.json'))
S3=json.load(open('Li2_2_2.json')); W,Cf=abstract_rho(S3)
ext=lambda v:[m(x) for x in v]+[0,0]
B.add(np.array([sk.abstract(W,Cf,{a:ext(v) for a,v in enumerate(T[i]['vecs'])}) for i in used]))
r3=len(B); print("+Li222 (466 triples)",r3-r62,f"{time.time()-t0:.0f}s",flush=True)
dv=data_vectors('',6,even_only=False,maxpts=10)
Dn=np.array([sk.data(v) for a,v in dv if a not in ALG]); Da=np.array([sk.data(v) for a,v in dv if a in ALG])
U=B.copy(); print("not covered (no spurious): non-alg",U.add(Dn),"alg",U.add(Da),flush=True)
# spurious-free depth-2 combinations
args,Wg,res,newf=pickle.load(open('spur.pkl','rb')); LETR=C['letters']
EXTRA=[9,10,17,18,19,20]+[('N',i) for i in range(len(newf))]
ALLL=BX+EXTRA; AI={l:i for i,l in enumerate(ALLL)}; NB=len(BX)
KS=400; R2=np.random.default_rng(77).integers(1,p,size=(6,len(ALLL),KS),dtype=np.int64); R2b=R2.copy(); R2b[:,NB:,:]=0
def G_of(vec,R,baseonly=False):
    out=np.zeros((R.shape[0],R.shape[2]),dtype=np.int64)
    for l,c in vec.items():
        i=AI[l]
        if baseonly and i>=NB: continue
        out=(out+(int(c)%p)*R[:,i,:])%p
    return out
def sk_abs(Wd,Cf,vecs,R,baseonly=False):
    G=np.stack([G_of(vecs[a],R,baseonly) for a in range(5)])
    prod=G[Wd[:,0],0]
    for i in range(1,6): prod=prod*G[Wd[:,i],i]%p
    return (prod*Cf[:,None]%p).sum(axis=0)%p
Rb=np.transpose(sk.R,(0,1,2))   # (6, NB, K) base letters for main sketch
argvec=lambda a:{l:e for l,e in zip(LETR,a[1]) if e}
amap={(a[0],tuple(a[1])):(b[0],tuple(b[1])) for a,b in C['args']}
fns={fn:abstract_rho(C['symbols'][fn]) for fn in ['Li4_2','Li3_3','Li2_4','Li5_1','Li1_5']}
def rref_small(M):
    M=M.copy()%p; n,K_=M.shape; piv=[]; r=0
    for c in range(K_):
        if r==n: break
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]; M[[r,i]]=M[[i,r]]; M[r]=M[r]*pow(int(M[r,c]),p-2,p)%p
        f=M[:,c].copy(); f[r]=0; nzr=np.nonzero(f)[0]
        if len(nzr): M[nzr]=(M[nzr]-(f[nzr,None]*M[r][None,:])%p)%p
        piv.append(c); r+=1
    return M[:r],piv
groups=sorted(Wg.items(),key=lambda kv:-len(kv[1]))[1:1501]
combos=[]; gi=0
for w,prs in groups:
    gi+=1
    v4=res[w]; meta=[]; rows_sp=[]
    for i,j in prs:
        x,y=args[i],args[j]; ox,oy=amap[x],amap[y]
        vecs={0:argvec(x),1:argvec(y),2:argvec(ox),3:argvec(oy),4:v4}
        for fn,(Wd,Cf_) in fns.items():
            rows_sp.append((sk_abs(Wd,Cf_,vecs,R2)-sk_abs(Wd,Cf_,vecs,R2b))%p); meta.append((vecs,fn))
    Msp=np.array(rows_sp); m_=len(Msp)
    if m_>KS-50: continue
    Rr,piv=rref_small(np.concatenate([Msp,np.eye(m_,dtype=np.int64)],axis=1))
    nulls=[r_[KS:] for r_,pc in zip(Rr,piv) if pc>=KS]
    if nulls:
        need=sorted({int(k) for c in nulls for k in np.nonzero(c)[0]})
        Mb={k:sk_abs(*fns[meta[k][1]],meta[k][0],Rb,baseonly=True) for k in need}
        for c in nulls:
            v=np.zeros(K,dtype=np.int64)
            for k in np.nonzero(c)[0]: v=(v+int(c[k])*Mb[int(k)])%p
            combos.append(v)
    if len(combos)>=1000:
        B.add(np.array(combos)); combos=[]
        U=B.copy(); print(gi,"groups: total rank",len(B),"(spurious adds",len(B)-r3,"); not covered non-alg",U.add(Dn),"alg",U.add(Da),f"{time.time()-t0:.0f}s",flush=True)
if combos: B.add(np.array(combos))
U=B.copy(); print("FINAL total rank",len(B),"K",K,"; not covered non-alg",U.add(Dn),"alg",U.add(Da),f"{time.time()-t0:.0f}s",flush=True)
