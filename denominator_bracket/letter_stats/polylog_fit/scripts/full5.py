exec(open('fast.py').read())
t0=time.time(); K=1500; sk=R1(K,5,seed=1)
ALG={a for a in by if a[:6] in ('C1F0Z0','C7F0Z0','C9F0Z0')}
B=Basis(K)
W,Cf=abstract_rho(C['symbols']['Li5'])
B.add(np.array([sk.abstract(W,Cf,{0:[m(x) for x in z],2:[m(x) for x in o]}) for z,o in E['args']])); r5=len(B); print("Li5 (ext args)",r5,f"{time.time()-t0:.0f}s",flush=True)
W,Cf=abstract_rho(C['symbols']['Li3_2']); rows=[]
for x1,x2,o1,o2,o12 in E['pairs']:
    rows.append(sk.abstract(W,Cf,{0:[m(x) for x in x1],1:[m(x) for x in x2],2:[m(x) for x in o1],3:[m(x) for x in o2],4:[m(x) for x in o12]}))
B.add(np.array(rows)); r52=len(B); print("+Li32 (ext pairs)",r52-r5,f"{time.time()-t0:.0f}s",flush=True)
dv=data_vectors('',5,even_only=False)
Dn=np.array([sk.data(v) for a,v in dv if a not in ALG]); Da=np.array([sk.data(v) for a,v in dv if a in ALG])
U=Basis(K); U.R=B.R.copy(); U.piv=list(B.piv); print("data not covered (Li5+depth2): non-alg",U.add(Dn),"alg",U.add(Da),flush=True)
T=json.load(open('triples.json')); D3=json.load(open('depth3_w5.json'))
W,Cf=abstract_rho(D3['Li2_2_1']); rows=[]; 
ext=lambda v:[m(x) for x in v]+[0,0]
for t in T: rows.append(sk.abstract(W,Cf,{a:ext(v) for a,v in enumerate(t['vecs'])}))
B.add(np.array(rows)); print("+Li221 (rational triples)",len(B)-r52,f"{time.time()-t0:.0f}s",flush=True)
U=Basis(K); U.R=B.R.copy(); U.piv=list(B.piv); print(">>> w5 data not covered: non-alg",U.add(Dn),"alg",U.add(Da),f"{time.time()-t0:.0f}s",flush=True)
