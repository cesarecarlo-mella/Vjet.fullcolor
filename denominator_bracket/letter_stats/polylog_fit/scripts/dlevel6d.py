exec(open('dlevel6.py').read().split("T=json.load")[0])
for fn in ['Li3_3','Li2_4','Li5_1','Li1_5']:
    for i in range(0,len(pv),500): B.add(np.array(cand(C['symbols'][fn],pv[i:i+500])))
    U=Basis(2*K); U.R=B.R.copy(); U.piv=list(B.piv)
    print("+",fn,"ext pairs: delta span",len(B),"; not covered non-alg",U.add(Dn),"alg",U.add(Da),f"{time.time()-t0:.0f}s",flush=True)
r2=len(B)
T=json.load(open('triples.json')); S3=json.load(open('Li2_2_2.json'))
ext=lambda v:[m(x) for x in v]+[0,0]
for i in range(0,len(T),400):
    B.add(np.array(cand(S3,[{a:ext(v) for a,v in enumerate(t['vecs'])} for t in T[i:i+400]])))
    U=Basis(2*K); U.R=B.R.copy(); U.piv=list(B.piv)
    print(i+400,"Li222 triples: span beyond depth-2",len(B)-r2,"; not covered non-alg",U.add(Dn),"alg",U.add(Da),f"{time.time()-t0:.0f}s",flush=True)
