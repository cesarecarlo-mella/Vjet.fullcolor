import sys
sys.argv=['fit_li.py','1','4','none']
src=open('fit_li.py').read().split("T=S.copy()")[0]
exec(src)
# obstruction rank of even-only data
O=Span(); [O.add(delta(v,2,2)) for v in data]; print("delta22 rank of even-only bare w4 data:",len(O.b))
# residual mod (Li4 + Li22)
C2=C; D=S.copy(); rA=rhoA(C['symbols']['Li2_2'])
for p in C['pairs']:
    vecs={0:p['x1'][1],1:p['x2'][1],2:p['1-x1'][1],3:p['1-x2'][1],4:p['1-x1x2'][1]}; D.add(subst(rA,vecs))
before=set(D.b)
for v in data: D.add(v)
new=[D.b[k] for k in D.b if k not in before]
print("residual directions:",len(new))
for r in new:
    dd=delta(r,2,2); print(" residual: #terms",len(r)," delta22 terms:",len(dd)," letters used:",sorted({l for w in r for l in w}))
    print("   sample:",list(r.items())[:6])
