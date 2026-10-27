import sys
sys.argv=['fit_li.py','1','4','none']
exec(open('fit_li.py').read().split("T=S.copy()")[0])
OL=Span(); rA=rhoA(C['symbols']['Li2_2'])
for p in C['pairs']:
    vecs={0:p['x1'][1],1:p['x2'][1],2:p['1-x1'][1],3:p['1-x2'][1],4:p['1-x1x2'][1]}; OL.add(delta(subst(rA,vecs),2,2))
print("delta22 span of Li22 candidates:",len(OL.b))
g=sum(OL.copy().add(delta(v,2,2)) for v in [])
X=OL.copy(); g=sum(X.add(delta(v,2,2)) for v in data); print("data obstruction outside Li22 obstruction span:",g)
# classical-part check: data minus ... : test whether ker(delta) part of data is in Li4 span
# Li4 obstruction is zero; compute data-combos with zero obstruction: dimension = (data rank mod nothing) - obstruction rank
R=Span(); r=sum(R.add(v) for v in data); print("rank of data (rho):",r)
# full space of B2^B2 for 18-dim B2 is 153; also delta22 rank of all products Li2(x)*Li2(y)? (zero)
