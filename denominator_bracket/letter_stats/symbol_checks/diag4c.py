import sys
sys.argv=['fit_li.py','1','4','none']
src=open('fit_li.py').read().split("T=S.copy()")[0]
src=src.replace("data.append({k:v%P for k,v in out.items() if v%P})","data.append((amp,{k:v%P for k,v in out.items() if v%P}))")
exec(src)
D=S.copy(); rA=rhoA(C['symbols']['Li2_2'])
for p in C['pairs']:
    vecs={0:p['x1'][1],1:p['x2'][1],2:p['1-x1'][1],3:p['1-x2'][1],4:p['1-x1x2'][1]}; D.add(subst(rA,vecs))
bad=collections.Counter()
for amp,v in data:
    if D.copy().add(v): bad[amp]+=1
print("amplitudes whose bare w4 (even part) is NOT in span Li4+Li22:",dict(bad))
