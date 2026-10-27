import re,glob,collections,pickle,os
D=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CW={'':0,'Pi^2':2,'Zeta[3]':3,'Pi^4':4,'I':0,'Zeta[5]':5,'Pi^2*Zeta[3]':5,'Zeta[3]^2':6,'Pi^6':6}
rows=[]  # (amp, const, word tuple, den, coef, poly)
for f in sorted(glob.glob(D+'/*.wl')):
    amp=os.path.basename(f)[:-3]
    L=open(f).read().splitlines()
    i=[k for k,l in enumerate(L) if l.startswith('+eps^(0)*')][0]
    buf=[]
    for l in L[i+1:]:
        m=re.match(r'^\+\(([-0-9/]+)\)\*\((.*)\)\*\((.*)\)$',l)
        if m:
            c=m.group(3); w=re.search(r'II\[([0-9,]*)\]',c).group(1)
            const=c[:c.index('II[')].rstrip('*')
            word=tuple(int(x) for x in w.split(',')) if w else ()
            buf.append((const,word,m.group(1),m.group(2)))
        elif l.startswith(')/('):
            den=l[3:-1]
            for b in buf: rows.append((amp,b[0],b[1],den,b[2],b[3]))
            buf=[]
    assert not buf, amp
pickle.dump(rows,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'rows.pkl'),'wb'))
print(len(rows))
tw=collections.Counter(len(r[2])+CW[r[1]] for r in rows); print("total weight",tw)
print("const x II-weight",collections.Counter((r[1],len(r[2])) for r in rows))
