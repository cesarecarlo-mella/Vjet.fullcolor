# integrability of the eps^0 symbols with the alphabet of letters.dict.m (even letters only)
import sys,os,random,collections
_argv=sys.argv[:]; sys.argv=['x','NONE']
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'classical.py')).read().split("TESTS=")[0]); sys.argv=_argv
# letters as polynomials in x,y: dict {(i,j):c}
POL={1:{(1,0):1},2:{(0,1):1},3:{(0,0):1,(1,0):-1,(0,1):-1},4:{(0,0):1,(1,0):-1},5:{(0,0):1,(0,1):-1},6:{(1,0):1,(0,1):1},
 11:{(1,0):1,(0,1):-1,(0,2):1},12:{(1,0):-1,(2,0):1,(0,1):1},13:{(2,0):1,(1,0):-2,(0,0):1,(0,1):-1},
 14:{(2,0):1,(1,1):2,(0,2):1,(0,1):-1},15:{(2,0):1,(1,1):2,(0,2):1,(1,0):-1},16:{(1,0):-1,(0,2):1,(0,1):-2,(0,0):1}}
def ev(p,X,Y): return sum(c*pow(X,i,P)*pow(Y,j,P) for (i,j),c in p.items())%P
def dx(p): return {(i-1,j):c*i for (i,j),c in p.items() if i}
def dy(p): return {(i,j-1):c*j for (i,j),c in p.items() if j}
pts=[(random.randrange(P),random.randrange(P)) for _ in range(40)]
def omega(a,b,X,Y):
    A,B=POL[a],POL[b]
    return (ev(dx(A),X,Y)*ev(dy(B),X,Y)-ev(dy(A),X,Y)*ev(dx(B),X,Y))%P*inv(ev(A,X,Y)*ev(B,X,Y)%P)%P
OM={(a,b,k):omega(a,b,*pts[k]) for a in POL for b in POL for k in range(len(pts))}
rnd=random.Random(5)
for amp in sys.argv[1:]:
    tot=bad=0
    for w in (4,5,6):
        g=collections.defaultdict(collections.Counter)
        for c,wd,d,co,p in [(r[1],r[2],r[3],r[4],r[5]) for r in by[amp]]:
            if c=='' and len(wd)==w: g[(p,d)][wd]+=evc(co)
        if not g: continue
        env={'t':rnd.randrange(P),'u':rnd.randrange(P)}
        F=collections.Counter()
        for (p,d),vec in g.items():
            r=ev_(p,d) if False else eval(p.replace('^','**'),env)%P*inv(eval(d.replace('^','**'),env)%P)%P
            for x_,c in vec.items(): F[x_]+=r*c
        for k in range(w-1):
            C=collections.defaultdict(collections.Counter)
            skip=set()
            for wd,c in F.items():
                key=(wd[:k],wd[k+2:])
                if {wd[k],wd[k+1]}&{7,8}: skip.add(key); continue
                C[key][(wd[k],wd[k+1])]+=c
            for key,pr in C.items():
                if key in skip: continue
                tot+=1
                if any(sum(c*OM[(a,b,j)] for (a,b),c in pr.items())%P for j in range(len(pts))): bad+=1
    print(amp,"integrability conditions tested:",tot,"violated:",bad)
