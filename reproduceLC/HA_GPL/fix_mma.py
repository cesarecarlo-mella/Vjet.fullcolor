import re,sys
def fix(s,name):
    s=re.sub(r'^\(\* HA(?:G|ser)\[[^:]*:', '(* %s:'%name, s, count=1)
    s,k=re.subn(r'^HA(?:G|ser)\[.*? = \(?', name+' = (', s, count=1, flags=re.M)
    assert k==1, name
    s=s.rstrip()
    if not s.endswith(');'):
        s=(s[:-1] if s.endswith(';') else s)+');'
    s+='\n'
    body=re.sub(r'\(\*.*?\*\)','',s,flags=re.S)
    assert body.count('(')==body.count(')') and body.count('[')==body.count(']'), name
    return s
if __name__=='__main__':
    for f in sys.argv[1:]:
        m=re.search(r'(C\d+F\d+Z\d+)_HA(\d)_(G|series15)\.m$',f)
        L,i,kind=m.groups()
        name=('HAG["%s", %s]' if kind=='G' else 'HAser["%s", %s]')%(L,i)
        s=open(f).read()
        if not s.strip(): print('EMPTY',f); continue
        s=fix(s,name)
        open(f,'w').write(s)
        print('ok',f, [l for l in s.split('\n') if l.startswith('HA')][0][:40], repr(s[-12:]))
