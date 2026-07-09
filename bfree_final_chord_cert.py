from fractions import Fraction
from math import comb, factorial
from functools import lru_cache
import json, time, os

rho=Fraction(152999,100000)  # leaves margin for multiplicity residuals
eps=Fraction(7,50)
a_layer=Fraction(150,107)
theta=Fraction(11,14)
lam=Fraction(10,13)
gamma=Fraction(7,13)
alpha_ip=Fraction(107,57); beta_ip=Fraction(157,57)
alpha_ug=Fraction(2,1); beta_ug=Fraction(3,1)
alpha_g=Fraction(107,50); beta_g=Fraction(3,1)
# certified upper endpoint e^{-10/13}<58/125; chord upper U(T)=1-(1-u)T.
exp1_upper=Fraction(58,125)
chord_slope=Fraction(1,1)-exp1_upper  # 67/125
c_supported = Fraction(1,1) - Fraction(1,1)/rho

def clean(p): return {k:v for k,v in p.items() if v}
def const(c,n):
    c=Fraction(c); return {} if c==0 else {tuple([0]*n):c}
def var(j,n):
    e=[0]*n; e[j]=1; return {tuple(e):Fraction(1)}
def add(p,q,scale=Fraction(1)):
    out=p.copy()
    for e,c in q.items():
        v=out.get(e,Fraction(0))+scale*c
        if v: out[e]=v
        elif e in out: del out[e]
    return out
def sub(p,q): return add(p,q,Fraction(-1))
def smul(c,p):
    c=Fraction(c)
    if c==0 or not p: return {}
    return {e:c*v for e,v in p.items() if c*v}
def mul(p,q,n):
    if not p or not q: return {}
    out={}
    for e1,c1 in p.items():
        for e2,c2 in q.items():
            e=tuple(e1[i]+e2[i] for i in range(n))
            out[e]=out.get(e,Fraction(0))+c1*c2
    return clean(out)
def pow_poly(p,k,n):
    res=const(1,n); base=p
    while k:
        if k&1: res=mul(res,base,n)
        k//=2
        if k: base=mul(base,base,n)
    return res
def one_minus(c,j,n): return add(const(1,n), smul(-Fraction(c), var(j,n)))
def degree(p): return max((sum(e) for e in p), default=0)

def group_factors(kind,param,idx,n):
    if kind=='UG':
        D=one_minus(theta,idx,n); E=smul(theta, one_minus(1,idx,n)); return D,E,alpha_ug,beta_ug
    if kind=='IP':
        m=int(param); D=pow_poly(one_minus(theta/Fraction(m),idx,n),m,n); E=pow_poly(smul(theta, one_minus(Fraction(1,m),idx,n)),m,n); return D,E,alpha_ip,beta_ip
    if kind=='GS':
        if param=='inf':
            # linear chord upper for e^{-lam x}. E=0.
            D=add(const(1,n), smul(-chord_slope, var(idx,n)))
            return D, {}, alpha_g, beta_g
        k=int(param); D=pow_poly(one_minus(lam/Fraction(k),idx,n),k,n); E=pow_poly(smul(gamma, one_minus(Fraction(1,k),idx,n)),k,n); return D,E,alpha_g,beta_g
    raise ValueError((kind,param))

def cases():
    K=[None]+list(range(1,14))+['inf']
    M=[None]+list(range(1,6))
    out=[]; idx=0
    for k in K:
        for ug in [False,True]:
            for m in M:
                blocks=[]; name=[]
                if k is not None: blocks.append(('GS',k,Fraction(1))); name.append(f'GS{k}@1')
                if ug: blocks.append(('UG',None,Fraction(1))); name.append('UG@1')
                if m is not None: blocks.append(('IP',m,Fraction(1))); name.append(f'IP{m}@1')
                blocks.append(('GS','inf',a_layer)); name.append('GSinf@a')
                out.append((idx,'+'.join(name),blocks)); idx+=1
    return out

def slack_poly(blocks):
    n=len(blocks)
    xs=[var(i,n) for i in range(n)]
    S={}
    for x in xs: S=add(S,x)
    e=sub(const(1,n),S)
    Ds=[]; Es=[]; al=[]; be=[]; qs=[]
    for i,(kind,param,q) in enumerate(blocks):
        D,E,A,B=group_factors(kind,param,i,n); Ds.append(D); Es.append(E); al.append(A); be.append(B); qs.append(q)
    X={}
    for i,q in enumerate(qs): X=add(X, smul(q,xs[i]))
    N={}; pref=const(1,n)
    for i in range(n):
        N=add(N, mul(pref, smul(qs[i], sub(const(1,n),Ds[i])), n))
        pref=mul(pref,Ds[i],n)
    for i in range(n):
        pre=const(1,n); post=const(1,n)
        for h in range(i): pre=mul(pre,Es[h],n)
        for h in range(i+1,n): post=mul(post,Ds[h],n)
        term=mul(pre,sub(Ds[i],Es[i]),n); term=mul(term,post,n)
        N=add(N,smul(al[i]*qs[i],term))
    fall=const(1,n)
    for E in Es: fall=mul(fall,E,n)
    N=add(N,smul(be[0],fall))
    c=sub(const(1,n), smul(rho,e))
    F=sub(smul(rho,X), mul(c,N,n))
    return F,n

@lru_cache(None)
def multi_leq(total,length):
    if length==0: return ((),)
    if length==1: return tuple((i,) for i in range(total+1))
    out=[]
    for i in range(total+1):
        for rest in multi_leq(total-i,length-1): out.append((i,)+rest)
    return tuple(out)
@lru_cache(None)
def binom(n,k): return 0 if k<0 or k>n else comb(n,k)
@lru_cache(None)
def multinomial_for_monom(N, alpha):
    s=sum(alpha)
    if s>N: return Fraction(0)
    out=Fraction(factorial(N),factorial(N-s))
    for a in alpha: out/=factorial(a)
    return out
@lru_cache(None)
def monom_coeff(alpha,beta,N):
    for a,b in zip(alpha,beta):
        if a>b: return Fraction(0)
    num=Fraction(1)
    for a,b in zip(alpha,beta): num*=binom(b,a)
    den=multinomial_for_monom(N,alpha)
    return Fraction(0) if den==0 else num/den

def bernstein_min_standard(poly,n,deg=None):
    if deg is None: deg=degree(poly)
    items=list(poly.items()); minc=None; neg=0; idx_min=None
    for beta in multi_leq(deg,n):
        s=Fraction(0)
        for alpha,cf in items: s += cf*monom_coeff(alpha,beta,deg)
        if minc is None or s<minc: minc=s; idx_min=beta
        if s<0: neg+=1
    return minc,neg,idx_min,deg

def affine_substitute(poly,verts):
    d=len(verts[0]); n=d
    L=[]
    for j in range(d):
        p=const(verts[0][j],n)
        for k in range(1,d+1):
            diff=verts[k][j]-verts[0][j]
            if diff: p=add(p, smul(diff,var(k-1,n)))
        L.append(p)
    maxdeg=degree(poly); powers=[]
    for j in range(d):
        pw=[const(1,n)]
        for r in range(1,maxdeg+1): pw.append(mul(pw[-1],L[j],n))
        powers.append(pw)
    out={}
    for exp,cf in poly.items():
        term=const(cf,n)
        for j,e in enumerate(exp):
            if e: term=mul(term,powers[j][e],n)
        out=add(out,term)
    return out

def split_simplex(verts):
    d=len(verts[0]); best=None
    for i in range(len(verts)):
        for j in range(i+1,len(verts)):
            dist=sum(float(verts[i][k]-verts[j][k])**2 for k in range(d))
            if best is None or dist>best[0]: best=(dist,i,j)
    _,i,j=best; mid=tuple((verts[i][k]+verts[j][k])/2 for k in range(d))
    v1=list(verts); v2=list(verts); v1[i]=mid; v2[j]=mid
    return [v1,v2]

def supported_region_triangulation(n):
    c=c_supported
    upper=[]; lower=[]
    for i in range(n):
        u=[Fraction(0)]*n; u[i]=Fraction(1); upper.append(tuple(u))
        l=[Fraction(0)]*n; l[i]=c; lower.append(tuple(l))
    sims=[]
    for k in range(n):
        verts=[]
        for i in range(k+1): verts.append(upper[i])
        for i in range(k,n): verts.append(lower[i])
        sims.append(tuple(verts))
    return sims

def check_poly(poly,n,max_depth=10,deg_raise=0):
    base_deg=degree(poly)+deg_raise
    stack=[(list(vs),0) for vs in supported_region_triangulation(n)]
    nodes=0; fails=[]; min_leaf=None; min_seen=None
    while stack:
        verts,dep=stack.pop(); nodes+=1
        loc=affine_substitute(poly,verts)
        mn,neg,idx,deg=bernstein_min_standard(loc,n,base_deg)
        if min_seen is None or mn<min_seen: min_seen=mn
        if neg:
            if dep<max_depth: stack.extend((sv,dep+1) for sv in split_simplex(verts))
            else:
                fails.append({'min':str(mn),'neg':neg,'depth':dep,'idx':idx})
                if len(fails)>3: break
        else:
            if min_leaf is None or mn<min_leaf: min_leaf=mn
    return len(fails)==0, (min_leaf if min_leaf is not None else min_seen), nodes, fails[:2], base_deg

def check_case(case,max_depth=10,deg_raise=0):
    idx,name,blocks=case; poly,n=slack_poly(blocks)
    ok,mn,nodes,fails,deg=check_poly(poly,n,max_depth,deg_raise)
    return {'idx':idx,'name':name,'blocks':str(blocks),'ok':ok,'min':str(mn),'nodes':nodes,'fails':fails,'nfails':len(fails),'degree':deg,'degree_raise':deg_raise,'M_upper':'global chord with exp(−10/13)<58/125','region':'sum_x >= '+str(c_supported)}

def load(path):
    if os.path.exists(path):
        d=json.load(open(path)); return {int(r['idx']):r for r in d.get('rows',[])}
    return {}
def save(path,rows):
    arr=[rows[i] for i in sorted(rows)]; bad=[r for r in arr if not r['ok']]
    data={'rho':str(rho),'M_upper':'M(T)<=1-(67/125)T using exp(-10/13)<58/125','automatic_region':'e >= 1/rho = '+str(Fraction(1,1)/rho),'checked_region':'sum_x >= '+str(c_supported),'total_cases':len(cases()),'rows':arr,'count':len(arr),'bad':bad,'bad_count':len(bad)}
    tmp=path+'.tmp'; json.dump(data,open(tmp,'w'),indent=2); os.replace(tmp,path)

def main():
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--start',type=int,default=0); ap.add_argument('--end',type=int,default=180); ap.add_argument('--out',default='/mnt/data/bfree_final_chord_cert_results.json'); ap.add_argument('--max-depth',type=int,default=10); ap.add_argument('--degree-raise',type=int,default=0); ap.add_argument('--limit',type=int)
    args=ap.parse_args(); allc=cases(); rows=load(args.out); todo=[c for c in allc[args.start:args.end] if c[0] not in rows]
    if args.limit: todo=todo[:args.limit]
    st=time.time()
    for case in todo:
        rec=check_case(case,args.max_depth,args.degree_raise); rows[rec['idx']]=rec; save(args.out,rows)
        print(('OK' if rec['ok'] else 'BAD'),rec['idx'],rec['name'],'nodes',rec['nodes'],'min',rec['min'],'done',len(rows),flush=True)
        if not rec['ok']: raise SystemExit('failure')
    save(args.out,rows); print('DONE',len(todo),'saved',len(rows),'elapsed',time.time()-st)
if __name__=='__main__': main()
