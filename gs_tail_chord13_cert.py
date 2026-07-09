import argparse
import sympy as sp

X=sp.symbols('X')
lam=sp.Rational(10,13); gam=sp.Rational(7,13)

def D(k): return (1-lam*X/sp.Rational(k))**k

def E(k): return (gam*(1-X/sp.Rational(k)))**k

def T_lower(deg=5): return sum(((-lam*X)**j)/sp.factorial(j) for j in range(deg+1))

def bern_coeffs(poly, deg=None):
    poly=sp.expand(poly)
    if deg is None: deg=sp.Poly(poly,X).degree()
    l0,l1=sp.symbols('l0 l1')
    expr=sp.expand(poly.subs(X,l1))
    P=sp.Poly(expr,l0,l1,domain=sp.QQ)
    S=l0+l1; H=0
    for mon,coef in P.terms():
        total=sum(mon)
        H += coef*l0**mon[0]*l1**mon[1]*S**(deg-total)
    HP=sp.Poly(sp.expand(H),l0,l1,domain=sp.QQ)
    out=[]
    for i in range(deg+1):
        j=deg-i
        out.append(sp.simplify(HP.coeff_monomial(l0**i*l1**j)/sp.binomial(deg,j)))
    return out

def check_range(start=14,end=60,verbose=False):
    D0=D(13); E0=E(13); T=T_lower(5)
    rows=[]
    for k in range(start,end+1):
        ce=bern_coeffs(E0-E(k)); mnE=min(ce)
        if mnE<0: raise RuntimeError(('E',k,mnE))
        P=sp.expand(T*(E0-E(k))+E(k)*D0-E0*D(k))
        coeffs=bern_coeffs(P); mn=min(coeffs)
        if mn<0: raise RuntimeError(('chord',k,mn))
        row={'k':k,'degree':sp.Poly(P,X).degree(),'min':str(mn)}
        rows.append(row)
        if verbose: print('k',k,'ok','deg',row['degree'],'min',row['min'],flush=True)
    print(f'GS chord certificate passed {start}..{end}')
    return rows

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--start',type=int,default=14)
    ap.add_argument('--end',type=int,default=60)
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args()
    check_range(args.start,args.end,args.verbose)
