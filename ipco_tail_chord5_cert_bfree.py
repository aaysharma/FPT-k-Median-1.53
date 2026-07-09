import argparse
import sympy as sp

X=sp.symbols('X')
theta=sp.Rational(11,14)

def D(m): return (1-theta*X/sp.Rational(m))**m

def E(m): return (theta*(1-X/sp.Rational(m)))**m

def T_lower(deg=5): return sum(((-theta*X)**j)/sp.factorial(j) for j in range(deg+1))

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

def check_range(start=6,end=60,verbose=False):
    D0=D(5); E0=E(5); T=T_lower(5)
    rows=[]
    for m in range(start,end+1):
        ce=bern_coeffs(E0-E(m)); mnE=min(ce)
        if mnE<0: raise RuntimeError(('E',m,mnE))
        P=sp.expand(T*(E0-E(m))+E(m)*D0-E0*D(m))
        coeffs=bern_coeffs(P); mn=min(coeffs)
        if mn<0: raise RuntimeError(('chord',m,mn))
        row={'m':m,'degree':sp.Poly(P,X).degree(),'min':str(mn)}
        rows.append(row)
        if verbose: print('m',m,'ok','deg',row['degree'],'min',row['min'],flush=True)
    print(f'IPCO chord certificate passed {start}..{end}')
    return rows

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--start',type=int,default=6)
    ap.add_argument('--end',type=int,default=60)
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args()
    check_range(args.start,args.end,args.verbose)
