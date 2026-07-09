from fractions import Fraction
from math import factorial

def exp_lower_pos(x, deg=12):
    # lower bound for e^x by truncating positive Taylor series
    x=Fraction(x)
    return sum(x**j/Fraction(factorial(j)) for j in range(deg+1))

# certify e^{-10/13}<58/125 by e^{10/13}>125/58
assert exp_lower_pos(Fraction(10,13),12) > Fraction(125,58)
# residual tails fit into rho margin
rho_margin = Fraction(153,100) - Fraction(152999,100000)
ip_tail = Fraction(153,53) * 3 * Fraction(11,14)**61
gs_tail = Fraction(153,53) * 3 * Fraction(7,13)**61
assert ip_tail < Fraction(1,250000)
assert ip_tail + gs_tail < rho_margin
print('exp and tail margin checks passed')
print('ip_tail_bound =', ip_tail)
print('gs_tail_bound =', gs_tail)
print('rho_margin =', rho_margin)
