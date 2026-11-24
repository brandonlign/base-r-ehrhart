"""Checks for the analytic-continuation step of Section 7 (bases 4 and 5).

1. With x = (1 - w)/(1 + w), the map z = (1 - x^q)/S_q(x) equals
   w(3 + w^2)/(2(1 + w^2)) for q = 3 and 8w(1 + w^2)/(5 + 10w^2 + w^4) for q = 4,
   and the numerators of dz/dw are proportional to w^4 + 3 and u^3 - 7u^2 - 5u - 5
   (u = w^2). Both sides are rational functions of bounded degree, so agreement at
   40 rational points proves each identity.
2. The critical values have modulus 3^(3/4)/2 > 9/8 (q = 3) and > 51/50 (q = 4),
   by rational interval bounds on the roots of u^3 - 7u^2 - 5u - 5.
Expected output: two "ok" lines.
"""
from fractions import Fraction as Q

pts = [Q(k, 7) for k in range(1, 41)]

def z_of_x(x, q):
    return (1 - x ** q) / sum(x ** j for j in range(q + 1))

for w in pts:
    x = (1 - w) / (1 + w)
    assert z_of_x(x, 3) == w * (3 + w * w) / (2 * (1 + w * w))
    assert z_of_x(x, 4) == 8 * w * (1 + w * w) / (5 + 10 * w ** 2 + w ** 4)
    # derivative numerators, from (N/D)' = (N'D - ND')/D^2
    N3, D3, dN3, dD3 = w * (3 + w * w), 2 * (1 + w * w), 3 + 3 * w * w, 4 * w
    assert dN3 * D3 - N3 * dD3 == 2 * (w ** 4 + 3)
    N4, D4, dN4, dD4 = 8 * w * (1 + w * w), 5 + 10 * w ** 2 + w ** 4, 8 + 24 * w * w, 20 * w + 4 * w ** 3
    u = w * w
    assert dN4 * D4 - N4 * dD4 == -8 * (u ** 3 - 7 * u ** 2 - 5 * u - 5)
print("rational maps in w and their critical-point equations: ok")

# q = 3: w^4 = -3, so w^2 = +-i sqrt3, |3 + w^2|^2 = 12, |1 + w^2|^2 = 4,
# |z|^4 = |w|^4 * 12^2 / (2^4 * 4^2) = 3 * 144 / 256 = 27/16.
assert Q(3 * 144, 256) == Q(27, 16) > Q(9, 8) ** 4

# q = 4: one real root U of g, the other two are a +- ib with 2a + U = 7 and (a^2 + b^2) U = 5.
g = lambda u: u ** 3 - 7 * u ** 2 - 5 * u - 5
lo, hi = Q("7.73045"), Q("7.73047")
assert g(lo) < 0 < g(hi)
# real critical value: |z|^2 = 64 U (1 + U)^2 / (5 + 10 U + U^2)^2, increasing numerator in U
assert 64 * lo * (1 + lo) ** 2 > Q("1.36") ** 2 * (5 + 10 * hi + hi ** 2) ** 2
# nonreal critical values; a = (7 - U)/2 lies in (-0.37, -0.36), m = |u|^2 = 5/U
a_lo, a_hi = (7 - hi) / 2, (7 - lo) / 2
m_lo, m_hi = 5 / hi, 5 / lo
assert a_lo < a_hi < 0
s_lo = Q("0.80423")                                   # lower bound for sqrt(m)
assert s_lo ** 2 < m_lo
num_lo = 64 * s_lo * (1 + 2 * a_lo + m_lo)           # 64 |u| |1 + u|^2, increasing in a and m
re_lo = 5 + 10 * a_lo + 2 * a_hi ** 2 - m_hi          # Re(5 + 10u + u^2) = 5 + 10a + 2a^2 - m
re_hi = 5 + 10 * a_hi + 2 * a_lo ** 2 - m_lo          # (a < 0, so a^2 is largest at a_lo)
assert 0 < re_lo <= re_hi
den_hi = re_hi ** 2 + (10 + 2 * a_hi) ** 2 * (m_hi - a_hi ** 2)   # |Im|^2 = (10 + 2a)^2 (m - a^2)
assert num_lo > Q("47.15") and den_hi < Q("45.06") and Q("47.15") / Q("45.06") > Q(51, 50) ** 2
print("critical values: modulus^4 = 27/16 > (9/8)^4 for q = 3; modulus > 1.36 (real) and > 51/50 (nonreal) for q = 4: ok")
