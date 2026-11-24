"""Checks for Theorems 1.2 and 1.3, Proposition 4.1, Section 5 and Table 1.

1. Theorem 1.2 (large base): for d <= 16, each coefficient [t^k] L_d(t) is a polynomial
   in q of degree <= d (recovered by exact interpolation at q = 1..d+1 and confirmed
   at q = d+2, d+3); its leading term has the sign stated in the theorem.
2. Proposition 4.1 and Theorem 1.3 (codimension four): the closed form 4320 q^4 A_(d,4) = 6 (1+q)^d P + (q-1) T
   agrees with the h*-polynomial for 5 <= d <= 20 and 1 <= q <= 35, and the finite
   rows 6 <= d <= 11, 1 <= q < 100 d have the claimed sign.
3. The low-dimensional formulas of Section 4 (as polynomials in r, checked at d+1 bases).
4. The first-failure table: for each base r in 6..14 and every d below the listed
   dimension all coefficients are positive, and at the listed dimension the linear
   coefficient is negative.
5. The fourth coefficient of exp(((v/2) coth(v/2) - 1)/c) is (5 - 2c)/(1440 c^2).
Expected output: five "ok" lines.
"""
from fractions import Fraction as F
from math import factorial
from ehrhart import ehr


def interpolate(xs, ys):
    """Coefficients (low to high) of the polynomial through the points (xs, ys)."""
    n = len(xs)
    out = [F(0)] * n
    for i in range(n):
        basis, den = [F(1)], F(1)
        for j in range(n):
            if j != i:
                basis = [a - xs[j] * b for a, b in zip([F(0)] + basis, basis + [F(0)])]
                den *= xs[i] - xs[j]
        for k in range(n):
            out[k] += ys[i] * basis[k] / den
    return out


def peval(c, x):
    return sum(a * x ** k for k, a in enumerate(c))


# 1. large base
for d in range(1, 17):
    qs = list(range(1, d + 2))
    rows = {q: ehr(d, q + 1) for q in qs + [d + 2, d + 3]}
    for k in range(d + 1):
        poly = interpolate(qs, [rows[q][k] for q in qs])
        assert all(peval(poly, q) == rows[q][k] for q in (d + 2, d + 3))
        lead = next(a for a in reversed(poly) if a != 0)
        gap = d - k
        negative = k > 0 and gap >= 4 and gap % 4 in (0, 1)
        assert (lead < 0) == negative, (d, k)
print("Theorem 1.2: eventual signs of all coefficients, d <= 16: ok")


# 2. codimension four
def closed_form(d, q):
    P = -1 + 4 * d * q + (d + 5 * (d * (d - 1) // 2)) * q * q - d * q ** 3
    T = (d * (d - 1) * (d - 2) * q ** 6 + d * (d - 1) * (d - 11) * q ** 5 + (-19 * d ** 3 + 18 * d ** 2 + 19 * d) * q ** 4
         + (26 * d ** 3 + 18 * d ** 2 + 4 * d - 6) * q ** 3 + (36 * d ** 2 + 12 * d - 6) * q ** 2 + (18 * d - 6) * q - 6)
    return 6 * (1 + q) ** d * P + (q - 1) * T


count = 0
for d in range(5, 21):
    for q in range(1, 36):
        # A_(d,4) = (d-4)! [t^(d-4)] L_d and ehr gives d! [t^k] L_d
        assert ehr(d, q + 1)[d - 4] * factorial(d - 4) * 4320 * q ** 4 == factorial(d) * closed_form(d, q), (d, q)
        count += 1
for d in range(6, 12):
    for q in range(1, 100 * d):
        assert closed_form(d, q) * (1 if 2 * q <= 5 * d - 3 else -1) > 0, (d, q)
# terms of T as {(power of q, power of d): coefficient}; |T| <= (sum of |coefficients|) d^3 q^6 for d, q >= 1
T_terms = {(6, 3): 1, (6, 2): -3, (6, 1): 2, (5, 3): 1, (5, 2): -12, (5, 1): 11, (4, 3): -19, (4, 2): 18, (4, 1): 19,
           (3, 3): 26, (3, 2): 18, (3, 1): 4, (3, 0): -6, (2, 2): 36, (2, 1): 12, (2, 0): -6, (1, 1): 18, (1, 0): -6,
           (0, 0): -6}
for d in range(5, 21):
    for q in range(1, 36):
        T = sum(c * q ** i * d ** j for (i, j), c in T_terms.items())
        assert 6 * (1 + q) ** d * (-1 + 4 * d * q + (d + 5 * (d * (d - 1) // 2)) * q * q - d * q ** 3) + (q - 1) * T \
            == closed_form(d, q)
abs_sum = sum(abs(c) for c in T_terms.values())
assert abs_sum == 224
assert F(abs_sum, 6) * F(144, 4096) < 3 and F(abs_sum, 3) < 75
print(f"Proposition 4.1 and Theorem 1.3: closed form matches {count} values; rows 6 <= d <= 11, q < 100d have the stated sign; "
      f"|(q-1)T| <= {abs_sum} d^3 q^7: ok")


# 3. low dimensions: d! [t^k] L_d(t)
def low(d, r):
    if d == 1:
        return [1, r]
    if d == 2:
        return [2, r + 2, r ** 2]
    if d == 3:
        return [6, F(3 * r * r + r + 18, 2), F(3 * (r * r + r + 2), 2), r ** 3]
    if d == 4:
        return [24, 2 * (2 * r * r + r + 22), 4 * r ** 3 + 5 * r * r + 2 * r + 24, 2 * (r ** 3 + r * r + r + 2), r ** 4]
    return [120, F(-5 * r ** 4 + 60 * r ** 3 + 65 * r * r + 24 * r + 1500, 6), F(25 * (r ** 3 + 2 * r * r + r + 14), 2),
            F(5 * (5 * r ** 4 + 6 * r ** 3 + 7 * r * r + 3 * r + 30), 3), F(5 * (r ** 4 + r ** 3 + r * r + r + 2), 2), r ** 5]


for d in range(1, 6):
    for r in range(2, d + 3):
        assert low(d, r) == ehr(d, r), (d, r)
g = lambda r: -5 * r ** 4 + 60 * r ** 3 + 65 * r * r + 24 * r + 1500
assert g(13) > 0 > g(14)
print("low-dimensional formulas (d <= 5) and the sign change of [t] L_5 between r = 13 and 14: ok")

# 4. first failures
table = {6: 43, 7: 19, 8: 15, 9: 7, 10: 7, 11: 7, 12: 7, 13: 7, 14: 5}
for r, first in table.items():
    for d in range(1, first + 1):
        c = ehr(d, r)
        if d < first:
            assert min(c) > 0, (r, d)
        else:
            assert c[1] < 0 and all(x > 0 for k, x in enumerate(c) if k != 1), (r, d)
print("first-failure table for r = 6..14: ok")

# 5. limiting coefficient: (v/2) coth(v/2) - 1 = v^2/12 - v^4/720 + O(v^6)
for c in (F(1), F(5, 2), F(7, 3), F(10)):
    g2, g4 = F(1, 12) / c, F(-1, 720) / c
    assert g4 + g2 ** 2 / 2 == (5 - 2 * c) / (1440 * c * c)
print("fourth coefficient of the limiting exponential is (5 - 2c)/(1440 c^2): ok")
