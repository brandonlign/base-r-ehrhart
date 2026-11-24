"""Checks for the proof that bases 2 and 3 are Ehrhart positive (Section 6).

1. q = 1: sum_d L_d(t) z^d = (1-z)^(-1) ((1+z)/(1-z))^t, compared with the
   h*-polynomial for d <= 14.
2. q = 2: sum_d L_d(t) z^d = b_+ e^(t l_+) + b_- e^(t l_-) with
   l_pm' = 2/(1-z^2) (1 pm (4-3z^2)^(-1/2)), l_pm(0) = 0 and
   b_pm = (1 pm (2-z)(4-3z^2)^(-1/2)) / (2(1-z)), compared for d <= 14.
3. The coefficient facts used in the proof, through degree 200: l_pm' have
   nonnegative coefficients, [z] l_+ > 0, b_+ has positive coefficients, and
   b_- has nonnegative coefficients.
4. Positivity of every coefficient of L_d(t) for r = 2, 3 and d <= 60.
All arithmetic is exact. Expected output: four "ok" lines.
"""
from fractions import Fraction as F
from math import comb, factorial
from ehrhart import coeffs

N = 200


def mul(a, b, n):
    out = [F(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                out[i + j] += x * y
    return out


def integrate(a, n):
    return [F(0)] + [a[k] / (k + 1) for k in range(n)]


def exp_t(L, n):
    """[e^(t L)] as a list over k of series in z: e^(tL) = sum_k t^k L^k / k!."""
    out, P = [], [F(1)] + [F(0)] * n
    for k in range(n + 1):
        out.append([x / factorial(k) for x in P])
        P = mul(P, L, n)
    return out


def check_gf(b_pm, l_pm, r, n):
    terms = [exp_t(l, n) for l in l_pm]
    for d in range(1, n + 1):
        got = [sum(mul(b, T[k], d)[d] for b, T in zip(b_pm, terms)) for k in range(d + 1)]
        assert got == coeffs(d, r), (r, d)


inv_sqrt = [F(0)] * (N + 1)                      # (4 - 3z^2)^(-1/2)
for i in range(N // 2 + 1):
    inv_sqrt[2 * i] = F(comb(2 * i, i)) * F(3, 16) ** i / 2
geom = [F(1)] * (N + 1)                          # 1/(1 - z)
geom2 = [F(1) if k % 2 == 0 else F(0) for k in range(N + 1)]   # 1/(1 - z^2)

# 1. q = 1
n = 14
log_ratio = [F(0)] + [F(2, k) if k % 2 else F(0) for k in range(1, n + 1)]   # log((1+z)/(1-z))
check_gf([geom[:n + 1]], [log_ratio], 2, n)
print(f"q=1: generating function matches L_d(t) for d <= {n}: ok")

# 2. q = 2
lp_d = [2 * x for x in mul(geom2, [(F(1) if k == 0 else F(0)) + inv_sqrt[k] for k in range(N + 1)], N)]
lm_d = [2 * x for x in mul(geom2, [(F(1) if k == 0 else F(0)) - inv_sqrt[k] for k in range(N + 1)], N)]
two_minus_z = [F(2), F(-1)] + [F(0)] * (N - 1)
s = mul(two_minus_z, inv_sqrt, N)
b_p = [x / 2 for x in mul(geom, [(F(1) if k == 0 else F(0)) + s[k] for k in range(N + 1)], N)]
b_m = [x / 2 for x in mul(geom, [(F(1) if k == 0 else F(0)) - s[k] for k in range(N + 1)], N)]
l_p, l_m = integrate(lp_d, N), integrate(lm_d, N)
check_gf([b_p[:n + 1], b_m[:n + 1]], [l_p[:n + 1], l_m[:n + 1]], 3, n)
print(f"q=2: generating function matches L_d(t) for d <= {n}: ok")

# 3. coefficient facts
assert all(x >= 0 for x in lp_d) and all(x >= 0 for x in lm_d) and l_p[1] > 0
assert all(x > 0 for x in b_p) and all(x >= 0 for x in b_m)
print(f"q=2: l_+', l_-' >= 0, [z]l_+ > 0, b_+ > 0, b_- >= 0 coefficientwise through degree {N}: ok")

# 4. direct positivity
for r in (2, 3):
    for d in range(1, 61):
        assert all(x > 0 for x in coeffs(d, r)), (r, d)
print("r = 2, 3: all coefficients of L_d(t) positive for d <= 60: ok")
