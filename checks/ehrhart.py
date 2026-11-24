"""Ehrhart polynomials of the base-r simplices, computed three independent ways.

B_(r,d) = conv(e_1, ..., e_d, -q(1, r, ..., r^(d-1))) with q = r - 1, and L_d(t) is the
number of lattice points in t B_(r,d).

  ehr(d, r)       d! times the coefficients of L_d(t), from Solus's h*-polynomial
                  (a digit recurrence) and the binomial basis.
  conv(d, q, t)   L_d(t) from the convolution formula of Section 2 of the paper.
  brute(d, r, t)  L_d(t) from barycentric coordinates: a lattice point x lies in
                  t B_(r,d) iff mu = (t - sum x)/r^d >= 0 and x_i >= -mu q r^(i-1).
"""
from fractions import Fraction
from math import comb, factorial


def h_star(d, r):
    """h*-vector of B_(r,d) by a digit recurrence (Solus, Thm. 4.2)."""
    up = [[0] * (d + 1) for _ in range(r)]
    down = [[0] * (d + 1) for _ in range(r)]
    down[0][0] = 1
    for n in range(d):
        u = [[0] * (d + 1) for _ in range(r)]
        v = [[0] * (d + 1) for _ in range(r)]
        for k in range(n + 1):
            less = 0
            more = sum(up[j][k] + down[j][k] for j in range(r))
            for j in range(r):
                tot = up[j][k] + down[j][k]
                more -= tot
                u[j][k + 1] = less + up[j][k]
                v[j][k] = more + down[j][k]
                less += tot
        up, down = u, v
    return [sum(up[j][k] + down[j][k] for j in range(r)) for k in range(d + 1)]


def ehr(d, r):
    """[d! [t^k] L_d(t) for k = 0..d], using L_d(t) = sum_j h*_j binom(t + d - j, d)."""
    H = h_star(d, r)
    A = [0] * (d + 1)
    P = [1]                                   # d! binom(t + d, d) as a polynomial in t
    for i in range(1, d + 1):
        P = [a + b for a, b in zip([i * p for p in P] + [0], [0] + P)]
    for j in range(d + 1):
        for k, p in enumerate(P):
            A[k] += H[j] * p
        if j < d:                             # d! binom(t + d - j - 1, d) from d! binom(t + d - j, d)
            Q = [0] * d
            Q[-1] = P[-1]
            for k in range(d - 2, -1, -1):
                Q[k] = P[k + 1] - (d - j) * Q[k + 1]
            assert P[0] == (d - j) * Q[0]
            P = [a + b for a, b in zip([-j * x for x in Q] + [0], [0] + Q)]
    return A


def coeffs(d, r):
    """The coefficients of L_d(t) as Fractions."""
    return [Fraction(a, factorial(d)) for a in ehr(d, r)]


def conv(d, q, t):
    b = lambda n, k: comb(n, k) if n >= k >= 0 else 0
    return b(q * t + d, d) + sum(comb(d, j) * sum(b(t - u - 1, j - 1) * b(q * u + d - j, d - j) for u in range(t))
                                 for j in range(1, d + 1))


def brute(d, r, t):
    q, rd = r - 1, r ** d
    total = 0
    for s in range(-t * (rd - 1), t + 1):
        # mu = (t - s)/rd; x_i >= L_i = ceil(-mu q r^(i-1))
        L = sum(-((t - s) * q * r ** i // rd) for i in range(d))
        total += comb(s - L + d - 1, d - 1)
    return total


def value(c, t):
    return sum(x * t ** k for k, x in enumerate(c))


if __name__ == "__main__":
    n = 0
    for d in range(1, 5):
        for r in range(2, 8):
            c = coeffs(d, r)
            for t in range(0, d + 2):
                assert value(c, t) == brute(d, r, t) == conv(d, r - 1, t), (d, r, t)
                n += 1
    c = coeffs(5, 14)
    assert c == [1, Fraction(-268, 15), Fraction(3955, 12), Fraction(5833, 2), Fraction(10343, 12), Fraction(67228, 15)]
    for t in range(0, 3):
        assert value(c, t) == brute(5, 14, t) == conv(5, 13, t)
        n += 1
    for d in range(5, 13):
        for r in (2, 3, 6, 14):
            c = coeffs(d, r)
            for t in range(0, d + 2):
                assert value(c, t) == conv(d, r - 1, t), (d, r, t)
                n += 1
    print(f"h*-recurrence, convolution formula and direct lattice-point counts agree ({n} values)")
    print("L_5(t) for r = 14: 1 - 268/15 t + 3955/12 t^2 + 5833/2 t^3 + 10343/12 t^4 + 67228/15 t^5")
