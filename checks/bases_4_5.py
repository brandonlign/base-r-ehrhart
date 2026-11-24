"""Finite certificate and tail inequalities for bases 4 and 5 (Section 7).

For q = r - 1 in {3, 4} and each q-th root of unity zeta, write x_zeta = zeta(1 - y),
so that y = z phi(y). With y = q u, A(u) = S_q(zeta(1 - q u)), B(u) = (1 - (1 - q u)^q)/(q u)
and F = (q A / B)^n, Lagrange inversion gives, with C_j = [u^j] F/(1 - q u),

    [z^n] b_1     = C_n / q^(2n),
    [z^n] v_zeta  = q^2 C_(n-1) / (n q^(2n)),
    [z^n] b_zeta  = beta D_(n-1) / q^(2n)   (zeta != 1),

where beta = q/(1 - zeta), D_(-1) = 0 and D_j = C_j - zeta beta D_(j-1). The coefficients
of F come from F(0) = A(0)^n and the differential equation A B F' = n (A'B - A B') F.
Everything is computed exactly in Z[zeta] (Gaussian or Eisenstein integers). Upper
bounds for absolute values are rounded up to multiples of 10^-9.

Checked for 1 <= n <= N (N = 180 for q = 3, N = 900 for q = 4):
  principal branch: [z^n] b_1 >= 1/2, [z^n] v_1 > 0 for odd n and = 0 for even n;
  other branches:   |[z^n] b_zeta| < 1/5, and sum_n |[z^n] v_zeta| < 7/4 (q = 3),
                    12/5 (q = 4, zeta = i), 5/4 (q = 4, zeta = -1).
Then the tail inequalities 100 n R^-n < 1/100 and sum_{n > N} 100 R^-n < 1/100 for
R = 9/8 and 51/50. Expected output: one line per branch and a final line.
"""
from fractions import Fraction as Fr
from math import comb, isqrt

M = 10 ** 9


class Ring:
    """Z[w] with w^2 = u + t w: Gaussian (w = i: u = -1, t = 0) or Eisenstein (w = omega: u = -1, t = -1)."""

    def __init__(self, u, t):
        self.u, self.t = u, t

    def mul(self, x, y):
        (a, b), (c, d) = x, y
        return (a * c + self.u * b * d, a * d + b * c + self.t * b * d)

    def norm(self, x):
        """|a + b w|^2 = a^2 + t a b + b^2 for w = i or omega (Re w = t/2, |w| = 1); a^2 in the real case."""
        a, b = x
        return a * a + self.t * a * b + b * b if self.u == -1 else a * a


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def sc(x, c):
    return (x[0] * c, x[1] * c)


def up_abs(R, x, den):
    """An integer >= M |x| / den."""
    sq = R.norm(x) * M * M
    d2 = den * den
    v = isqrt(sq // d2)
    return v if v * v * d2 == sq else v + 1


def setup(mode):
    # zeta and beta = q/(1 - zeta) as elements of Z[w]
    if mode == "eis":
        R, zeta, beta = Ring(-1, -1), (0, 1), (2, 1)      # 3/(1 - omega) = 2 + omega
    elif mode == "i":
        R, zeta, beta = Ring(-1, 0), (0, 1), (2, 2)       # 4/(1 - i) = 2 + 2i
    elif mode == "neg":
        R, zeta, beta = Ring(0, 0), (-1, 0), (2, 0)       # 4/2
    else:
        R, zeta, beta = Ring(0, 0), (1, 0), None
    return R, zeta, beta


def series(q, mode, N):
    """Yield (n, C_n, C_(n-1), D_(n-1), q^(2n)) for n = 1..N, exactly in Z[w]."""
    R, zeta, beta = setup(mode)
    zero, one = (0, 0), (1, 0)
    # A(u) = sum_j zeta^j (1 - q u)^j and B(u) = sum_j (-q)^j C(q, j+1) u^j
    A = [zero] * (q + 1)
    power = one
    for j in range(q + 1):
        for i in range(j + 1):
            A[i] = add(A[i], sc(power, comb(j, i) * (-q) ** i))
        power = R.mul(power, zeta)
    B = [((-q) ** j * comb(q, j + 1), 0) for j in range(q)]

    def pmul(a, b):
        out = [zero] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] = add(out[i + j], R.mul(x, y))
        return out

    der = lambda a: [sc(x, i) for i, x in enumerate(a)][1:]
    D = pmul(A, B)
    U, V = pmul(der(A), B), pmul(A, der(B))
    J = [sub(U[i], V[i] if i < len(V) else zero) for i in range(len(U))]
    gam = R.mul(zeta, beta) if beta else None
    assert A[0][1] == 0 and D[0][1] == 0

    den = 1
    for n in range(1, N + 1):
        den *= q * q
        F = [(A[0][0] ** n, 0)]
        S = T = zero
        for p in range(n + 1):
            if p > 0:
                acc = zero
                for j in range(min(len(J) - 1, p - 1) + 1):
                    acc = add(acc, sc(R.mul(J[j], F[p - 1 - j]), n))
                for j in range(1, min(len(D) - 1, p - 1) + 1):
                    acc = sub(acc, sc(R.mul(D[j], F[p - j]), p - j))
                div = D[0][0] * p
                assert acc[0] % div == 0 and acc[1] % div == 0, "division"
                F.append((acc[0] // div, acc[1] // div))
            Sp = add(F[p], sc(S, q))                 # C_p
            if p == n:
                Cn = Sp
                break
            S = Sp                                   # C_(n-1) after the last pass
            if gam:
                T = sub(S, R.mul(gam, T))            # D_p
        yield n, Cn, S, T, den


def run(q, mode, N):
    R, zeta, beta = setup(mode)
    sum_v = 0
    max_b = 0
    min_b1 = None
    for n, Cn, S, T, den in series(q, mode, N):
        if mode == "real":
            ratio = Fr(Cn[0], den)                   # [z^n] b_1
            assert 2 * Cn[0] >= den, ("principal b", n)
            if n % 2:
                assert S[0] > 0, ("principal v", n)
            else:
                assert S[0] == 0, ("principal v not odd", n)
            min_b1 = ratio if min_b1 is None else min(min_b1, ratio)
        else:
            vb = up_abs(R, R.mul(beta, T), den)
            max_b = max(max_b, vb)
            sum_v += up_abs(R, sc(S, q * q), den * n)
    if mode == "real":
        print(f"q={q}, zeta=1, N={N}: [z^n] b_1 >= {float(min_b1):.4f} >= 1/2; [z^n] v_1 > 0 (n odd), = 0 (n even)")
        return
    limit = {"eis": Fr(7, 4), "i": Fr(12, 5), "neg": Fr(5, 4)}[mode]
    assert 5 * max_b < M, "b bound"
    assert Fr(sum_v, M) < limit, "v bound"
    name = {"eis": "omega", "i": "i", "neg": "-1"}[mode]
    print(f"q={q}, zeta={name}, N={N}: max |[z^n] b| <= {max_b / M:.4f} < 1/5; "
          f"sum |[z^n] v| <= {sum_v / M:.4f} < {limit}")


def branch_series(q, mode, nmax):
    """Exact coefficients [z^n] b_zeta and [z^n] v_zeta, n = 0..nmax, as pairs of Fractions in Q(w)."""
    R, zeta, beta = setup(mode)
    b = [(Fr(1), Fr(0)) if mode == "real" else (Fr(0), Fr(0))]
    v = [(Fr(0), Fr(0))]
    for n, Cn, S, T, den in series(q, mode, nmax):
        b.append((Fr(Cn[0], den), Fr(0)) if mode == "real"
                 else tuple(Fr(x, den) for x in R.mul(beta, T)))
        v.append((Fr(q * q * S[0], n * den), Fr(q * q * S[1], n * den)))
    return R, b, v


def generating_check(dmax=12):
    """Recover L_d(t) = sum_zeta b_zeta exp(t v_zeta) [z^d] from the branch series, d <= dmax,
    and compare with the h*-polynomial (ehrhart.py). Conjugate branches enter as 2 Re."""
    from math import factorial
    from ehrhart import coeffs
    for q, modes in ((3, (("real", 1), ("eis", 2))), (4, (("real", 1), ("i", 2), ("neg", 1)))):
        total = [[Fr(0)] * (dmax + 1) for _ in range(dmax + 1)]   # total[d][k] = [t^k z^d]
        for mode, weight in modes:
            R, b, v = branch_series(q, mode, dmax)
            smul = lambda x, y: [convolve(R, x, y, n) for n in range(dmax + 1)]
            term = b
            for k in range(dmax + 1):
                for d in range(dmax + 1):
                    a, c = term[d]
                    re = a + c * Fr(R.t, 2) if R.u == -1 else a   # Re(a + c w), Re w = t/2
                    total[d][k] += weight * re / factorial(k)
                term = smul(term, v)
        for d in range(1, dmax + 1):
            assert total[d][:d + 1] == coeffs(d, q + 1), (q, d)
            assert all(x == 0 for x in total[d][d + 1:])
    print(f"sum over branches of b exp(t v) reproduces L_d(t) for r = 4, 5 and d <= {dmax}")


def convolve(R, x, y, n):
    acc = (Fr(0), Fr(0))
    for i in range(n + 1):
        acc = add(acc, R.mul(x[i], y[n - i]))
    return acc


def tails():
    for Rr, N in ((Fr(9, 8), 180), (Fr(51, 50), 900)):
        c = N + 1
        assert 100 * c / Rr ** c < Fr(1, 100), "b tail"
        assert 100 / Rr ** c * Rr / (Rr - 1) < Fr(1, 100), "v tail"
        assert Fr(c + 1, c) < Rr, "monotonicity"
    print("tails: 100 n R^-n < 1/100 and sum_{n>N} 100 R^-n < 1/100 for (R, N) = (9/8, 180), (51/50, 900)")


if __name__ == "__main__":
    for q, mode, N in ((3, "real", 180), (3, "eis", 180), (4, "real", 900), (4, "i", 900), (4, "neg", 900)):
        run(q, mode, N)
    tails()
    generating_check()
    print("bases 4 and 5: certificate passed")
