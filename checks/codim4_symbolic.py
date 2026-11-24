"""Symbolic derivation of the codimension-four formula (needs sympy).

Starting from the residue formula (Lemma 2.2),
    A_{d,4} = [v^4] sum_{zeta^q = 1} (v/q) a_zeta(v) (v (a_zeta(v) + b(v)))^d,
it expands the principal term with q and d as symbols, writes each other term
as a polynomial in alpha = 1/(1 - zeta), and replaces the powers of alpha by
their sums over the nontrivial q-th roots of unity, read off from the
logarithmic derivative of ((1+s)^q - 1)/s. The result is compared with
4320 q^4 A_{d,4} = 6 (1+q)^d P + (q-1) T of Proposition 4.1, identically in q and d.
"""
import sympy as sp

v, s, q, d, al, Lam = sp.symbols("v s q d alpha Lambda")
ORDER = 6


def ser(expr, n=ORDER):
    return sp.series(expr, v, 0, n).removeO()


def power_d(w):
    """(1 + w)^d as a series in v, for a series w with w(0) = 0."""
    L = ser(sp.log(1 + w))
    return sp.expand(ser(sp.exp(d * L)))


def coeff4(expr):
    return sp.simplify(sp.expand(expr).coeff(v, 4))


b = 1 / (sp.exp(v) - 1)

# Principal root: v a_1 = v/(1 - e^{-v/q}) and v(a_1 + b) = (q+1)(1 + w).
a1 = 1 / (1 - sp.exp(-v / q))
U = ser(v * a1 / q)
w1 = sp.expand(ser(v * (a1 + b)) / (q + 1) - 1)
principal = coeff4(sp.expand(U * power_d(w1)))  # times (1+q)^d

# Other roots: zeta = 1 - 1/alpha, so a_zeta = 1/(1 - zeta e^{-v/q}).
zeta = 1 - 1 / al
az = 1 / (1 - zeta * sp.exp(-v / q))
Uz = ser(v * az / q)
wz = sp.expand(ser(v * (az + b)) - 1)
other = sp.expand(coeff4(sp.expand(Uz * power_d(wz))))
other = sp.Poly(other, al)

# Power sums p_k = sum over zeta != 1 of alpha^k.
# d/ds log(((1+s)^q - 1)/s) = sum_{zeta != 1} alpha/(1 + alpha s) = sum_k (-s)^k p_{k+1}.
K = other.degree() + 1
base = sum(sp.expand_func(sp.binomial(q, j)) * s ** (j - 1) for j in range(1, K + 3))
dlog = sp.series(sp.diff(sp.log(base), s), s, 0, K + 1).removeO()
p = {k: sp.factor(sp.expand(dlog).coeff(s, k - 1) * (-1) ** (k - 1)) for k in range(1, K + 1)}
p[0] = q - 1
assert sp.simplify(p[1] - (q - 1) / 2) == 0

other_sum = sum(c * p[m[0]] for m, c in zip(other.monoms(), other.coeffs()))

P = -1 + 4 * d * q + (d + 5 * d * (d - 1) / 2) * q ** 2 - d * q ** 3
T = (d * (d - 1) * (d - 2) * q ** 6 + d * (d - 1) * (d - 11) * q ** 5 + (-19 * d ** 3 + 18 * d ** 2 + 19 * d) * q ** 4
     + (26 * d ** 3 + 18 * d ** 2 + 4 * d - 6) * q ** 3 + (36 * d ** 2 + 12 * d - 6) * q ** 2 + (18 * d - 6) * q - 6)

assert sp.simplify(4320 * q ** 4 * principal - 6 * P) == 0
print("principal root: 4320 q^4 [v^4](...) = 6 (1+q)^d P")
assert sp.simplify(4320 * q ** 4 * other_sum - (q - 1) * T) == 0
print("other roots:    4320 q^4 sum [v^4](...) = (q-1) T")
print("codimension-four formula holds identically in q and d")
