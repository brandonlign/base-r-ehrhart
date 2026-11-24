# Ehrhart positivity of base-$r$ simplices

Brandon Li

For an integer $r\ge2$, Solus's base-$r$ simplex $B_{r,d}$ is the convex hull
of the standard basis vectors of $\mathbb R^d$ and the point
$-(r-1)(1,r,\ldots,r^{d-1})$. Solus observed that these simplices appeared to
be Ehrhart positive; Liu recorded this as a conjecture, and Liu and Solus
stated it for all $r\ge2$ and $d\ge1$. P. Graham showed that $B_{14,5}$ is a
counterexample
([Hexagon 2611.00793](https://hexagonmath.org/2611.00793v1)). This paper
determines which bases are Ehrhart positive.

**Theorem.** $B_{r,d}$ is Ehrhart positive for every $d\ge1$ if and only if
$r\in\{2,3,4,5\}$. For $r\ge6$ the first dimension with a negative
coefficient is

| $r$ | 6 | 7 | 8 | 9–13 | ≥ 14 |
|---|---|---|---|---|---|
| first $d$ | 43 | 19 | 15 | 7 | 5 |

and in that dimension the negative coefficient is the linear one.

The paper also proves:

- at fixed $d$ and large $r$, the coefficient of $t^k$ is negative exactly
  when $k\ge1$, $d-k\ge4$ and $d-k\equiv0,1\pmod 4$;
- for $d\ge6$ the coefficient of $t^{d-4}$ is negative exactly when
  $r-1>(5d-3)/2$;
- at fixed $r$, coefficients in fixed codimension are eventually positive,
  with an explicit limit when $r/d$ tends to a constant.

Everything starts from $L_d(t)=[x^{(r-1)t}]F_d(x)$ for an explicit rational
function $F_d$, and a residue formula over the $(r-1)$st roots of unity. For
bases 4 and 5 the generating function over all dimensions splits into
branches of an algebraic function. The principal branch is shown to
dominate, using an exact Lagrange-inversion certificate for the first 180 or
900 coefficients and Cauchy estimates on a disk for the rest.

- [Paper (PDF)](paper/main.pdf)
- [LaTeX source](paper/)
- [Verification code](checks/)

## Reproduce

```sh
bash checks/run_all.sh
```

This needs only Python 3 (exact integer and rational arithmetic) and takes
about fifteen seconds. It runs:

- `ehrhart.py`: $L_d(t)$ three ways (Solus's $h^*$-polynomial, the
  convolution formula of Section 2, and direct lattice-point counts),
  including Graham's $L_5$ for $r=14$.
- `rows.py`: the large-base signs (Theorem 1.2) for $d\le16$, the
  codimension-four formula, its finite rows and the bound in its proof (Theorem 1.3), the formulas for
  $d\le5$, and the first-failure table.
- `bases_2_3.py`: the generating functions and coefficient facts for $r=2,3$
  (Section 6).
- `bases_4_5.py`: the Lagrange-inversion certificate and tail inequalities for
  $r=4,5$ (Section 7).
- `disks.py`: the critical-value bounds for the analytic continuation in
  Section 7.

The output is in [`checks/output.txt`](checks/output.txt).
`codim4_symbolic.py` (needs `sympy`, a few minutes) derives the
codimension-four formula of Proposition 4.1 symbolically in $q=r-1$ and $d$. The
PDF builds with Tectonic (`cd paper && tectonic main.tex`) or `pdflatex`.

## License

Code: MIT. Paper: CC BY 4.0.
