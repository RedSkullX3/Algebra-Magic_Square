"""
check_3x3.py

Checks every formula and theorem about 3 x 3 magic squares in the paper.

Each formula is tested on hundreds of randomly chosen magic squares, using
exact fractions, so there is no rounding. A formula that holds in all these
random cases (with no exceptions) is the computer's confirmation of the
algebra done by hand in the paper.

Run with:   python3 check_3x3.py
"""

import random
from fractions import Fraction
from magic_tools import (zero_matrix, identity, all_ones, reversal, add,
                         subtract, scale, combine, multiply, multiply_all,
                         power, equal, trace, is_magic, is_semi_magic,
                         is_zero_magic, line_sum, determinant, inverse,
                         characteristic_polynomial, multiply_polynomials,
                         minimal_polynomial_degree, zero_magic_basis,
                         check, summary, show)

rng = random.Random(2026)          # fixed seed, so the run is repeatable
J = all_ones(3)
I = identity(3)
P = reversal(3)
TRIALS = 300


def Z(x, y):
    """The zero-magic square with Ward coordinates (x, y)."""
    return [[x, y, -x - y],
            [-2 * x - y, 0, 2 * x + y],
            [x + y, -y, -x]]


U = Z(1, 0)          # Ward's basis
V = Z(0, 1)


def magic(alpha, beta, c):
    """The magic square alpha*U + beta*V + c*J (line-sum 3c)."""
    return add(Z(alpha, beta), scale(c, J))


def r():
    """A random whole number between -9 and 9."""
    return rng.randint(-9, 9)


def r_nonzero():
    x = 0
    while x == 0:
        x = r()
    return x


print("Checks for 3 x 3 magic squares")
print("==============================")

# ----------------------------------------------------------------------
print("\nWard's theorem (dimension of the zero-magic squares)")
check("dimension of 0MS(n) is n^2 - 2n - 1 for n = 3, ..., 8",
      all(len(zero_magic_basis(n)) == n * n - 2 * n - 1 for n in range(3, 9)))
check("every 3x3 magic square is alpha*U + beta*V + c*J "
      "(0MS(3) has dimension 2)", len(zero_magic_basis(3)) == 2)

# ----------------------------------------------------------------------
print("\nThe centre and the canonical form M(m, a, b)")
ok_centre = True
for _ in range(TRIALS):
    A = magic(r(), r(), r())
    m = A[1][1]
    if line_sum(A) != 3 * m:
        ok_centre = False
    for i in range(3):
        for j in range(3):
            if A[i][j] + A[2 - i][2 - j] != 2 * m:
                ok_centre = False
check("line-sum = 3 x centre, and opposite cells add to 2 x centre", ok_centre)


def M_canonical(m, a, b):
    return [[m + a, m - a - b, m + b],
            [m - a + b, m, m + a - b],
            [m - b, m + a + b, m - a]]


check("M(m, a, b) is always magic with line-sum 3m",
      all(is_magic(M_canonical(m, a, b)) and line_sum(M_canonical(m, a, b)) == 3 * m
          for m, a, b in [(r(), r(), r()) for _ in range(TRIALS)]))
lo_shu = [[2, 9, 4], [7, 5, 3], [6, 1, 8]]
check("the Lo Shu square is M(5, -3, -1)", equal(M_canonical(5, -3, -1), lo_shu))

# ----------------------------------------------------------------------
print("\nThe general basis U', V' (Theorem 3) and the table of bases")
ok = True
for _ in range(TRIALS):
    p, q, rr, s, a, b = r(), r(), r(), r(), r(), r()
    Up, Vp = Z(p, q), Z(rr, s)
    if not (is_zero_magic(Up) and is_zero_magic(Vp)):
        ok = False
    if not equal(Up, combine([(p, U), (q, V)])):
        ok = False
    if not equal(combine([(a, Up), (b, Vp)]), Z(a * p + b * rr, a * q + b * s)):
        ok = False
check("U' = pU + qV is zero-magic, and aU' + bV' has Ward coordinates "
      "(ap + br, aq + bs)", ok)
thompson = ([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], [[1, -1, 0], [-1, 0, 1], [0, 1, -1]])
hill_elzaidi = ([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], [[-1, 1, 0], [1, 0, -1], [0, -1, 1]])
gauthier = ([[1, -1, 0], [-1, 0, 1], [0, 1, -1]], [[0, -1, 1], [1, 0, -1], [-1, 1, 0]])
check("table: Ward (1,0),(0,1)", equal(Z(1, 0), U) and equal(Z(0, 1), V))
check("table: Thompson (0,1),(1,-1)",
      equal(Z(0, 1), thompson[0]) and equal(Z(1, -1), thompson[1]))
check("table: Hill-Elzaidi (0,1),(-1,1)",
      equal(Z(0, 1), hill_elzaidi[0]) and equal(Z(-1, 1), hill_elzaidi[1]))
check("table: Gauthier (1,-1),(0,-1), up to his scaling by 1/sqrt(3)",
      equal(Z(1, -1), gauthier[0]) and equal(Z(0, -1), gauthier[1]))

# ----------------------------------------------------------------------
print("\nThe 2 x 2 picture (equation (5))")
S = [[1, 1, 1],
     [1, -2, 0],
     [1, 1, -1]]      # columns are (1,1,1), u = (1,-2,1), w = (1,0,-1)
S_inv = inverse(S)


def conjugate(A):
    """The matrix S^(-1) A S."""
    return multiply(multiply(S_inv, A), S)


ok = True
for _ in range(TRIALS):
    alpha, beta, c = r(), r(), r()
    x, y = 2 * alpha + beta, -3 * beta
    expected = [[3 * c, 0, 0], [0, 0, x], [0, y, 0]]
    if not equal(conjugate(magic(alpha, beta, c)), expected):
        ok = False
check("S^(-1) A S = (3c) together with [[0, x], [y, 0]], "
      "x = 2alpha + beta, y = -3beta", ok)
check("S^(-1) J S = diag(3, 0, 0)", equal(conjugate(J), [[3, 0, 0], [0, 0, 0], [0, 0, 0]]))
check("Lo Shu square: the 2 x 2 block is [[0, -2], [-12, 0]]",
      equal([row[1:] for row in conjugate(lo_shu)[1:]], [[0, -2], [-12, 0]]))
ok = all(is_magic(multiply(multiply(S, [[t1, 0, 0], [0, 0, t2], [0, t3, 0]]), S_inv))
         for t1, t2, t3 in [(r(), r(), r()) for _ in range(TRIALS)])
check("converse: every 'number + anti-diagonal block' comes from a magic square", ok)
diagonals = [conjugate(M) for M in (I, J, P)]
is_diag = all(M[i][j] == 0 for M in diagonals for i in range(3) for j in range(3) if i != j)
independent = determinant([[M[0][0], M[1][1], M[2][2]] for M in diagonals]) != 0
check("I, J and P become diagonal, and together give all diagonal matrices",
      is_diag and independent)

# the bullet points
ok = True
for _ in range(TRIALS):
    x, y = r(), r()
    block = [[0, x], [y, 0]]
    if not equal(multiply(block, block), [[x * y, 0], [0, x * y]]):
        ok = False
check("[[0, x], [y, 0]] squared is xy times the identity", ok)
ok = True
for _ in range(100):
    alpha, beta, c = r(), r(), r()
    A = magic(alpha, beta, c)
    A0 = Z(alpha, beta)
    xy = (2 * alpha + beta) * (-3 * beta)
    for k in range(5):
        if not equal(power(A, 2 * k + 1), add(scale(xy ** k, A0), scale(3 ** (2 * k) * c ** (2 * k + 1), J))):
            ok = False
check("A^(2k+1) = (xy)^k A0 + 3^(2k) c^(2k+1) J  for k = 0, ..., 4", ok)
ok = True
for _ in range(TRIALS):
    a1, b1, a2, b2, c1, c2 = r(), r(), r(), r(), r(), r()
    x1, y1, x2, y2 = 2 * a1 + b1, -3 * b1, 2 * a2 + b2, -3 * b2
    A, B = magic(a1, b1, c1), magic(a2, b2, c2)
    commute = equal(multiply(A, B), multiply(B, A))
    if commute != (x1 * y2 == x2 * y1):
        ok = False
check("two magic squares commute exactly when xy' = x'y", ok)
ok = True
for _ in range(TRIALS):
    alpha, beta = r(), r()
    nilpotent = equal(power(Z(alpha, beta), 3), zero_matrix(3))
    if nilpotent != (beta == 0 or beta == -2 * alpha):
        ok = False
check("the zero-magic part is nilpotent exactly when beta = 0 or beta = -2alpha", ok)

# Cook, Bacon and Hillman's form of the odd-power formula (their Prop. 2.1)
ok = True
for _ in range(100):
    a, b = r(), r()
    Sc = [[a, -(a + b), b], [-(a - b), 0, a - b], [-b, a + b, -a]]
    for n in range(4):
        if not equal(power(Sc, 2 * n + 1), scale(3 ** n * (a * a - b * b) ** n, Sc)):
            ok = False
check("Cook-Bacon-Hillman: S^(2n+1) = 3^n (a^2 - b^2)^n S (the same formula)", ok)

# ----------------------------------------------------------------------
print("\nCharacteristic polynomial and determinant (Theorem 4)")
ok_char = ok_det = ok_amir = True
for _ in range(TRIALS):
    p, q, rr, s, a, b, c = r(), r(), r(), r(), r(), r(), r()
    A = combine([(a, Z(p, q)), (b, Z(rr, s)), (c, J)])
    alpha, beta = a * p + b * rr, a * q + b * s
    expected = multiply_polynomials([1, -3 * c], [1, 0, 3 * beta * (beta + 2 * alpha)])
    if characteristic_polynomial(A) != expected:
        ok_char = False
    if determinant(A) != 9 * c * beta * (beta + 2 * alpha):
        ok_det = False
    Aw = magic(a, b, c)
    amir = [1, -3 * c, 3 * b * (b + 2 * a), -9 * b * c * (b + 2 * a)]
    if characteristic_polynomial(Aw) != amir:
        ok_amir = False
check("det(xI - A) = (x - 3c)(x^2 + 3beta(beta + 2alpha))", ok_char)
check("Amir-Moez and Fredricks' form in Ward's basis", ok_amir)
check("det(A) = 9c beta(beta + 2alpha)", ok_det)
ok = True
for _ in range(100):
    m, u, v = r(), r(), r()
    lucas = [[m + u, m + v - u, m - v], [m - u - v, m, m + u + v], [m + v, m + u - v, m - u]]
    if not equal(lucas, M_canonical(m, u, -v)):
        ok = False
    turned = [[lucas[2 - j][i] for j in range(3)] for i in range(3)]   # quarter turn
    if characteristic_polynomial(lucas) != multiply_polynomials([1, -3 * m], [1, 0, -3 * u * u + 3 * v * v]):
        ok = False
    if characteristic_polynomial(turned) != multiply_polynomials([1, -3 * m], [1, 0, 3 * u * u - 3 * v * v]):
        ok = False
check("Lucas's form is M(m, u, -v), and a quarter turn swaps real and "
      "imaginary eigenvalues (Loly et al.)", ok)

# ----------------------------------------------------------------------
print("\nProducts of two magic squares")
ok_sums = ok_magic = True
nilpotent_cases = 0
for _ in range(2000):
    a, b, c, d, e, f = r(), r(), r(), r(), r(), r()
    if rng.random() < 0.3:          # make some special cases where AB can be magic
        choice = rng.randint(0, 2)
        if choice == 0:
            b, e = 0, 0
        elif choice == 1:
            b, e = -2 * a, -2 * d
        else:
            d, e = 0, 0
    A, B = magic(a, b, c), magic(d, e, f)
    AB = multiply(A, B)
    if not (is_semi_magic(AB) and line_sum(AB) == 9 * c * f):
        ok_sums = False
    if trace(AB) != 9 * c * f - 6 * (a * e + b * d + b * e):
        ok_sums = False
    if AB[0][2] + AB[1][1] + AB[2][0] != 9 * c * f + 6 * (b * d - a * e):
        ok_sums = False
    if is_magic(AB):
        nilpotent_cases += 1
        if not equal(AB, scale(3 * c * f, J)):
            ok_magic = False
check("rows and columns of AB sum to 9cf, and the diagonals are "
      "9cf - 6(ae + bd + be) and 9cf + 6(bd - ae)", ok_sums)
check("whenever AB is magic, AB = 3cf J (" + str(nilpotent_cases) + " such cases found)",
      ok_magic and nilpotent_cases > 0)

# ----------------------------------------------------------------------
print("\nCubes and odd products (Theorems 5 and 6)")
ok_cube = True
for _ in range(TRIALS):
    alpha, beta, c = r(), r(), r()
    A0 = Z(alpha, beta)
    if not equal(power(A0, 3), scale(-3 * beta * (beta + 2 * alpha), A0)):
        ok_cube = False
    if not is_magic(power(magic(alpha, beta, c), 3)):
        ok_cube = False
check("A0^3 = -3beta(beta + 2alpha) A0, and every cube is magic", ok_cube)
check("J^3 = 9J", equal(power(J, 3), scale(9, J)))
ok = all(is_magic(multiply_all([magic(r(), r(), r()) for _ in range(k)]))
         for k in (3, 5, 7) for _ in range(100))
check("products of 3, 5 and 7 magic squares are always magic", ok)
not_magic = sum(1 for k in (2, 4) for _ in range(100)
                if not is_magic(multiply_all([magic(r_nonzero(), r_nonzero(), r()) for _ in range(k)])))
check("products of 2 or 4 nonconstant squares are usually not magic "
      "(" + str(not_magic) + " of 200 were not)", not_magic >= 190)

# ----------------------------------------------------------------------
print("\nInverses (Theorem 8)")
ok6 = ok_sq = ok_A0 = ok_simple = True
tested = 0
while tested < TRIALS:
    alpha, beta, c = r(), r(), r()
    k = 3 * beta * (beta + 2 * alpha)
    if c == 0 or k == 0:
        continue
    tested += 1
    A = magic(alpha, beta, c)
    A0 = Z(alpha, beta)
    A_inv = inverse(A)
    formula6 = scale(Fraction(1, 9 * c * beta * (beta + 2 * alpha)),
                     add(subtract(power(A, 2), scale(3 * c, A)), scale(k, I)))
    if not equal(A_inv, formula6):
        ok6 = False
    if not equal(power(A, 2), add(scale(-k, I), scale(3 * c * c + Fraction(k, 3), J))):
        ok_sq = False
    if not equal(power(A0, 2), add(scale(-k, I), scale(Fraction(k, 3), J))):
        ok_A0 = False
    simple = subtract(scale(Fraction(1, 9 * c), J), scale(Fraction(1, k), A0))
    if not (equal(A_inv, simple) and is_magic(A_inv)):
        ok_simple = False
check("equation (6): A^(-1) = (A^2 - 3cA + 3beta(beta+2alpha) I) / (9c beta(beta+2alpha))", ok6)
check("A^2 = -3b(2a + b) I + (3c^2 + b(2a + b)) J", ok_sq)
check("A0^2 = -kI + (k/3) J  with k = 3beta(beta + 2alpha)", ok_A0)
check("A^(-1) = J/(9c) - A0/k, and A^(-1) is magic", ok_simple)

# ----------------------------------------------------------------------
print("\nThe minimal polynomial")
cases = [
    ("a typical square", magic(1, 1, 1), 3),
    ("another typical square", magic(2, 3, 5), 3),
    ("beta = 0, c not 0, not constant (still the full cubic)", magic(1, 0, 1), 3),
    ("beta = -2alpha, c not 0 (still the full cubic)", magic(1, -2, 2), 3),
    ("c = 0 and beta = 0", magic(1, 0, 0), 2),
    ("c = 0 and beta = -2alpha", magic(1, -2, 0), 2),
    ("a constant matrix", magic(0, 0, 3), 2),
    ("beta(beta + 2alpha) = -3c^2, here -3 with c = 1", magic(2, -1, 1), 2),
    ("c = 0 but beta(beta + 2alpha) not 0", magic(1, 1, 0), 3),
]
for description, A, degree in cases:
    check("degree " + str(degree) + " for " + description,
          minimal_polynomial_degree(A) == degree)

# ----------------------------------------------------------------------
print("\nCommuting squares (Theorem 9)")
K = [[1, 1, -2], [1, -2, 1], [-2, 1, 1]]
ok = True
for _ in range(TRIALS):
    p, q, rr, s, a, b, c, d, e, f = [r() for _ in range(10)]
    A = combine([(a, Z(p, q)), (b, Z(rr, s)), (c, J)])
    B = combine([(d, Z(p, q)), (e, Z(rr, s)), (f, J)])
    if not equal(subtract(multiply(A, B), multiply(B, A)),
                 scale(2 * (p * s - q * rr) * (a * e - b * d), K)):
        ok = False
check("AB - BA = 2(ps - qr)(ae - bd) K", ok)

# ----------------------------------------------------------------------
print("\nOrderings of a triple product (Theorem 10)")


def Mtri(x, y):
    return [[x + y, -y, -x], [-2 * x - y, 0, 2 * x + y], [x, y, -x - y]]


ok1 = ok2 = ok3 = True
for _ in range(TRIALS):
    a, b, c, d, e, f, g, h, i = [r() for _ in range(9)]
    A, B, C = magic(a, b, c), magic(d, e, f), magic(g, h, i)
    ABC = multiply_all([A, B, C])
    if not equal(subtract(ABC, multiply_all([A, C, B])), scale(6 * (d * h - e * g), Mtri(a, b))):
        ok1 = False
    if not equal(subtract(ABC, multiply_all([B, A, C])), scale(-6 * (a * e - b * d), Mtri(g, h))):
        ok2 = False
    diff3 = subtract(multiply_all([A, C, B]), multiply_all([B, A, C]))
    if not (equal(diff3, scale(6 * (a * h - b * g), Mtri(d, e))) or
            equal(diff3, scale(-6 * (a * h - b * g), Mtri(d, e)))):
        ok3 = False
check("ABC - ACB = 6(dh - eg) M(A)", ok1)
check("ABC - BAC = -6(ae - bd) M(C)", ok2)
check("ACB - BAC is a multiple of (ah - bg) times M(B)", ok3)
check("M(x, y) is zero only when x = y = 0",
      all((equal(Mtri(x, y), zero_matrix(3))) == (x == 0 and y == 0)
          for x in range(-3, 4) for y in range(-3, 4)))
values = [-1, 0, 1, 2]
ok = True
count = 0
for a in values:
    for b in values:
        for d in values:
            for e in values:
                for g in values:
                    for h in values:
                        A, B, C = magic(a, b, 1), magic(d, e, 2), magic(g, h, 3)
                        products = [multiply_all(order) for order in
                                    ([A, B, C], [A, C, B], [B, A, C], [B, C, A], [C, A, B], [C, B, A])]
                        all_equal = all(equal(X, products[0]) for X in products)
                        proportional = (a * e == b * d and a * h == b * g and d * h == e * g)
                        constant = ((a, b) == (0, 0) or (d, e) == (0, 0) or (g, h) == (0, 0))
                        if all_equal != (proportional or constant):
                            ok = False
                        count += 1
check("the 'if and only if' statement holds for all " + str(count) + " small test triples", ok)

summary()
