"""
check_odd_powers.py

Checks the last section of the paper, "Odd powers in every order":

  * Theorem 12 (a finite test for when every odd power is magic),
  * Corollary 13 (order four: a magic cube exactly when the square is singular,
    and the formula for all odd powers),
  * every example used there (D, the Durer-type square F175, Freitag's
    Fibonacci square, Hill and Elzaidi's matrix, the 5 x 5 square C, and the
    squares V, W, Z of Cook, Bacon and Hillman),
  * the high-precision check in orders 6 and 7.

Run with:   python3 check_odd_powers.py      (takes a few seconds)
"""

import random
from decimal import Decimal, getcontext
from fractions import Fraction
from magic_tools import (all_ones, add, scale, combine, power, equal, trace,
                         anti_trace, line_sums, is_magic, is_zero_magic,
                         is_associative, is_pandiagonal, line_sum,
                         determinant, characteristic_polynomial,
                         polynomial_to_text, zero_magic_basis,
                         associative_zero_magic_basis, line_sum_equations,
                         solve_homogeneous, vector_to_matrix,
                         place_in_rows_and_columns, check, summary, show)

rng = random.Random(12)

print("Odd powers in every order")
print("=========================")


def odd_powers_magic(A, up_to=9):
    """True if A^1, A^3, ..., A^up_to are all magic."""
    return all(is_magic(power(A, k)) for k in range(1, up_to + 1, 2))


def odd_coefficients_vanish(A):
    """True if the eigenvalues other than the line-sum come in pairs +mu, -mu.
    We divide det(xI - A) by (x - m) and check that the result has only
    even powers of x, or only odd powers of x."""
    m = line_sum(A)
    p = characteristic_polynomial(A)
    quotient = []                       # synthetic division by (x - m)
    carry = Fraction(0)
    for c in p[:-1]:
        carry = carry * m + c
        quotient.append(carry)
    degree = len(quotient) - 1
    return all(c == 0 for i, c in enumerate(quotient) if (degree - i) % 2 == 1)


# ----------------------------------------------------------------------
print("\nCorollary 13 (order four)")
f175 = [[1, 12, 8, 13], [14, 7, 11, 2], [15, 6, 10, 3], [4, 9, 5, 16]]
D = [[1, 6, 15, 12], [11, 16, 5, 2], [8, 3, 10, 13], [14, 9, 4, 7]]
J4 = all_ones(4)
for name, A in (("F175", f175), ("D", D)):
    m = 34
    lam2 = Fraction(trace(power(A, 2)) - m * m, 2)
    ok = all(equal(power(A, 2 * k + 1),
                   add(scale(lam2 ** k, A), scale(Fraction(m ** (2 * k + 1) - lam2 ** k * m, 4), J4)))
             for k in range(0, 6))
    check(name + ": singular, lambda^2 = " + str(lam2) + ", and "
          "A^(2k+1) = lambda^(2k) A + (m^(2k+1) - lambda^(2k) m)/4 J for k = 0..5",
          determinant(A) == 0 and ok)

H = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]]
check("Hill and Elzaidi's matrix: determinant 1, cube not magic",
      determinant(H) == 1 and not is_magic(power(H, 3)))

F5 = [[13, 89, 97, 34], [110, 21, 63, 39], [68, 94, 55, 16], [42, 29, 18, 144]]
check("Freitag's Fibonacci square: magic, not associative, determinant 88,134,813, "
      "cube not magic",
      is_magic(F5) and not is_associative(F5) and determinant(F5) == 88134813
      and not is_magic(power(F5, 3)))

# random 4 x 4 magic squares: ordinary ones, associative ones and
# pandiagonal ones (the last two kinds are always singular)
basis4 = zero_magic_basis(4)
assoc4 = associative_zero_magic_basis(4)
equations = line_sum_equations(4)
for shift in range(4):              # all broken diagonals sum to 0
    equations.append([1 if k % 4 == (k // 4 + shift) % 4 else 0 for k in range(16)])
    equations.append([1 if k % 4 == (shift - k // 4) % 4 else 0 for k in range(16)])
pandiagonal4 = [vector_to_matrix(v, 4) for v in solve_homogeneous(equations, 16)]
check("found " + str(len(pandiagonal4)) + " basic zero-magic pandiagonal 4 x 4 squares",
      all(is_zero_magic(B) and is_pandiagonal(B) for B in pandiagonal4))
agree = True
singular_count = 0
total = 0
for kind in (basis4, assoc4, pandiagonal4):
    for _ in range(100):
        A = combine([(rng.randint(-5, 5), B) for B in kind] + [(rng.randint(1, 5), J4)])
        singular = (determinant(A) == 0)
        singular_count += singular
        total += 1
        if singular != is_magic(power(A, 3)) or singular != odd_powers_magic(A, 7):
            agree = False
check("300 random 4 x 4 magic squares (" + str(singular_count) + " singular): "
      "singular <=> cube magic <=> all odd powers magic", agree)

# ----------------------------------------------------------------------
print("\nMattingly: associative squares of even order are singular")
ok = True
for n in (4, 6):
    basis = associative_zero_magic_basis(n)
    for _ in range(20):
        A = combine([(rng.randint(-5, 5), B) for B in basis] + [(rng.randint(1, 5), all_ones(n))])
        if determinant(A) != 0:
            ok = False
check("20 random associative squares each of orders 4 and 6 have determinant 0", ok)

# ----------------------------------------------------------------------
print("\nTheorem 12 in odd orders: associative squares (Nordgren)")
ok = True
for n in (5, 7):
    basis = associative_zero_magic_basis(n)
    for _ in range(10):
        A = combine([(rng.randint(-3, 3), B) for B in basis] + [(rng.randint(1, 3), all_ones(n))])
        if not (odd_powers_magic(A, 7) and odd_coefficients_vanish(A)):
            ok = False
check("random associative squares of orders 5 and 7: eigenvalues in pairs and "
      "all odd powers (to the 7th) magic", ok)

# ----------------------------------------------------------------------
print("\nRemark 14: order five")
V = [[4, 13, 19, 28, 36], [26, 35, 8, 11, 20], [15, 18, 27, 33, 7],
     [34, 5, 14, 22, 25], [21, 29, 32, 6, 12]]
check("Cook-Bacon-Hillman's 5 x 5 square V: magic, not associative, cube magic, "
      "and then A^5, A^7, A^9 magic too",
      is_magic(V) and not is_associative(V) and odd_powers_magic(V, 9))
print("      det(xI - V) = " + polynomial_to_text(characteristic_polynomial(V)))

C = [[7, 8, 5, 6, 9],
     [6, 8, 9, 5, 7],
     [8, 5, 5, 9, 8],
     [6, 8, 7, 9, 5],
     [8, 6, 9, 6, 6]]
show(C)
C3 = power(C, 3)
char = characteristic_polynomial(C)
print("      det(xI - C) = " + polynomial_to_text(char))
check("C is magic with line-sum 35", is_magic(C) and line_sum(C) == 35)
check("det(xI - C) = (x - 35)(x^4 - 3x^2 - 50), so the other eigenvalues "
      "come in pairs", char == [1, -35, -3, 105, -50, 1750])
check("main diagonal of C^3 sums to 35^3, but the anti-diagonal sums to 35^3 + 40",
      trace(C3) == 35 ** 3 and anti_trace(C3) == 35 ** 3 + 40)
C0 = add(C, scale(-7, all_ones(5)))
for n in (5, 7, 9, 11):
    positions = [0, 1, (n - 1) // 2, n - 2, n - 1]
    Cn = place_in_rows_and_columns(C0, positions, n)
    check("order " + str(n) + ": C - 7J placed in rows and columns 1, 2, (n+1)/2, n-1, n "
          "is zero-magic, eigenvalues in pairs, cube anti-diagonal = 40",
          is_zero_magic(Cn) and odd_coefficients_vanish(Cn) and anti_trace(power(Cn, 3)) == 40)

# ----------------------------------------------------------------------
print("\nCook, Bacon and Hillman's squares of orders 6 and 7")
W = [[31, 9, 2, 22, 27, 20], [3, 32, 7, 21, 23, 25], [8, 28, 6, 26, 19, 24],
     [35, 1, 33, 17, 10, 15], [30, 5, 34, 12, 14, 16], [4, 36, 29, 13, 18, 11]]
check("W (6 x 6): magic, eigenvalues not in pairs, cube not magic",
      is_magic(W) and not odd_coefficients_vanish(W) and not is_magic(power(W, 3)))
Z7 = [[59, 99, 167, 11, 79, 147, 187], [91, 131, 199, 43, 111, 151, 23],
      [95, 163, 35, 75, 143, 183, 55], [127, 195, 39, 107, 175, 19, 87],
      [159, 31, 71, 139, 179, 51, 119], [191, 63, 103, 171, 15, 83, 123],
      [27, 67, 135, 203, 47, 115, 155]]
check("Z (7 x 7): magic and associative, so all odd powers (to the 9th) are magic",
      is_magic(Z7) and is_associative(Z7) and odd_powers_magic(Z7, 9))

# ----------------------------------------------------------------------
print("\nOrders 6 and 7: a magic cube but a non-magic fifth power")
print("      (solved with 80-digit decimal arithmetic)")
getcontext().prec = 80


def find_example(n, seed):
    """Look for a zero-magic B = B0 + s*D1 + t*D2 whose cube is zero-magic.
    Newton's method is used to solve trace(B^3) = 0 and anti_trace(B^3) = 0."""
    basis = zero_magic_basis(n)
    picker = random.Random(seed)
    B0 = combine([(picker.randint(-3, 3), M) for M in basis])
    D1, D2 = basis[0], basis[1]

    def build(s, t):
        return add(B0, add(scale(s, D1), scale(t, D2)))

    def F(s, t):
        cube = power(build(s, t), 3)
        return trace(cube), anti_trace(cube)

    h = Decimal("1e-30")
    for start in [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1), (1, 1)]:
        s, t = Decimal(start[0]), Decimal(start[1])
        for _ in range(100):
            f1, f2 = F(s, t)
            if abs(f1) + abs(f2) < Decimal("1e-60"):
                return build(s, t)
            a1, a2 = F(s + h, t)
            b1, b2 = F(s, t + h)
            j11, j21 = (a1 - f1) / h, (a2 - f2) / h
            j12, j22 = (b1 - f1) / h, (b2 - f2) / h
            det = j11 * j22 - j12 * j21
            if det == 0:
                break
            s = s - (f1 * j22 - f2 * j12) / det
            t = t - (j11 * f2 - j21 * f1) / det
            if abs(s) > 10 ** 6 or abs(t) > 10 ** 6:
                break
    return None


for n in (6, 7):
    B = None
    seed = 1
    while B is None:
        B = find_example(n, seed)
        seed += 1
    cube = power(B, 3)
    fifth = power(B, 5)
    biggest = max(abs(x) for row in B for x in row)
    print("      order " + str(n) + ": largest entry of B is about " + format(biggest, ".4g"))
    print("        cube:        main diagonal " + format(trace(cube), ".2e") +
          ",  anti-diagonal " + format(anti_trace(cube), ".2e"))
    print("        fifth power: main diagonal " + format(trace(fifth), ".6g") +
          ",  anti-diagonal " + format(anti_trace(fifth), ".6g"))
    tiny = Decimal("1e-50")
    check("order " + str(n) + ": B is zero-magic, B^3 is zero-magic, but B^5 is not",
          all(abs(x) < tiny for x in line_sums(B)) and
          abs(trace(cube)) < tiny and abs(anti_trace(cube)) < tiny and
          (abs(trace(fifth)) > 1 or abs(anti_trace(fifth)) > 1))

summary()
