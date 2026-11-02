"""
classical_4x4.py

Finds every classical 4 x 4 magic square (the numbers 1, 2, ..., 16 each
used once, every row, column and diagonal adding up to 34) and checks the
paper's main example:

    A 4 x 4 magic square has a magic cube exactly when it is singular
    (determinant 0). Exactly 640 of the 880 classical squares are singular,
    so exactly 640 have magic cubes.

Run with:   python3 classical_4x4.py        (takes a few seconds)
"""

import time
from magic_tools import (power, equal, all_ones, add, scale, is_magic,
                         is_associative, is_pandiagonal, determinant,
                         characteristic_polynomial, polynomial_to_text,
                         check, summary, show)

# Name the 16 cells like this:
#     a b c d
#     e f g h
#     i j k l
#     m n o p


def find_standard_squares():
    """All classical 4 x 4 magic squares in 'standard position': the
    smallest corner is at the top left (a), and d < m. Every classical
    square can be turned or flipped into exactly one such square."""
    squares = []
    numbers = range(1, 17)

    def ok(x, used):
        return 1 <= x <= 16 and x not in used

    for a in numbers:
        for p in range(a + 1, 17):
            for f in numbers:
                k = 34 - a - f - p                     # main diagonal
                if f in (a, p) or not ok(k, {a, p, f}):
                    continue
                for d in range(a + 1, 17):
                    if d in (p, f, k):
                        continue
                    for m in range(d + 1, 17):
                        if m in (p, f, k):
                            continue
                        used = {a, f, k, p, d, m}
                        for g in numbers:
                            j = 34 - d - g - m             # anti-diagonal
                            if g in used or not ok(j, used | {g}):
                                continue
                            used2 = used | {g, j}
                            for b in numbers:
                                c = 34 - a - b - d         # row 1
                                if b in used2 or not ok(c, used2 | {b}):
                                    continue
                                used3 = used2 | {b, c}
                                n = 34 - b - f - j         # column 2
                                o = 34 - c - g - k         # column 3
                                if not ok(n, used3) or not ok(o, used3 | {n}):
                                    continue
                                if m + n + o + p != 34:    # row 4
                                    continue
                                used4 = used3 | {n, o}
                                for e in numbers:
                                    i = 34 - a - e - m     # column 1
                                    h = 34 - e - f - g     # row 2
                                    l = 34 - d - h - p     # column 4
                                    if e in used4 or not ok(i, used4 | {e}):
                                        continue
                                    if not ok(h, used4 | {e, i}):
                                        continue
                                    if not ok(l, used4 | {e, i, h}):
                                        continue
                                    if i + j + k + l != 34:  # row 3
                                        continue
                                    squares.append([[a, b, c, d],
                                                    [e, f, g, h],
                                                    [i, j, k, l],
                                                    [m, n, o, p]])
    return squares


def rotate(A):
    """Turn the square a quarter turn clockwise."""
    return [[A[3 - j][i] for j in range(4)] for i in range(4)]


def flip(A):
    """Reflect the square in its main diagonal."""
    return [[A[j][i] for j in range(4)] for i in range(4)]


def all_eight_versions(A):
    """The 8 rotations and reflections of a square."""
    versions = []
    B = A
    for _ in range(4):
        versions.append(B)
        versions.append(flip(B))
        B = rotate(B)
    return versions


print("All classical 4 x 4 magic squares")
print("=================================")
start = time.time()
standard = find_standard_squares()
print("\nSearch finished in " + str(round(time.time() - start, 1)) + " seconds.")

everything = []
for A in standard:
    everything.extend(all_eight_versions(A))
distinct = {tuple(x for row in A for x in row) for A in everything}

check("there are 880 classical squares up to rotations and reflections",
      len(standard) == 880)
check("there are 7040 classical squares in total (880 x 8), all different",
      len(everything) == 7040 and len(distinct) == 7040)
check("every square found really is magic with line-sum 34",
      all(is_magic(A) and sum(A[0]) == 34 for A in everything))

# ----------------------------------------------------------------------
print("\nSingular squares and magic cubes")
singular_standard = 0
agree = True
for A in everything:
    singular = (determinant(A) == 0)
    cube_magic = is_magic(power(A, 3))
    fifth_magic = is_magic(power(A, 5))
    if not (singular == cube_magic == fifth_magic):
        agree = False
for A in standard:
    if determinant(A) == 0:
        singular_standard += 1
singular_all = sum(1 for A in everything if determinant(A) == 0)

print("      singular squares: " + str(singular_standard) + " of 880 "
      "(" + str(singular_all) + " of 7040 counting rotations and reflections)")
check("640 of the 880 classical squares are singular", singular_standard == 640)
check("for all 7040 squares: singular <=> the cube is magic <=> the fifth power is magic",
      agree)
check("so exactly 240 of the 880 do not have a magic cube", 880 - singular_standard == 240)

# ----------------------------------------------------------------------
print("\nThe example D in the paper")
D = [[1, 6, 15, 12],
     [11, 16, 5, 2],
     [8, 3, 10, 13],
     [14, 9, 4, 7]]
show(D)
J = all_ones(4)
char = characteristic_polynomial(D)
print("      det(xI - D) = " + polynomial_to_text(char))
check("D is a classical magic square", tuple(x for row in D for x in row) in distinct)
check("D is neither associative nor pandiagonal",
      not is_associative(D) and not is_pandiagonal(D))
check("the eigenvalues of D are 34, 0, 8, -8, i.e. det(xI - D) = x(x - 34)(x^2 - 64)",
      char == [1, -34, -64, 2176, 0])
check("D^3 = 64 D + 9282 J", equal(power(D, 3), add(scale(64, D), scale(9282, J))))

summary()
