"""
check_counterexamples.py

For 3 x 3 magic squares, every odd power is magic. This program checks the
counterexamples in the paper that show this fails in every larger order
n = 4, 5, 6, ...

Run with:   python3 check_counterexamples.py
"""

from magic_tools import (identity, power, equal, trace, anti_trace,
                         line_sums, is_magic, is_zero_magic, determinant,
                         place_in_rows_and_columns, check, summary, show)

print("Counterexamples in larger orders")
print("================================")

# ----------------------------------------------------------------------
print("\nHill and Elzaidi's 4 x 4 permutation matrix")
H = [[1, 0, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [0, 1, 0, 0]]
check("H is magic (line-sum 1)", is_magic(H))
check("H^3 is the identity matrix", equal(power(H, 3), identity(4)))
check("the identity matrix is not magic", not is_magic(identity(4)))
check("det(H) = 1", determinant(H) == 1)

# ----------------------------------------------------------------------
print("\nThe zero-magic square A4")
A4 = [[1, 2, 3, -6],
      [4, 5, -1, -8],
      [0, 12, -16, 4],
      [-5, -19, 14, 10]]
show(A4)
cube = power(A4, 3)
check("A4 is zero-magic (every line sums to 0)", is_zero_magic(A4))
check("the rows and columns of A4^3 sum to 0", all(s == 0 for s in line_sums(cube)[:8]))
check("trace of A4^3 is 2448", trace(cube) == 2448)
check("anti-diagonal sum of A4^3 is -816", anti_trace(cube) == -816)

# ----------------------------------------------------------------------
print("\nOrders n = 5, 6, ..., 12: A4 placed in rows and columns 1, 2, n-1, n")
for n in range(5, 13):
    positions = [0, 1, n - 2, n - 1]          # rows 1, 2, n-1, n (counted from 0)
    An = place_in_rows_and_columns(A4, positions, n)
    cube = power(An, 3)
    check("n = " + str(n) + ": A_n is zero-magic and the trace of A_n^3 is 2448",
          is_zero_magic(An) and trace(cube) == 2448)

summary()
