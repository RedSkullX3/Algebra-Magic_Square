"""
magic_tools.py

Small helper functions shared by all the other programs in this folder.

All arithmetic is exact. We only use whole numbers and fractions (Python's
built-in Fraction type), so there is never any rounding error. Only Python's
standard library is used, so nothing needs to be installed.

A matrix is stored as a list of rows. For example, the Lo Shu square is

    [[2, 9, 4],
     [7, 5, 3],
     [6, 1, 8]]
"""

from fractions import Fraction
from math import lcm


# ----------------------------------------------------------------------
# Building matrices
# ----------------------------------------------------------------------

def zero_matrix(n):
    """The n x n matrix of zeros."""
    return [[0] * n for _ in range(n)]


def identity(n):
    """The n x n identity matrix I."""
    I = zero_matrix(n)
    for i in range(n):
        I[i][i] = 1
    return I


def all_ones(n):
    """The n x n matrix J in which every entry is 1."""
    return [[1] * n for _ in range(n)]


def reversal(n):
    """The reversal matrix P, with ones on the anti-diagonal."""
    P = zero_matrix(n)
    for i in range(n):
        P[i][n - 1 - i] = 1
    return P


# ----------------------------------------------------------------------
# Arithmetic with matrices
# ----------------------------------------------------------------------

def add(A, B):
    n = len(A)
    return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]


def subtract(A, B):
    n = len(A)
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def scale(c, A):
    """The number c times the matrix A."""
    return [[c * x for x in row] for row in A]


def combine(terms):
    """A sum c1*A1 + c2*A2 + ... given as a list of pairs (c, A)."""
    n = len(terms[0][1])
    total = zero_matrix(n)
    for c, A in terms:
        total = add(total, scale(c, A))
    return total


def multiply(A, B):
    """The matrix product AB."""
    n = len(A)
    C = zero_matrix(n)
    for i in range(n):
        for k in range(n):
            a = A[i][k]
            if a != 0:
                for j in range(n):
                    C[i][j] += a * B[k][j]
    return C


def multiply_all(matrices):
    """The product of a list of matrices, in order."""
    result = matrices[0]
    for M in matrices[1:]:
        result = multiply(result, M)
    return result


def power(A, k):
    """A to the power k (k = 0, 1, 2, ...)."""
    result = identity(len(A))
    for _ in range(k):
        result = multiply(result, A)
    return result


def equal(A, B):
    n = len(A)
    return all(A[i][j] == B[i][j] for i in range(n) for j in range(n))


def trace(A):
    """Sum of the main diagonal."""
    return sum(A[i][i] for i in range(len(A)))


def anti_trace(A):
    """Sum of the anti-diagonal (top-right to bottom-left)."""
    n = len(A)
    return sum(A[i][n - 1 - i] for i in range(n))


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


# ----------------------------------------------------------------------
# Magic squares
# ----------------------------------------------------------------------

def line_sums(A):
    """All 2n+2 line sums: the rows, the columns, and the two diagonals."""
    n = len(A)
    rows = [sum(A[i]) for i in range(n)]
    columns = [sum(A[i][j] for i in range(n)) for j in range(n)]
    return rows + columns + [trace(A), anti_trace(A)]


def is_magic(A):
    """True if all rows, columns and both diagonals have the same sum."""
    sums = line_sums(A)
    return all(s == sums[0] for s in sums)


def is_semi_magic(A):
    """True if all rows and columns have the same sum (diagonals ignored)."""
    n = len(A)
    sums = line_sums(A)[:2 * n]
    return all(s == sums[0] for s in sums)


def is_zero_magic(A):
    """True if every row, column and diagonal sums to 0."""
    return all(s == 0 for s in line_sums(A))


def line_sum(A):
    return sum(A[0])


def is_associative(A):
    """True if every two cells opposite each other through the centre
    add up to the same number."""
    n = len(A)
    target = A[0][0] + A[n - 1][n - 1]
    return all(A[i][j] + A[n - 1 - i][n - 1 - j] == target
               for i in range(n) for j in range(n))


def is_pandiagonal(A):
    """True if all the 'broken' diagonals also have the magic sum."""
    n = len(A)
    s = line_sum(A)
    for shift in range(n):
        down = sum(A[i][(i + shift) % n] for i in range(n))
        up = sum(A[i][(shift - i) % n] for i in range(n))
        if down != s or up != s:
            return False
    return True


# ----------------------------------------------------------------------
# Determinant, inverse and characteristic polynomial (exact)
# ----------------------------------------------------------------------

def to_fractions(A):
    return [[Fraction(x) for x in row] for row in A]


def determinant(A):
    """The determinant, computed exactly by row reduction."""
    n = len(A)
    M = to_fractions(A)
    det = Fraction(1)
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
            det = -det
        det *= M[col][col]
        for r in range(col + 1, n):
            factor = M[r][col] / M[col][col]
            if factor != 0:
                for c in range(col, n):
                    M[r][c] -= factor * M[col][c]
    return det


def inverse(A):
    """The inverse matrix, computed exactly (A must be invertible)."""
    n = len(A)
    M = to_fractions(A)
    I = to_fractions(identity(n))
    for col in range(n):
        pivot = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[pivot] = M[pivot], M[col]
        I[col], I[pivot] = I[pivot], I[col]
        p = M[col][col]
        M[col] = [x / p for x in M[col]]
        I[col] = [x / p for x in I[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [M[r][j] - f * M[col][j] for j in range(n)]
                I[r] = [I[r][j] - f * I[col][j] for j in range(n)]
    return I


def characteristic_polynomial(A):
    """Coefficients of det(xI - A), highest power of x first.

    For example [1, -15, 24, -360] means x^3 - 15x^2 + 24x - 360.
    This uses the Faddeev-LeVerrier method, which needs only matrix
    products and traces, so it is exact.
    """
    n = len(A)
    coefficients = [Fraction(1)]
    M = zero_matrix(n)
    for k in range(1, n + 1):
        M = add(multiply(A, M), scale(coefficients[-1], identity(n)))
        coefficients.append(-Fraction(trace(multiply(A, M))) / k)
    return coefficients


def multiply_polynomials(p, q):
    """Multiply two polynomials given as coefficient lists (highest first)."""
    result = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i + j] += a * b
    return result


def polynomial_to_text(coefficients):
    """Write a coefficient list as readable text, e.g. 'x^2 - 3x + 2'."""
    degree = len(coefficients) - 1
    parts = []
    for i, c in enumerate(coefficients):
        c = Fraction(c)
        if c == 0:
            continue
        d = degree - i
        sign = "-" if c < 0 else "+"
        size = abs(c)
        number = str(size.numerator) if size.denominator == 1 else str(size)
        if d == 0:
            term = number
        else:
            power_text = "x" if d == 1 else "x^" + str(d)
            term = power_text if size == 1 else number + power_text
        parts.append((sign, term))
    if not parts:
        return "0"
    text = ("-" if parts[0][0] == "-" else "") + parts[0][1]
    for sign, term in parts[1:]:
        text += " " + sign + " " + term
    return text


def minimal_polynomial_degree(A):
    """The smallest d such that I, A, A^2, ..., A^d are linearly dependent."""
    n = len(A)
    powers = [identity(n)]
    for d in range(1, n + 1):
        powers.append(multiply(powers[-1], A))
        vectors = [[x for row in P for x in row] for P in powers]
        if rank(vectors) < len(vectors):
            return d
    return n


def rank(vectors):
    """Rank of a list of vectors (exact row reduction)."""
    rows = [[Fraction(x) for x in v] for v in vectors]
    r = 0
    columns = len(rows[0]) if rows else 0
    for c in range(columns):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [rows[i][j] - f * rows[r][j] for j in range(columns)]
        r += 1
        if r == len(rows):
            break
    return r


# ----------------------------------------------------------------------
# Bases of zero-magic squares
# ----------------------------------------------------------------------

def solve_homogeneous(equations, number_of_unknowns):
    """A basis of all solutions of a system 'equations = 0', as lists of
    whole numbers. Each equation is a list of coefficients."""
    rows = [[Fraction(x) for x in eq] for eq in equations]
    pivot_columns = []
    r = 0
    for c in range(number_of_unknowns):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        p = rows[r][c]
        rows[r] = [x / p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [rows[i][j] - f * rows[r][j]
                           for j in range(number_of_unknowns)]
        pivot_columns.append(c)
        r += 1
        if r == len(rows):
            break
    free_columns = [c for c in range(number_of_unknowns) if c not in pivot_columns]
    basis = []
    for free in free_columns:
        v = [Fraction(0)] * number_of_unknowns
        v[free] = Fraction(1)
        for i, pc in enumerate(pivot_columns):
            v[pc] = -rows[i][free]
        common = lcm(*[x.denominator for x in v])
        basis.append([int(x * common) for x in v])
    return basis


def line_sum_equations(n):
    """The 2n+2 equations 'line sum = 0' for an n x n matrix."""
    equations = []
    for i in range(n):
        equations.append([1 if k // n == i else 0 for k in range(n * n)])
        equations.append([1 if k % n == i else 0 for k in range(n * n)])
    equations.append([1 if k // n == k % n else 0 for k in range(n * n)])
    equations.append([1 if k // n + k % n == n - 1 else 0 for k in range(n * n)])
    return equations


def vector_to_matrix(v, n):
    return [v[i * n:(i + 1) * n] for i in range(n)]


def zero_magic_basis(n):
    """A basis (with whole-number entries) of all n x n zero-magic squares."""
    vectors = solve_homogeneous(line_sum_equations(n), n * n)
    return [vector_to_matrix(v, n) for v in vectors]


def associative_zero_magic_basis(n):
    """A basis of the n x n zero-magic squares in which opposite cells
    through the centre are negatives of each other."""
    equations = line_sum_equations(n)
    for i in range(n):
        for j in range(n):
            eq = [0] * (n * n)
            eq[i * n + j] += 1
            eq[(n - 1 - i) * n + (n - 1 - j)] += 1
            equations.append(eq)
    vectors = solve_homogeneous(equations, n * n)
    return [vector_to_matrix(v, n) for v in vectors]


def place_in_rows_and_columns(small, positions, n):
    """Put the entries of a small square matrix into the given rows and
    columns (numbered from 0) of an n x n zero matrix."""
    big = zero_matrix(n)
    for a, r in enumerate(positions):
        for b, c in enumerate(positions):
            big[r][c] = small[a][b]
    return big


# ----------------------------------------------------------------------
# Printing results
# ----------------------------------------------------------------------

passed = 0
failed = 0


def check(description, condition):
    """Print one result line and keep count."""
    global passed, failed
    if condition:
        passed += 1
        print("  PASS  " + description)
    else:
        failed += 1
        print("  FAIL  " + description)


def summary():
    print()
    print("Result: " + str(passed) + " checks passed, " + str(failed) + " failed.")
    if failed == 0:
        print("Everything agrees with the paper.")


def show(A, indent="      "):
    """Print a matrix neatly."""
    width = max(len(str(x)) for row in A for x in row)
    for row in A:
        print(indent + "  ".join(str(x).rjust(width) for x in row))
