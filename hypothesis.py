"""
Magic Square Generator and A^3 Magic-Property Checker
======================================================
Generates an n x n magic square using a chosen starting value and step,
displays the matrix, computes A^3 (matrix cube), and checks whether
A^3 is also a magic square (equal row, column, and both diagonal sums).
"""

import numpy as np


# ---------- Magic square construction ----------

def _odd_magic(n):
    """Siamese (De La Loubere) method for odd n. Fills with 1..n^2."""
    sq = np.zeros((n, n), dtype=int)
    i, j = 0, n // 2
    for k in range(1, n * n + 1):
        sq[i, j] = k
        ni, nj = (i - 1) % n, (j + 1) % n
        if sq[ni, nj]:
            i = (i + 1) % n
        else:
            i, j = ni, nj
    return sq


def _doubly_even_magic(n):
    """X-pattern method for n divisible by 4. Fills with 1..n^2."""
    sq = np.arange(1, n * n + 1).reshape(n, n)
    for i in range(n):
        for j in range(n):
            if (i % 4 == j % 4) or ((i % 4) + (j % 4) == 3):
                sq[i, j] = n * n + 1 - sq[i, j]
    return sq


def _singly_even_magic(n):
    """Strachey method for n = 4k + 2 (n >= 6). Fills with 1..n^2."""
    m = n // 2          # m is odd
    k = (n - 2) // 4

    A = _odd_magic(m)              # quadrant values: 1..m^2
    B = A + m * m                  # m^2+1..2m^2
    C = A + 2 * m * m              # 2m^2+1..3m^2
    D = A + 3 * m * m              # 3m^2+1..4m^2

    sq = np.zeros((n, n), dtype=int)
    sq[:m, :m] = A      # top-left
    sq[:m, m:] = C      # top-right
    sq[m:, :m] = D      # bottom-left
    sq[m:, m:] = B      # bottom-right

    mid = m // 2

    # Left swaps between A (rows 0..m-1) and D (rows m..2m-1)
    for i in range(m):
        for c in range(k):
            j = c + 1 if i == mid else c     # middle row: shift swap by one
            sq[i, j], sq[i + m, j] = sq[i + m, j], sq[i, j]

    # Right swaps between C (top right) and B (bottom right)
    for i in range(m):
        for c in range(k - 1):
            j = n - 1 - c
            sq[i, j], sq[i + m, j] = sq[i + m, j], sq[i, j]

    return sq


def generate_magic_square(n, start=1, step=1):
    """
    Build an n x n magic square whose entries are an arithmetic
    progression: start, start+step, start+2*step, ..., start+(n^2-1)*step.

    Different (start, step) values give different magic squares with
    different magic sums:
        magic_sum = n * start + step * n * (n^2 - 1) / 2

    Args:
        n     : size of the square (n != 2)
        start : starting value of the arithmetic progression
        step  : common difference of the arithmetic progression

    Returns:
        numpy.ndarray of shape (n, n)
    """
    if n < 1:
        raise ValueError("n must be at least 1.")
    if n == 2:
        raise ValueError("No magic square of order 2 exists.")

    if n % 2 == 1:
        base = _odd_magic(n)
    elif n % 4 == 0:
        base = _doubly_even_magic(n)
    else:
        base = _singly_even_magic(n)

    return start + (base - 1) * step


# ---------- Magic-square verification ----------

def is_magic_square(M):
    """
    Check whether square matrix M is a magic square: every row, every
    column, and both diagonals share the same sum.

    Returns: (is_magic, magic_sum_or_None, diagnostics_dict)
    """
    M = np.asarray(M)
    n = M.shape[0]
    if M.shape != (n, n):
        return False, None, {"error": "not square"}

    row_sums  = M.sum(axis=1)
    col_sums  = M.sum(axis=0)
    diag_main = int(np.trace(M))
    diag_anti = int(np.trace(np.fliplr(M)))

    target = int(row_sums[0])
    diagnostics = {
        "row_sums":  row_sums.tolist(),
        "col_sums":  col_sums.tolist(),
        "diag_main": diag_main,
        "diag_anti": diag_anti,
    }

    ok = (np.all(row_sums == target) and
          np.all(col_sums == target) and
          diag_main == target and
          diag_anti == target)
    return (ok, target if ok else None, diagnostics)


# ---------- Pretty printing ----------

def print_matrix(M, label=""):
    if label:
        print(label)
    width = len(str(int(np.max(np.abs(M))))) + 2
    for row in M:
        print("".join(f"{int(v):>{width}d}" for v in row))


# ---------- Custom matrix input ----------

def enter_custom_matrix():
    """Prompt the user to enter an n x n matrix row by row."""
    while True:
        try:
            n = int(input("Enter matrix size n (n x n): "))
            if n < 1:
                print("  n must be at least 1.")
                continue
            break
        except ValueError:
            print("  Please enter a valid integer.")

    print(f"Enter the matrix row by row ({n} values per row, space-separated):")
    rows = []
    for i in range(n):
        while True:
            try:
                raw = input(f"  Row {i + 1}: ").strip().split()
                if len(raw) != n:
                    print(f"  Expected {n} values, got {len(raw)}. Try again.")
                    continue
                rows.append([int(x) for x in raw])
                break
            except ValueError:
                print("  Please enter integers only.")

    return np.array(rows, dtype=int)


def check_and_cube(A, label="A"):
    """Print the matrix, check if it is a magic square, then compute and check A^3."""
    print_matrix(A, label=f"\nMatrix {label}:\n")

    ok, s, diag = is_magic_square(A)
    print(f"\nIs {label} a magic square?  {ok}")
    if ok:
        print(f"Magic sum of {label}:        {s}")
    else:
        print(f"  Row sums:       {diag['row_sums']}")
        print(f"  Column sums:    {diag['col_sums']}")
        print(f"  Main diagonal:  {diag['diag_main']}")
        print(f"  Anti-diagonal:  {diag['diag_anti']}")

    A3 = A @ A @ A
    print_matrix(A3, label=f"\n{label} cubed  ({label}^3 = {label} . {label} . {label}):\n")

    ok3, s3, diag3 = is_magic_square(A3)
    print("\n" + "-" * 64)
    print(f"Is {label}^3 a magic square?")
    print("-" * 64)
    if ok3:
        print(f"YES — every row, column, and diagonal of {label}^3 sums to {s3}.")
    else:
        print(f"NO — the row, column, and diagonal sums are not all equal.")
        print(f"  Row sums:       {diag3['row_sums']}")
        print(f"  Column sums:    {diag3['col_sums']}")
        print(f"  Main diagonal:  {diag3['diag_main']}")
        print(f"  Anti-diagonal:  {diag3['diag_anti']}")
        rows_equal = len(set(diag3['row_sums'])) == 1
        cols_equal = len(set(diag3['col_sums'])) == 1
        if rows_equal and cols_equal and ok:
            r = diag3['row_sums'][0]
            print(f"\n  Note: {label}^3 IS 'semi-magic' — every row and column sums to {r}.")
            print(f"        (In fact this equals (magic_sum_of_{label})^3 = {s}^3 = {s**3}.)")
            print(f"        Only the two diagonals fail to match.")


# ---------- Main demo ----------

def main():
    print("=" * 64)
    print("  Magic Square Tool")
    print("=" * 64)
    print("  1. Generate a magic square")
    print("  2. Enter a custom matrix")
    print("=" * 64)

    while True:
        choice = input("Choose an option (1 or 2): ").strip()
        if choice in ("1", "2"):
            break
        print("  Please enter 1 or 2.")

    if choice == "2":
        A = enter_custom_matrix()
        check_and_cube(A, label="A")
        return

    # --- option 1: generate ---
    while True:
        try:
            n = int(input("Enter size n (n != 2, n >= 1): "))
            if n == 2:
                print("  No magic square of order 2 exists. Try another value.")
                continue
            if n < 1:
                print("  n must be at least 1.")
                continue
            break
        except ValueError:
            print("  Please enter a valid integer.")

    while True:
        try:
            start = int(input("Enter starting value (default 1): ") or "1")
            break
        except ValueError:
            print("  Please enter a valid integer.")

    while True:
        try:
            step = int(input("Enter step/common difference (default 1): ") or "1")
            break
        except ValueError:
            print("  Please enter a valid integer.")

    print("=" * 64)
    print(f"Generating magic square:  n={n},  start={start},  step={step}")
    print("=" * 64)

    A = generate_magic_square(n, start=start, step=step)

    expected_sum = n * (2 * start + step * (n * n - 1)) // 2
    print_matrix(A, label=f"\nMagic square A ({n} x {n}):\n")

    ok, s, _ = is_magic_square(A)
    print(f"\nIs A a magic square?    {ok}")
    print(f"Magic sum of A:         {s}   (formula gives {expected_sum})")

    # ---- compute A^3 (matrix cube, not element-wise) ----
    A3 = A @ A @ A
    print_matrix(A3, label=f"\nA cubed   (A^3 = A . A . A):\n")

    # ---- check whether A^3 is a magic square ----
    ok3, s3, diag = is_magic_square(A3)
    print("\n" + "-" * 64)
    print("Is A^3 a magic square?")
    print("-" * 64)
    if ok3:
        print(f"YES — every row, column, and diagonal of A^3 sums to {s3}.")
    else:
        print("NO — the row, column, and diagonal sums are not all equal.")
        print(f"  Row sums:       {diag['row_sums']}")
        print(f"  Column sums:    {diag['col_sums']}")
        print(f"  Main diagonal:  {diag['diag_main']}")
        print(f"  Anti-diagonal:  {diag['diag_anti']}")
        rows_equal = len(set(diag['row_sums'])) == 1
        cols_equal = len(set(diag['col_sums'])) == 1
        if rows_equal and cols_equal:
            r = diag['row_sums'][0]
            print(f"\n  Note: A^3 IS 'semi-magic' — every row and column sums "
                  f"to {r}.")
            print(f"        (In fact this equals (magic_sum_of_A)^3 = "
                  f"{s}^3 = {s**3}.)")
            print(f"        Only the two diagonals fail to match.")


if __name__ == "__main__":
    main()