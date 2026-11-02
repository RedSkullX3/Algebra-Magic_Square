# Code for "The Hidden Algebra of Magic Squares"

These short programs check, by computer, every formula, example and number
in the paper. The check programs are written in plain Python and use only
what comes with Python, so there is nothing to install. (Two extra
exploration tools, described at the end, need a few extra packages.)

All the arithmetic is **exact**. The programs work with whole numbers and
fractions, never with rounded decimals, so a result such as "the trace is 2448"
is exactly 2448. The one exception is the last check in `check_odd_powers.py`,
which solves an equation numerically using 80 decimal digits.

## How to run

You need Python 3.9 or newer. Open a terminal in this folder and type

```
python3 run_all.py
```

to run everything (it takes a few seconds), or run one program at a time,
for example

```
python3 classical_4x4.py
```

Each program prints a list of checks. `PASS` means the computer confirms the
statement in the paper. At the end it prints how many checks passed. If
anything were wrong, the line would say `FAIL` instead.

## What each file is for

| File | What it does | Where in the paper |
|---|---|---|
| `magic_tools.py` | A toolbox used by all the other programs: adding and multiplying matrices, testing whether a square is magic, determinants, characteristic polynomials, and so on. It does not print anything by itself. | (used everywhere) |
| `check_3x3.py` | Checks everything about 3 x 3 magic squares: the canonical form and the Lo Shu square, the family of bases and the table, the "2 x 2 picture", the characteristic polynomial and determinant, why products of two squares fail, cubes and odd products, the inverse formulas, the minimal polynomial, commuting squares, and the orderings of a triple product. Each formula is tested on hundreds of randomly chosen squares. | Theorems 1 to 10 |
| `check_counterexamples.py` | Checks that odd powers stop being magic in larger squares: Hill and Elzaidi's 4 x 4 example, the square A4 (its cube has trace 2448 and anti-diagonal sum -816), and the same example placed inside squares of every order from 5 to 12. | Question 11 and the section "Why order three is special" |
| `classical_4x4.py` | Finds **all** 4 x 4 magic squares that use the numbers 1 to 16 (there are 880, or 7040 counting rotations and reflections). It then checks that exactly 640 of the 880 are singular, and that a square has a magic cube exactly when it is singular. It also checks the example D. | Corollary 13 and the 640 of 880 result |
| `check_odd_powers.py` | Checks the last section: the formula for the odd powers of a singular 4 x 4 square, the examples of Freitag, Hill and Elzaidi, and Cook, Bacon and Hillman, random associative squares, the 5 x 5 square C and its copies in every odd order up to 11, and the high-precision examples in orders 6 and 7 whose cube is magic but whose fifth power is not. | Theorem 12, Corollary 13 and Remark 14 |
| `run_all.py` | Runs the four check programs above, one after another. | |

## A few words on the methods

* **Testing a formula.** For identities such as A^(2k+1) = (xy)^k A0 + 3^(2k) c^(2k+1) J,
  the programs pick hundreds of random magic squares and check that both
  sides are exactly equal every time.
* **Finding all 4 x 4 squares.** The program fills in the two diagonals
  first. The rows and columns must add up to 34, and that then fixes most of
  the other cells, so the search is quick. It only looks for squares whose
  smallest corner is in the top left, and then adds the other seven
  rotations and reflections of each one.
* **Characteristic polynomials.** These are computed with the
  Faddeev-LeVerrier method, which needs only matrix products and traces, so
  it stays exact.

## Two extra exploration tools

These two files are for trying things out by hand. They are not needed to
check the paper, and they use the `numpy` package (and the notebook also uses
`prettytable`). Install them with `pip install numpy prettytable` if needed.

### `hypothesis.py`

An interactive tool. Run it with `python3 hypothesis.py` and choose:

1. **Generate a magic square.** You give the size n (any n except 2), a
   starting value and a step. It builds an n x n magic square whose entries
   are start, start + step, start + 2 x step, and so on, using the classical
   construction methods (the Siamese method for odd n, the "X pattern" for n
   divisible by 4, and Strachey's method for n = 6, 10, 14, ...).
2. **Enter your own matrix.** You type in any square matrix row by row.

It then prints the matrix, computes its cube A^3, and says whether A^3 is
magic. If A^3 is not magic but its rows and columns still agree, it points
out that A^3 is semi-magic, with only the diagonals failing.

Good to know: the squares this tool generates all turn out to have magic
cubes (we checked n = 3 to 10). The odd-order squares and those with n
divisible by 4 are associative, and the Strachey squares of orders 6 and 10
also have magic cubes. To see a cube that is **not** magic, use option 2
with, for example, the square A4 or Freitag's square from the paper.
`numpy` uses fixed-size whole numbers, so for very large squares (hundreds of
rows) the entries of the cube can become too big and overflow.

### `Magic Sqaure Matrix Github.ipynb`

A Jupyter notebook with small experiments on 3 x 3 magic squares. Open it
with Jupyter (for example `jupyter notebook`). Its cells do the following,
in order:

1. Take a sample 3 x 3 magic square, compute its cube, and check that the
   cube is magic.
2. Decide whether a matrix is magic, semi-magic with equal diagonals, or
   semi-magic with unequal diagonals. The example is the product of two 3 x 3
   magic squares written with letters, which comes out only semi-magic.
3. Build the special family of 3 x 3 squares
   [[n-d, n-2d, n+3d], [n+4d, n, n-4d], [n-3d, n+2d, n+d]]
   from n and d (a family used in an earlier draft of the paper).
4. Make a random 3 x 3 magic square with a magic sum you choose, using the
   formula of equation (4) in the paper with b = 2a.
5. The same with whole-number entries (the magic sum must be a multiple of 3).
6. The same with real-number entries.
7. The same with a = 1 and b = 2, which gives a square of consecutive
   numbers (for magic sum 15 it gives a turned version of the Lo Shu square).
8. Read in nine whole numbers and decide whether they can be arranged into a
   3 x 3 magic square, printing the arrangement if they can.

