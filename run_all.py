"""
run_all.py

Runs all four check programs, one after another.

Run with:   python3 run_all.py
"""

import os
import subprocess
import sys

folder = os.path.dirname(os.path.abspath(__file__))
programs = ["check_3x3.py",
            "check_counterexamples.py",
            "classical_4x4.py",
            "check_odd_powers.py"]

all_good = True
for name in programs:
    print()
    print("#" * 70)
    print("# " + name)
    print("#" * 70)
    result = subprocess.run([sys.executable, name], cwd=folder,
                            capture_output=True, text=True)
    print(result.stdout, end="")
    if result.returncode != 0 or "FAIL" in result.stdout:
        all_good = False
        print(result.stderr, end="")

print()
if all_good:
    print("All programs finished and every check passed.")
else:
    print("Some check failed. Look for FAIL above.")
