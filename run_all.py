"""
run_all.py

Run every builder in dependency order. Exits non-zero on the first failure.

Usage:
  python3 run_all.py
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

# steps that write a table a later builder reads
PREREQUISITES = [
    "figures/main_05/scripts/module_group_permutation.py",
    "figures/main_06/scripts/module_profile_prepare.py",
    "figures/main_06/scripts/module_permutation_bootstrap.py",
    "figures/suppl_06/scripts/module_coherence_permutation.py",
    "figures/suppl_06/scripts/module_coherence_bh.py",
    "figures/suppl_07/scripts/cortex_group_permutation.py",
]


def builders():
    found = []
    for base in ("figures", "tables"):
        for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, base)):
            dirnames[:] = [d for d in dirnames if d != "__pycache__"]
            if os.path.basename(dirpath) != "scripts":
                continue
            for name in sorted(filenames):
                if name.endswith(".py"):
                    found.append(os.path.relpath(os.path.join(dirpath, name), ROOT))
    return found


def main():
    order = PREREQUISITES + [b for b in builders() if b not in PREREQUISITES]
    failed = []
    for rel in order:
        path = os.path.join(ROOT, rel)
        if not os.path.isfile(path):
            continue
        r = subprocess.run([sys.executable, path], capture_output=True, text=True)
        print("%-4s %s" % ("OK" if r.returncode == 0 else "FAIL", rel))
        if r.returncode != 0:
            print(r.stdout[-2000:])
            print(r.stderr[-2000:])
            failed.append(rel)
    print("\n%d builders, %d failed" % (len(order), len(failed)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
