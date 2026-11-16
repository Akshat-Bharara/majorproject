#!/usr/bin/env python3
"""Run every experiment of the paper and regenerate its tables and plots.

    python3 run_convex_experiments.py [--quick]

Steps: build, self-test, counterexample check, correctness splits (with the
Python and FMM cross-checks), property hunts, timing splits, tables.
--quick skips the hunts and the instances with n >= 300.
Timing runs are sequential on one core; do not run other jobs alongside.
"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
quick = "--quick" in sys.argv


def sh(cmd, out=None):
    print("+", cmd if len(cmd) < 150 else cmd[:150] + " ...", flush=True)
    with (open(out, "w") if out else open(os.devnull, "w")) as f:
        subprocess.run(cmd, shell=True, check=True, cwd=HERE, stdout=f if out else None)


def files(split, pred=lambda f: True):
    return " ".join(sorted(f for f in glob.glob(os.path.join("data", split, "*.txt"), root_dir=HERE) if pred(f)))


def main():
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    sh("make -s all")
    if not os.path.isdir(os.path.join(HERE, "data", "verify")):
        sh("python3 gen_convex_dataset.py data")
    sh("bin/selftest 3000 1")
    sh("python3 experiments/check_counterexamples.py")
    sh(f"bin/bench --reps 1 {files('verify')}", "results/verify.csv")
    sh(f"bin/bench --reps 1 --brute-max 0 --brute-kmax 5 {files('medium')}", "results/medium.csv")
    sh("python3 experiments/verify_python.py results/verify.csv bin/fmm", "results/verify_python.txt")
    if not quick:
        sh("bin/hunt exh 700000 101 12 6", "results/hunt_exh.txt")
        sh("bin/hunt big 70000 202 30 12", "results/hunt_big.txt")
        sh("bin/hunt pert 800000 303 data/counterexample_bases.txt", "results/hunt_pert.txt")
    big = lambda f: "_n300_" in f or "_n400_" in f
    sh(f"bin/bench --reps 3 --brute-max 0 --brute-kmax 0 {files('scaling', lambda f: not big(f))}", "results/scaling_small.csv")
    if not quick:
        sh(f"bin/bench --reps 1 --brute-max 0 --brute-kmax 0 {files('scaling', big)}", "results/scaling_large.csv")
    sh(f"bin/bench --reps 3 --brute-max 0 --brute-kmax 0 {files('ksweep')}", "results/ksweep.csv")
    sh("python3 make_tables.py")


if __name__ == "__main__":
    main()
