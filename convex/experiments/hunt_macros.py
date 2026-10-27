#!/usr/bin/env python3
"""Summarise results/hunts/*.log into paper/generated/c1_macros.tex."""

import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def total(prefix):
    tests = fails = inst = 0
    for path in glob.glob(os.path.join(ROOT, "results", "hunts", prefix + "_*.log")):
        m = re.search(r"instances (\d+) tests (\d+) failures (\d+)", open(path).read())
        if m:
            inst += int(m.group(1)); tests += int(m.group(2)); fails += int(m.group(3))
    return inst, tests, fails


def fmt(x):
    return f"{x:,}".replace(",", "{,}")


def main():
    vi, vt, vf = total("v5")
    ci, ct, cf = total("c1")
    os.makedirs(os.path.join(ROOT, "paper", "generated"), exist_ok=True)
    with open(os.path.join(ROOT, "paper", "generated", "c1_macros.tex"), "w") as fh:
        fh.write(f"\\newcommand{{\\VFiveInst}}{{{fmt(vi)}}}\n")
        fh.write(f"\\newcommand{{\\VFiveTests}}{{{fmt(vt)}}}\n")
        fh.write(f"\\newcommand{{\\VFiveFails}}{{{vf}}}\n")
        fh.write(f"\\newcommand{{\\ConeInst}}{{{fmt(ci)}}}\n")
        fh.write(f"\\newcommand{{\\ConeTests}}{{{fmt(ct)}}}\n")
        fh.write(f"\\newcommand{{\\ConeFails}}{{{cf}}}\n")
    print(f"V5: {vi} instances, {vt} tests, {vf} failures; C1: {ci} instances, {ct} tests, {cf} failures")


if __name__ == "__main__":
    main()
