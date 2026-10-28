#!/usr/bin/env python3
"""Summarise results/hunts/*.log into paper/generated/c1_macros.tex."""

import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def total(prefix):
    tests = fails = inst = 0
    for path in glob.glob(os.path.join(ROOT, "results", "hunts", prefix + "_*.log")):
        found = re.findall(r"instances (\d+) tests (\d+) failures (\d+)", open(path).read())
        if found:
            m = re.match(r"(\d+) (\d+) (\d+)", " ".join(found[-1]))
            inst += int(m.group(1)); tests += int(m.group(2)); fails += int(m.group(3))
    return inst, tests, fails


def fmt(x):
    return f"{x:,}".replace(",", "{,}")


def main():
    vi, vt, vf = total("v5")
    ci, ct, cf = total("c1")
    fi, ft, ff = total("v4")
    pi, pt, pf = total("v4p")
    qi, qt, qf = total("c1p")
    m1 = total("v4pm"); m2 = total("v4pmp")
    mi, mt, mf = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
    os.makedirs(os.path.join(ROOT, "paper", "generated"), exist_ok=True)
    with open(os.path.join(ROOT, "paper", "generated", "c1_macros.tex"), "w") as fh:
        fh.write(f"\\newcommand{{\\VFiveInst}}{{{fmt(vi)}}}\n")
        fh.write(f"\\newcommand{{\\VFiveTests}}{{{fmt(vt)}}}\n")
        fh.write(f"\\newcommand{{\\VFiveFails}}{{{vf}}}\n")
        fh.write(f"\\newcommand{{\\ConeInst}}{{{fmt(ci)}}}\n")
        fh.write(f"\\newcommand{{\\ConeTests}}{{{fmt(ct)}}}\n")
        fh.write(f"\\newcommand{{\\ConeFails}}{{{cf}}}\n")
        for name, (a, b, c) in {"VFour": (fi, ft, ff), "VFourPert": (pi, pt, pf),
                                "ConePert": (qi, qt, qf),
                                "VFourPM": (mi, mt, mf)}.items():
            fh.write(f"\\newcommand{{\\{name}Inst}}{{{fmt(a)}}}\n")
            fh.write(f"\\newcommand{{\\{name}Tests}}{{{fmt(b)}}}\n")
            fh.write(f"\\newcommand{{\\{name}Fails}}{{{c}}}\n")
        fh.write(f"\\newcommand{{\\ConeAllTests}}{{{fmt(ct + qt)}}}\n")
    print(f"V5 {vi}/{vt}/{vf}  C1 {ci}/{ct}/{cf}  V4 {fi}/{ft}/{ff}  "
          f"V4-perturbed {pi}/{pt}/{pf}  C1-perturbed {qi}/{qt}/{qf}")


if __name__ == "__main__":
    main()
