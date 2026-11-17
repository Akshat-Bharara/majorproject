#!/usr/bin/env python3
"""results/ -> LaTeX tables and macros for Section 7 of the paper
(../paper/generated/arc_*.tex) and Markdown tables for NOTES.md (results/tables.md)."""
import collections
import csv
import math
import os
import re
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "..", "paper", "generated")
num = lambda v: f"{v:,}".replace(",", "{,}")


def fmt(x):
    if x is None:
        return "--"
    return f"{x:.0f}" if x >= 100 else f"{x:.1f}" if x >= 10 else f"{x:.2f}" if x >= 0.1 else f"{x:.3f}"


def main():
    os.makedirs(OUT, exist_ok=True)
    macros, md = [], []
    mac = lambda k, v: macros.append(f"\\newcommand{{\\{k}}}{{{v}}}")

    # ---- correctness ----
    txt = open(os.path.join(RES, "correctness.txt")).read()
    a, c = txt.split("PART C")
    small = [l.split() for l in a.splitlines()[2:] if l.strip()]
    big = [l.split() for l in c.splitlines()[2:] if l.strip() and not l.startswith("FAILURES")]
    fails = int(re.search(r"FAILURES (\d+)", txt).group(1))
    assert fails == 0, "correctness failures: no tables"
    mac("arcSmallInst", num(sum(int(r[1]) for r in small)))
    mac("arcDchecks", num(sum(int(r[5]) for r in small)))
    mac("arcOptChecks", num(sum(int(r[6]) for r in small)))
    mac("arcBigInst", num(sum(int(r[1]) for r in big)))
    mac("arcBigChecks", num(sum(int(r[3]) + int(r[4]) for r in big)))
    thin = [r for r in small if r[0] in ("ellipse", "twochain")]
    mac("arcUncheckedWrong", f"{100 * sum(int(r[10]) for r in thin) / sum(int(r[9]) for r in thin):.0f}")
    py = open(os.path.join(RES, "pybrute.txt")).read()
    mac("arcPyAtwo", num(int(re.search(r"A2 checks (\d+)", py).group(1))))
    assert "failures 0" in py
    md.append("### Correctness, n <= 16 (exhaustive)\n\n| family | instances | unimodal | radii tested | radii with arc property | D(r) vs brute force | optimum checks (all k) | adaptive answered | adaptive n/a | unchecked D tests | unchecked D wrong |\n|---|---|---|---|---|---|---|---|---|---|---|")
    md += ["| " + " | ".join(r) + " |" for r in small]
    md.append("\n### Correctness, 20 <= n <= 150 (capped Delaunay DP as oracle)\n\n| family | instances | unimodal | O1 and O2 checks | adaptive answered (= oracle) | adaptive n/a |\n|---|---|---|---|---|---|")
    md += ["| " + " | ".join(r) + " |" for r in big]
    md.append(f"\nFailures: {fails}. Python: `{py.strip()}`\n")

    # ---- hypotheses ----
    hyp = collections.defaultdict(dict)
    for line in open(os.path.join(RES, "probe_small_summary.txt")):
        f = line.split()
        hyp[f[0]][int(f[1][2:])] = (int(f[2]), int(f[3]))
    rows = []
    md.append("### Weakened hypotheses (n <= 13, exhaustive)\n\n| hypothesis | k=3 | k=4 | k=5 | k=6 |\n|---|---|---|---|---|")
    for h in ("H1", "H2", "H2u", "H3"):
        cells = [f"{num(hyp[h][k][0])} / {num(hyp[h][k][1])}" for k in (3, 4, 5, 6)]
        rows.append(f"{h} & " + " & ".join(cells) + " \\\\")
        md.append(f"| {h} | " + " | ".join(c.replace("{,}", ",") for c in cells) + " |")
    open(os.path.join(OUT, "arc_hyp.tex"), "w").write("\n".join(rows) + "\n")
    mac("arcHoneFourFive", num(hyp["H1"][4][0] + hyp["H1"][5][0]))
    mac("arcHtwo", num(sum(hyp["H2"][k][0] for k in (3, 4, 5, 6))))

    # ---- structure at r* ----
    rows = []
    md.append("\n### Structure at r* (20 <= n <= 60, k in {4, 8})\n\n| group | tests | mean fraction of arc-shaped points | mean beta | beta = 0 | greedy chain r*-far | adaptive answers | unchecked O1 = r* |\n|---|---|---|---|---|---|---|---|")
    lines = open(os.path.join(RES, "probe_stats.txt")).read().splitlines()[1:]
    cut = lines.index("--")
    order = ["beta=0", "beta 1-2", "beta 3-10", "beta 11-40", "beta >40"]
    lines = lines[:cut + 1] + sorted(lines[cut + 1:], key=lambda l: order.index(l.split("|")[0].strip()))
    for line in lines:
        if line.startswith("--"):
            rows.append("\\midrule")
            continue
        f = [x.strip() for x in line.split("|")]
        name = f[0].replace("beta=0", "$\\beta=0$").replace("beta ", "$\\beta$ ").replace(">", "$>$")
        if not name.startswith("$"):
            name = f"\\texttt{{{name}}}"
        rows.append(f"{name} & {f[1]} & {float(f[2]):.2f} & {float(f[3]):.1f} & {f[4]} & {f[5]} & {f[6]} \\\\")
        md.append("| " + " | ".join(f) + " |")
    open(os.path.join(OUT, "arc_stats.tex"), "w").write("\n".join(rows) + "\n")

    # ---- scaling ----
    data = collections.defaultdict(list)
    allok = True
    for r in csv.DictReader(open(os.path.join(RES, "scaling.csv"))):
        data[(int(r["n"]), int(r["k"]), r["mode"], r["algo"])].append(float(r["time_ms"]))
        allok = allok and r["ok"] == "1"
    assert allok, "a scaling run failed validation"
    med = lambda key: st.median(data[key]) if key in data else None
    keys = sorted({(n, k) for (n, k, _, _) in data})
    rows, last = [], None
    md.append("\n### Scaling on unimodal instances (median ms over 3 instances)\n\n| n | k | DP (64-bit) | O1 (64-bit) | O2 (64-bit) | O2 (128-bit) |\n|---|---|---|---|---|---|")
    for n, k in keys:
        if last is not None and n != last:
            rows.append("\\midrule")
        last = n
        v = [med((n, k, "narrow", "DP")), med((n, k, "narrow", "O1")), med((n, k, "narrow", "O2")), med((n, k, "wide", "O2"))]
        rows.append(f"{n} & {k} & " + " & ".join(fmt(x) for x in v) + " \\\\")
        md.append(f"| {n} | {k} | " + " | ".join(fmt(x) for x in v) + " |")
    open(os.path.join(OUT, "arc_scaling.tex"), "w").write("\n".join(rows) + "\n")

    def slope(mode, algo, n0, n1):
        s = []
        for k in (4, 16, 64):
            a, b = med((n0, k, mode, algo)), med((n1, k, mode, algo))
            if a and b:
                s.append(math.log(b / a) / math.log(n1 / n0))
        return st.median(s)

    mac("arcSlope", f"{slope('wide', 'O2', 1000, 1000000):.2f}")
    mac("arcSlopeOone", f"{slope('narrow', 'O1', 1000, 10000):.2f}")
    mac("arcMillion", f"{st.median([med((1000000, k, 'wide', 'O2')) for k in (4, 16, 64, 100000)]) / 1000:.0f}")
    mac("arcSpeedupDP", f"{med((400, 16, 'narrow', 'DP')) / med((400, 16, 'narrow', 'O2')):.0f}")
    md.append(f"\nLog-log slopes: O2 (128-bit, n = 10^3..10^6): {slope('wide', 'O2', 1000, 1000000):.2f}; "
              f"O2 (64-bit, n = 10^3..10^4): {slope('narrow', 'O2', 1000, 10000):.2f}; "
              f"O1 (n = 10^3..10^4): {slope('narrow', 'O1', 1000, 10000):.2f}; "
              f"DP (n = 100..400): {slope('narrow', 'DP', 100, 400):.2f}.\n")
    open(os.path.join(OUT, "arc_macros.tex"), "w").write("\n".join(macros) + "\n")
    open(os.path.join(RES, "tables.md"), "w").write("\n".join(md) + "\n")
    print("wrote", os.path.normpath(OUT), "and results/tables.md")


if __name__ == "__main__":
    main()
