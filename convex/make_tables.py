#!/usr/bin/env python3
"""Turn the raw experiment outputs in results/ into the LaTeX tables, plot
data and number macros used by the paper (written to ../paper/generated/).

    python3 make_tables.py
"""
import collections
import csv
import os
import re
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "..", "paper", "generated")
ORDER = ["square", "gauss", "disk", "annulus", "thin", "circle", "cocircular", "twochain"]


def load(name):
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        return []
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        m = re.search(r"([a-z]+)_n(\d+)_k(\d+)_s(\d+)", r["file"])
        r["dist"] = m.group(1)
        for key in r:
            if key.startswith(("t_", "dec_", "cand_", "far_", "incircle", "greedy_calls")) or key in ("n", "k"):
                r[key] = float(r[key]) if "." in r[key] or key.startswith("t_") else int(r[key])
    return rows


def med(xs):
    xs = [x for x in xs if x >= 0]
    return st.median(xs) if xs else None


def fmt(x):
    if x is None:
        return "--"
    if x >= 100:
        return f"{x:.0f}"
    if x >= 10:
        return f"{x:.1f}"
    return f"{x:.2f}"


def write(name, text):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w") as f:
        f.write(text)


def macro(lines, name, value):
    lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")


def main():
    macros = []
    verify, medium = load("verify.csv"), load("medium.csv")
    scaling = load("scaling_small.csv") + load("scaling_large.csv")
    ksweep = load("ksweep.csv")

    # ---- correctness -----------------------------------------------------
    pyfmm = {}
    p = os.path.join(RES, "verify_python.txt")
    if os.path.exists(p):
        m = re.search(r"(\d+)/(\d+) instances agree with Python brute force; FMM solver agrees on (\d+)/(\d+)", open(p).read())
        if m:
            pyfmm = {"py": int(m.group(1)), "fmm": int(m.group(3)), "n": int(m.group(2))}
    lines = []
    for split, rows in (("verify", verify), ("medium", medium)):
        by = collections.defaultdict(list)
        for r in rows:
            by[r["dist"]].append(r)
        for d in ORDER:
            if d not in by:
                continue
            rs = by[d]
            ok = sum(r["ok"] == "1" for r in rs)
            brute = sum(r["r2_brute"] != "-1" for r in rs)
            ns = [r["n"] for r in rs]
            lines.append(f"\\textit{{{split}}} & \\texttt{{{d}}} & {len(rs)} & {min(ns)}--{max(ns)} & {ok}/{len(rs)} & {brute} \\\\")
        if split == "verify":
            lines.append("\\midrule")
    write("table_correctness.tex", "\n".join(lines) + "\n")
    macro(macros, "nVerify", len(verify))
    macro(macros, "nVerifyOk", sum(r["ok"] == "1" for r in verify))
    macro(macros, "nMedium", len(medium))
    macro(macros, "nMediumOk", sum(r["ok"] == "1" for r in medium))
    macro(macros, "nMediumBrute", sum(r["r2_brute"] != "-1" for r in medium))
    macro(macros, "nVerifyPy", pyfmm.get("py", "??"))
    macro(macros, "nVerifyFmm", pyfmm.get("fmm", "??"))

    # ---- running time ----------------------------------------------------
    macro(macros, "nScaling", len(scaling))
    macro(macros, "nScalingOk", sum(r["ok"] == "1" for r in scaling))
    cells = collections.defaultdict(list)
    for r in scaling:
        cells[(r["k"], r["n"])].append(r)
    lines = []
    plots = collections.defaultdict(lambda: collections.defaultdict(list))
    lastk = None
    for (k, n) in sorted(cells):
        rs = cells[(k, n)]
        tt, tu, tn, tc, tz = (med([r[c] for r in rs]) for c in ("t_triple", "t_uncapped", "t_nowindow", "t_capped", "t_certified"))
        dn, dw, dz = (st.mean([r[c] for r in rs]) for c in ("dec_nowindow", "dec_window", "dec_certified"))
        if lastk is not None and k != lastk:
            lines.append("\\midrule")
        lastk = k
        lines.append(f"{k} & {n} & {fmt(tt)} & {fmt(tu)} & {fmt(tn)} & {fmt(tc)} & {fmt(tz)} & "
                     f"{dn:.1f} & {dw:.1f} & {dz:.1f} \\\\")
        for name, v in (("triple", tt), ("uncapped", tu), ("capped", tc), ("certified", tz)):
            if v is not None:
                plots[k][name].append((n, v))
    write("table_scaling.tex", "\n".join(lines) + "\n")
    for k in plots:
        for name, pts in plots[k].items():
            write(f"plot_k{k}_{name}.dat", "n t\n" + "".join(f"{n} {v:.4f}\n" for n, v in pts))

    # empirical growth exponent of the capped solver between n=100 and n=400
    import math
    for k in (4, 8, 16, 32):
        a, b = dict(plots[k]["capped"]).get(100), dict(plots[k]["capped"]).get(400)
        if a and b:
            macro(macros, "expCapped" + {4: "Four", 8: "Eight", 16: "Sixteen", 32: "ThirtyTwo"}[k], f"{math.log(b / a) / math.log(4):.2f}")

    # ---- window and decisions -------------------------------------------
    if scaling:
        frac = [100.0 * r["cand_window"] / r["cand_total"] for r in scaling]
        macro(macros, "windowMedian", f"{st.median(frac):.1f}")
        macro(macros, "windowMax", f"{max(frac):.1f}")
        macro(macros, "windowSaved", f"{st.mean([r['dec_nowindow'] - r['dec_window'] for r in scaling]):.1f}")
        sp = [r["t_nowindow"] / r["t_capped"] for r in scaling if r["t_capped"] > 0.5]
        macro(macros, "windowSpeedup", f"{st.median(sp):.2f}")
        tri = [100.0 * r["far_triples_at_opt"] / (r["n"] * (r["n"] - 1) * (r["n"] - 2) / 6) for r in scaling]
        macro(macros, "farTripleMedian", f"{st.median(tri):.0f}")
        macro(macros, "farTripleMin", f"{min(tri):.0f}")
        macro(macros, "farTripleMax", f"{max(tri):.0f}")
        sp = [r["t_uncapped"] / r["t_capped"] for r in scaling if r["t_uncapped"] > 0 and r["n"] == 200 and r["k"] == 4]
        macro(macros, "capSpeedupKFour", f"{st.median(sp):.1f}")
        sp = [r["t_triple"] / r["t_capped"] for r in scaling if r["t_triple"] > 0 and r["n"] == 100]
        macro(macros, "tripleSpeedupHundred", f"{st.median(sp):.1f}")
        sp = [r["t_capped"] / r["t_certified"] for r in scaling if r["certified_by_greedy"] == "1" and r["t_certified"] > 0.5]
        macro(macros, "certSpeedup", f"{st.median(sp):.1f}")
        q = st.quantiles(sp, n=4)
        macro(macros, "certSpeedupLo", f"{q[0]:.1f}")
        macro(macros, "certSpeedupHi", f"{q[2]:.1f}")
        ov = [r["t_certified"] / r["t_capped"] for r in scaling if r["certified_by_greedy"] != "1" and r["t_capped"] > 0.5]
        macro(macros, "certOverhead", f"{100 * (st.median(ov) - 1):.0f}" if ov else "0")
        macro(macros, "certScalingNo", sum(r["certified_by_greedy"] != "1" for r in scaling))

    # ---- certification rate per family ------------------------------------
    lines = []
    allrows = scaling + medium + verify
    by = collections.defaultdict(list)
    for r in allrows:
        by[r["dist"]].append(r)
    tot_c = tot = 0
    for d in ORDER:
        if d not in by:
            continue
        rs = by[d]
        c = sum(r["certified_by_greedy"] == "1" for r in rs)
        ff = sum(r["greedy_ff_r2"] == r["r2"] for r in rs)
        tot_c += c
        tot += len(rs)
        gap = [100.0 * (1 - (int(r["greedy_ff_r2"]) / int(r["r2"])) ** 0.5) for r in rs]
        lines.append(f"\\texttt{{{d}}} & {len(rs)} & {ff} & {c} & {100.0 * c / len(rs):.0f}\\% & {st.mean(gap):.1f}\\% & {max(gap):.1f}\\% \\\\")
    write("table_certified.tex", "\n".join(lines) + "\n")
    macro(macros, "certTotal", tot)
    macro(macros, "certYes", tot_c)
    tw = [r for r in scaling if r["dist"] == "twochain"]
    if tw:
        macro(macros, "certTwochainNo", sum(r["certified_by_greedy"] != "1" for r in tw))
        macro(macros, "certTwochainAll", len(tw))

    # ---- k sweep ------------------------------------------------------------
    lines = []
    by = collections.defaultdict(list)
    for r in ksweep:
        by[r["k"]].append(r)
    for k in sorted(by):
        rs = by[k]
        tu, tc, tz = (med([r[c] for r in rs]) for c in ("t_uncapped", "t_capped", "t_certified"))
        inc = st.median([r["incircle_per_decision"] for r in rs]) / 1e6
        lines.append(f"{k} & {fmt(tu)} & {fmt(tc)} & {tu / tc:.2f}$\\times$ & {fmt(tz)} & {inc:.1f} \\\\")
    write("table_ksweep.tex", "\n".join(lines) + "\n")

    # ---- random property tests ---------------------------------------------
    tot = collections.defaultdict(lambda: [0, 0])
    kfail = collections.defaultdict(set)
    krange = collections.defaultdict(set)
    for mode in ("exh", "big", "pert"):
        p = os.path.join(RES, f"hunt_{mode}.txt")
        if not os.path.exists(p):
            continue
        for line in open(p):
            f = line.split()
            if f[0] != "STAT":
                continue
            prop, k, tests, fails = f[2], int(f[3]), int(f[4]), int(f[5])
            tot[prop][0] += tests
            tot[prop][1] += fails
            krange[prop].add(k)
            if fails:
                kfail[prop].add(k)
    names = [("G", "$(\\mathrm{G}_k)$ arc characterisation"),
             ("T3", "single replacement, \\cref{thm:k3}(ii)"),
             ("C1p", "$(\\mathrm{TT}_k)$ tight triple"),
             ("V1", "$(\\mathrm{AG}'_k)$ alternating greedy, 3-point seed"),
             ("V3", "$(\\mathrm{ST}'_k)$ search tree, 3-point seed"),
             ("V5", "$(\\mathrm{AG}_k)$ alternating greedy"),
             ("C1", "$(\\mathrm{TP}_k)$ tight pair"),
             ("V4", "$(\\mathrm{ST}_k)$ search tree")]
    lines = []
    for key, label in names:
        if key not in tot:
            continue
        t, f = tot[key]
        ks = krange[key]
        kr = f"{min(ks)}--{max(ks)}" if len(ks) > 1 else f"{min(ks)}"
        kf = "--" if not kfail[key] else (f"{min(kfail[key])}--{max(kfail[key])}" if len(kfail[key]) > 1 else f"{min(kfail[key])}")
        num = lambda v: f"{v:,}".replace(",", "{,}")
        lines.append(f"{label} & {kr} & {num(t)} & {num(f)} & {kf} \\\\")
    write("table_hunt.tex", "\n".join(lines) + "\n")
    if "C1" in tot:
        macro(macros, "huntCone", f"{tot['C1'][0]:,}".replace(",", "{,}"))
        macro(macros, "huntVfour", f"{tot['V4'][0]:,}".replace(",", "{,}"))
    write("macros.tex", "\n".join(macros) + "\n")
    print(f"wrote tables to {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
