#!/usr/bin/env python3
"""Cross-verification and timing for k-dispersion in convex position.

Solvers compared (all read the same "n 2 k" files):
  capped  convex/convex_dispersion --algo=capped   (this project, Algorithm 1)
  certified  --algo=certified: seeded greedy + one exact decision (Algorithm 2)
  tw      convex/convex_dispersion --algo=tw       (Tkachenko-Wang DP, O(n^4)/decision)
  brute   verify/brute_force                       (exhaustive oracle)
  fmm     ./dispersion                             (general FMM solver, small n only)

Checks per instance: all solvers report the same optimum, and every returned
index set has k distinct points whose minimum pairwise squared distance,
recomputed here exactly, equals that optimum.

Usage: python3 convex/run_convex_experiments.py [--quick]
"""

import csv
import itertools
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BUILD = os.path.join(HERE, "build")
RESULTS = os.path.join(HERE, "results")
GEN = os.path.join(HERE, "paper", "generated")

CVX = os.path.join(BUILD, "convex_dispersion")
BRUTE = os.path.join(BUILD, "brute_force")
FMM = os.path.join(BUILD, "dispersion")


def build():
    os.makedirs(BUILD, exist_ok=True)
    cc = ["g++", "-O2", "-std=c++17", "-Wall", "-Wextra"]
    subprocess.run(cc + [os.path.join(HERE, f) for f in
                         ("geometry.cpp", "convex_dispersion.cpp", "main_convex.cpp")]
                   + ["-o", CVX], check=True)
    subprocess.run(cc + [os.path.join(REPO, "verify", "brute_force.cpp"), "-o", BRUTE],
                   check=True)
    subprocess.run(cc + [os.path.join(REPO, f) for f in
                         ("main.cpp", "dispersion.cpp", "clique.cpp", "matrix.cpp")]
                   + ["-o", FMM], check=True)


def run(cmd, text, timeout=None):
    p = subprocess.run(cmd, input=text, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(f"{cmd} failed: {p.stderr}")
    out = {}
    for line in (p.stdout + p.stderr).splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "optimal_squared_distance":
            out["opt"] = int(parts[1])
        elif parts[0] == "selected_indices":
            out["sel"] = [int(x) for x in parts[1:]]
        elif len(parts) == 2:
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return out


SOLVERS = {
    "capped": lambda: [CVX, "--algo=capped", "--stats"],
    "capped_nobounds": lambda: [CVX, "--algo=capped", "--no-bounds", "--stats"],
    "certified": lambda: [CVX, "--algo=certified", "--stats"],
    "tw": lambda: [CVX, "--algo=tw", "--stats"],
    "tw_nobounds": lambda: [CVX, "--algo=tw", "--no-bounds", "--stats"],
    "brute": lambda: [BRUTE],
    "fmm": lambda: [FMM, "--stats"],
}


def load(split):
    with open(os.path.join(HERE, "dataset", split, "manifest.csv")) as fh:
        man = list(csv.DictReader(fh))
    for m in man:
        with open(os.path.join(HERE, "dataset", split, "instances", m["id"] + ".txt")) as fh:
            text = fh.read()
        tok = text.split()
        n, k = int(tok[0]), int(tok[2])
        v = list(map(int, tok[3:]))
        pts = [(v[2 * i], v[2 * i + 1]) for i in range(n)]
        yield m, text, n, k, pts


def value(pts, sel):
    return min((pts[a][0] - pts[b][0]) ** 2 + (pts[a][1] - pts[b][1]) ** 2
               for a, b in itertools.combinations(sel, 2))


def valid(n, k, pts, out):
    s = out.get("sel", [])
    return (len(s) == k and len(set(s)) == k and all(0 <= x < n for x in s)
            and value(pts, s) == out["opt"])


def verify_split(split, solvers):
    rows = []
    for m, text, n, k, pts in load(split):
        outs = {name: run(SOLVERS[name](), text) for name in solvers}
        ref = outs["brute"]["opt"]
        row = {"id": m["id"], "dist": m["dist"], "n": n, "k": k, "opt": ref}
        ok = True
        for name, o in outs.items():
            row[f"{name}_match"] = o["opt"] == ref
            row[f"{name}_valid"] = valid(n, k, pts, o)
            ok &= row[f"{name}_match"] and row[f"{name}_valid"]
        row["all_ok"] = ok
        rows.append(row)
        print(f"[{split}] {m['id']:<26} opt={ref:<14} {'ok' if ok else 'MISMATCH'}")
    return rows


def scaling(quick):
    rows = []
    for m, text, n, k, pts in load("scaling"):
        if quick and n > 100:
            continue
        row = {"id": m["id"], "dist": m["dist"], "n": n, "k": k}
        c = run(SOLVERS["capped"](), text)
        cn = run(SOLVERS["capped_nobounds"](), text)
        ce = run(SOLVERS["certified"](), text)
        assert ce["opt"] == c["opt"], m["id"]
        row.update(opt=c["opt"], capped_ms=c["solve_ms"], capped_dec=c["decisions"],
                   capped_steps=c["inner_steps"], cand_total=c["candidates_total"],
                   cand_window=c["candidates_window"],
                   capped_nb_ms=cn["solve_ms"], capped_nb_dec=cn["decisions"],
                   cert_ms=ce["solve_ms"], cert_dec=ce["decisions"],
                   cert_greedy=ce["greedy_decisions"], cert_ok=ce["certified_by_greedy"])
        assert cn["opt"] == c["opt"]
        if n <= 100:
            t = run(SOLVERS["tw"](), text)
            assert t["opt"] == c["opt"], m["id"]
            row.update(tw_ms=t["solve_ms"], tw_steps=t["inner_steps"])
        else:
            row.update(tw_ms="", tw_steps="")
        row["capped_valid"] = valid(n, k, pts, c) and valid(n, k, pts, ce)
        rows.append(row)
        print(f"[scaling] {m['id']:<26} capped={row['capped_ms']:.1f}ms cert={row['cert_ms']:.1f}ms "
              f"tw={row['tw_ms'] if row['tw_ms'] == '' else round(row['tw_ms'], 1)}ms "
              f"window {row['cand_window']}/{row['cand_total']}")
    return rows


def write_csv(name, rows):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def tex_rows(rows, solvers):
    lines = []
    for dist in dict.fromkeys(r["dist"] for r in rows):
        g = [r for r in rows if r["dist"] == dist]
        cells = [f"\\texttt{{{dist}}}", str(len(g))]
        for s in solvers:
            cells.append(f"{sum(r[f'{s}_match'] and r[f'{s}_valid'] for r in g)}/{len(g)}")
        lines.append(" & ".join(cells) + " \\\\")
    return "\n".join(lines) + "\n"


def main():
    quick = "--quick" in sys.argv
    build()
    os.makedirs(GEN, exist_ok=True)

    v = verify_split("verify", ["capped", "certified", "tw", "brute", "fmm"])
    write_csv("verify.csv", v)
    md = verify_split("medium", ["capped", "certified", "tw", "brute"])
    write_csv("medium.csv", md)
    sc = scaling(quick)
    write_csv("scaling.csv", sc)

    with open(os.path.join(GEN, "verify_by_dist.tex"), "w") as fh:
        fh.write(tex_rows(v, ["capped", "certified", "tw", "fmm"]))
    with open(os.path.join(GEN, "medium_by_dist.tex"), "w") as fh:
        fh.write(tex_rows(md, ["capped", "certified", "tw"]))
    with open(os.path.join(GEN, "macros.tex"), "w") as fh:
        fh.write(f"\\newcommand{{\\NVerify}}{{{len(v)}}}\n")
        fh.write(f"\\newcommand{{\\NVerifyOk}}{{{sum(r['all_ok'] for r in v)}}}\n")
        fh.write(f"\\newcommand{{\\NMedium}}{{{len(md)}}}\n")
        fh.write(f"\\newcommand{{\\NMediumOk}}{{{sum(r['all_ok'] for r in md)}}}\n")
        fh.write(f"\\newcommand{{\\NScaling}}{{{len(sc)}}}\n")
        win = [r["cand_window"] / r["cand_total"] for r in sc]
        dec = [r["capped_nb_dec"] - r["capped_dec"] for r in sc]
        fh.write(f"\\newcommand{{\\WindowMedianPct}}{{{100 * statistics.median(win):.1f}}}\n")
        fh.write(f"\\newcommand{{\\WindowMaxPct}}{{{100 * max(win):.1f}}}\n")
        fh.write(f"\\newcommand{{\\DecisionsSavedMean}}{{{statistics.mean(dec):.1f}}}\n")
        fh.write(f"\\newcommand{{\\CertRate}}{{{sum(r['cert_ok'] == 1 for r in sc)}}}\n")
        allc = [r for r in v + md]
        fh.write(f"\\newcommand{{\\CertExactSmall}}{{{sum(r['certified_match'] and r['certified_valid'] for r in allc)}}}\n")
        fh.write(f"\\newcommand{{\\NSmallAll}}{{{len(allc)}}}\n")

    # Timing table and plot data: one row per (n, k), median over distributions.
    with open(os.path.join(GEN, "scaling_table.tex"), "w") as fh:
        for k in sorted({r["k"] for r in sc}):
            for n in sorted({r["n"] for r in sc if r["k"] == k}):
                g = [r for r in sc if r["k"] == k and r["n"] == n]
                cm = statistics.median(r["capped_ms"] for r in g)
                cem = statistics.median(r["cert_ms"] for r in g)
                tw = [r["tw_ms"] for r in g if r["tw_ms"] != ""]
                twm = f"{statistics.median(tw):.1f}" if tw else "--"
                sp = f"{statistics.median(t / max(r['capped_ms'], 1e-3) for r, t in ((r, r['tw_ms']) for r in g if r['tw_ms'] != '')):.1f}$\\times$" if tw else "--"
                fh.write(f"{k} & {n} & {twm} & {cm:.1f} & {cem:.1f} & {sp} & "
                         f"{statistics.median(r['capped_nb_dec'] for r in g):.0f} & "
                         f"{statistics.median(r['capped_dec'] for r in g):.0f} & "
                         f"{statistics.median(r['cert_dec'] for r in g):.0f} \\\\\n")
    for k in sorted({r["k"] for r in sc}):
        with open(os.path.join(GEN, f"scaling_k{k}.dat"), "w") as fh:
            fh.write("n capped tw certified\n")
            for n in sorted({r["n"] for r in sc if r["k"] == k}):
                g = [r for r in sc if r["k"] == k and r["n"] == n]
                tw = [r["tw_ms"] for r in g if r["tw_ms"] != ""]
                fh.write(f"{n} {statistics.median(r['capped_ms'] for r in g):.3f} "
                         f"{statistics.median(tw) if tw else 'nan'} "
                         f"{statistics.median(r['cert_ms'] for r in g):.3f}\n")

    bad = [r["id"] for r in v + md if not r["all_ok"]]
    bad += [r["id"] for r in sc if not r["capped_valid"]]
    print(f"\nverify {sum(r['all_ok'] for r in v)}/{len(v)}, "
          f"medium {sum(r['all_ok'] for r in md)}/{len(md)}, scaling {len(sc)} instances")
    if bad:
        print("FAILURES:", bad)
        sys.exit(1)


if __name__ == "__main__":
    main()
