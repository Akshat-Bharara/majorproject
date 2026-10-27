#!/usr/bin/env python3
"""Cross-verify the FMM k-dispersion solver against the exhaustive oracle.

Steps
  1. build ../dispersion (FMM solver) and brute_force (oracle) with g++
  2. run both on every instance listed in dataset/manifest.csv
  3. for each instance check
       - both report the same optimal squared distance,
       - each returned set has k distinct valid indices,
       - the minimum pairwise squared distance of each returned set,
         recomputed here in Python, equals the reported optimum
  4. scaling experiment: running time vs n for fixed k
  5. write results/*.csv and LaTeX fragments in ../report/generated/

Usage: python3 verify/run_verification.py [--skip-scaling]
"""

import csv
import itertools
import os
import random
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BUILD = os.path.join(HERE, "build")
RESULTS = os.path.join(HERE, "results")
GENERATED = os.path.join(REPO, "report", "generated")

FMM_EXE = os.path.join(BUILD, "dispersion")
BF_EXE = os.path.join(BUILD, "brute_force")

PAPER_SPLIT = {2: (0, 1, 1), 3: (1, 1, 1), 4: (1, 1, 2), 5: (2, 1, 2),
               6: (2, 2, 2), 7: (2, 3, 2), 8: (2, 2, 4)}


def build():
    os.makedirs(BUILD, exist_ok=True)
    srcs = [os.path.join(REPO, f) for f in
            ("main.cpp", "dispersion.cpp", "clique.cpp", "matrix.cpp")]
    flags = ["g++", "-O2", "-std=c++17", "-Wall", "-Wextra"]
    subprocess.run(flags + srcs + ["-o", FMM_EXE], check=True)
    subprocess.run(flags + [os.path.join(HERE, "brute_force.cpp"), "-o", BF_EXE],
                   check=True)


def run(exe, text, args=()):
    """Run a solver; return (parsed output, in-process solve time in ms)."""
    p = subprocess.run([exe, *args], input=text, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{exe} failed: {p.stderr}")
    out = {}
    for line in p.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "optimal_squared_distance":
            out["opt"] = int(parts[1])
        elif parts and parts[0] == "selected_indices":
            out["sel"] = [int(x) for x in parts[1:]]
    for line in p.stderr.splitlines():
        parts = line.split()
        if parts:
            out[parts[0]] = [float(x) if "." in x or "e" in x else int(x)
                             for x in parts[1:]]
    return out, float(out["solve_ms"][0])


def parse_instance(text):
    tok = text.split()
    n, d, k = int(tok[0]), int(tok[1]), int(tok[2])
    vals = list(map(int, tok[3:]))
    pts = [vals[i * d:(i + 1) * d] for i in range(n)]
    return n, d, k, pts


def sq(p, q):
    return sum((a - b) ** 2 for a, b in zip(p, q))


def set_value(pts, sel):
    return min(sq(pts[a], pts[b]) for a, b in itertools.combinations(sel, 2))


def valid_set(n, k, sel):
    return len(sel) == k and len(set(sel)) == k and all(0 <= v < n for v in sel)


def verify_dataset():
    with open(os.path.join(HERE, "dataset", "manifest.csv")) as fh:
        manifest = list(csv.DictReader(fh))

    rows = []
    for m in manifest:
        with open(os.path.join(HERE, "dataset", "instances", m["id"] + ".txt")) as fh:
            text = fh.read()
        n, d, k, pts = parse_instance(text)
        fmm, t_fmm = run(FMM_EXE, text, ["--stats"])
        bf, t_bf = run(BF_EXE, text)

        fmm_ok = valid_set(n, k, fmm["sel"]) and set_value(pts, fmm["sel"]) == fmm["opt"]
        bf_ok = valid_set(n, k, bf["sel"]) and set_value(pts, bf["sel"]) == bf["opt"]
        split = tuple(fmm["split"])
        rows.append({
            "id": m["id"], "family": m["family"], "n": n, "d": d, "k": k,
            "split": "-".join(map(str, split)),
            "split_matches_paper": PAPER_SPLIT.get(k, split) == split,
            "opt_fmm": fmm["opt"], "opt_bf": bf["opt"],
            "opt_match": fmm["opt"] == bf["opt"],
            "fmm_set_valid": fmm_ok, "bf_set_valid": bf_ok,
            "same_set": sorted(fmm["sel"]) == sorted(bf["sel"]),
            "candidate_distances": fmm["candidate_distances"][0],
            "feasibility_tests": fmm["feasibility_tests"][0],
            "max_V1": fmm["max_V"][0], "max_V2": fmm["max_V"][1],
            "max_V3": fmm["max_V"][2],
            "t_fmm_ms": round(t_fmm, 3), "t_bf_ms": round(t_bf, 3),
        })
        status = "ok" if rows[-1]["opt_match"] and fmm_ok and bf_ok else "MISMATCH"
        print(f"{m['id']:<28} opt={fmm['opt']:<22} {status}")
    return rows


def scaling():
    """Median solve time over 3 runs, uniform points in the plane."""
    rng = random.Random(7)
    plan = {3: range(25, 201, 25), 4: range(10, 81, 10), 6: range(10, 35, 4)}
    rows = []
    for k, ns in plan.items():
        for n in ns:
            pts = [[rng.randint(0, 10 ** 6), rng.randint(0, 10 ** 6)] for _ in range(n)]
            text = f"{n} 2 {k}\n" + "\n".join(f"{x} {y}" for x, y in pts) + "\n"
            tf, tb = [], []
            for _ in range(3):
                a, t1 = run(FMM_EXE, text, ["--stats"])
                b, t2 = run(BF_EXE, text)
                assert a["opt"] == b["opt"], (k, n)
                tf.append(t1)
                tb.append(t2)
            rows.append({"k": k, "n": n,
                         "t_fmm_ms": round(statistics.median(tf), 3),
                         "t_bf_ms": round(statistics.median(tb), 3),
                         "feasibility_tests": a["feasibility_tests"][0]})
            print(f"scaling k={k} n={n} fmm={rows[-1]['t_fmm_ms']}ms "
                  f"bf={rows[-1]['t_bf_ms']}ms")
    return rows


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def yes(b):
    return "\\checkmark" if b else "$\\times$"


def write_latex(rows, scal):
    os.makedirs(GENERATED, exist_ok=True)
    total = len(rows)
    match = sum(r["opt_match"] for r in rows)
    fv = sum(r["fmm_set_valid"] for r in rows)
    bv = sum(r["bf_set_valid"] for r in rows)
    same = sum(r["same_set"] for r in rows)
    with open(os.path.join(GENERATED, "summary_macros.tex"), "w") as fh:
        fh.write(f"\\newcommand{{\\NumInstances}}{{{total}}}\n")
        fh.write(f"\\newcommand{{\\NumOptMatch}}{{{match}}}\n")
        fh.write(f"\\newcommand{{\\NumFmmValid}}{{{fv}}}\n")
        fh.write(f"\\newcommand{{\\NumBfValid}}{{{bv}}}\n")
        fh.write(f"\\newcommand{{\\NumSameSet}}{{{same}}}\n")
        fh.write(f"\\newcommand{{\\NumFamilies}}{{{len({r['family'] for r in rows})}}}\n")

    # Per-k table.
    with open(os.path.join(GENERATED, "table_by_k.tex"), "w") as fh:
        for k in sorted({r["k"] for r in rows}):
            g = [r for r in rows if r["k"] == k]
            fh.write(
                f"{k} & ({g[0]['split'].replace('-', ',')}) & {len(g)} & "
                f"{min(r['n'] for r in g)}--{max(r['n'] for r in g)} & "
                f"{sum(r['opt_match'] for r in g)}/{len(g)} & "
                f"{sum(r['fmm_set_valid'] for r in g)}/{len(g)} & "
                f"{statistics.mean(r['feasibility_tests'] for r in g):.1f} & "
                f"{max(max(r['max_V1'], r['max_V2'], r['max_V3']) for r in g)} & "
                f"{statistics.mean(r['t_fmm_ms'] for r in g):.1f} & "
                f"{statistics.mean(r['t_bf_ms'] for r in g):.1f} \\\\\n")

    # Per-family table.
    with open(os.path.join(GENERATED, "table_by_family.tex"), "w") as fh:
        for fam in dict.fromkeys(r["family"] for r in rows):
            g = [r for r in rows if r["family"] == fam]
            dims = ",".join(map(str, sorted({r["d"] for r in g})))
            fh.write(
                f"\\texttt{{{fam}}} & {dims} & {len(g)} & "
                f"{sum(r['opt_match'] for r in g)}/{len(g)} & "
                f"{sum(r['same_set'] for r in g)}/{len(g)} & "
                f"{statistics.mean(r['candidate_distances'] for r in g):.0f} \\\\\n")

    # A sample of individual instances (one per k, uniform d=2, largest n).
    with open(os.path.join(GENERATED, "table_sample.tex"), "w") as fh:
        for k in sorted({r["k"] for r in rows}):
            g = [r for r in rows if r["k"] == k and r["family"] == "uniform" and r["d"] == 2]
            r = max(g, key=lambda r: r["n"])
            fh.write(
                f"\\texttt{{{r['id'].replace('_', '\\_')}}} & {r['n']} & {r['k']} & "
                f"{r['opt_fmm']} & {r['opt_bf']} & {yes(r['opt_match'])} & "
                f"{yes(r['fmm_set_valid'])} & {r['feasibility_tests']} \\\\\n")

    for k in sorted({r["k"] for r in scal}):
        with open(os.path.join(GENERATED, f"scaling_k{k}.dat"), "w") as fh:
            fh.write("n fmm bf\n")
            for r in scal:
                if r["k"] == k:
                    fh.write(f"{r['n']} {r['t_fmm_ms']} {r['t_bf_ms']}\n")


def main():
    build()
    os.makedirs(RESULTS, exist_ok=True)
    rows = verify_dataset()
    write_csv(os.path.join(RESULTS, "results.csv"), rows)
    scal = [] if "--skip-scaling" in sys.argv else scaling()
    if scal:
        write_csv(os.path.join(RESULTS, "scaling.csv"), scal)
    write_latex(rows, scal)

    bad = [r["id"] for r in rows
           if not (r["opt_match"] and r["fmm_set_valid"] and r["bf_set_valid"])]
    print(f"\n{len(rows) - len(bad)}/{len(rows)} instances verified")
    if bad:
        print("mismatches:", ", ".join(bad))
        sys.exit(1)


if __name__ == "__main__":
    main()
