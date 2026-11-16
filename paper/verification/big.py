import sys, random, itertools, subprocess
sys.path.insert(0, sys.argv[1])
from dispersion_check import *
rng = random.Random(7); t = m = fm = fn = 0; tri = 0
for it in range(400):
    kind = ('square','disk','circle','thin')[it % 4]
    n = rng.randint(10, 16); k = rng.randint(3, min(9, n))
    P = rand_convex(rng, n, kind)
    bf, _ = brute_opt(P, k); cp, w = capped_opt(P, k); t += 1; m += bf != cp
    if it % 4 == 0:
        inp = f"{n} 2 {k}\n" + "\n".join(f"{x} {y}" for x, y in P) + "\n"
        o = subprocess.run(['./fmm'], input=inp, capture_output=True, text=True).stdout.split()
        fn += 1; fm += int(o[1]) != bf
print(f"n in [10,16], k in [3,9]: {t} instances, capped-vs-brute mismatches {m}; repo FMM solver checked on {fn}, mismatches {fm}")
