import sys, random, itertools
sys.path.insert(0, sys.argv[1])
from dispersion_check import *
rng = random.Random(99); T = {}; F = {}
for it in range(400):
    P = rand_convex(rng, rng.randint(6, 11), ('square','disk','circle','thin')[it % 4])
    for k in range(3, min(6, len(P)) + 1):
        for r2 in sorted({d2(P[a], P[b]) for a, b in itertools.combinations(range(len(P)), 2)}):
            T[k] = T.get(k, 0) + 1
            F[k] = F.get(k, 0) + (bool(far_sets(P, k, r2)) != G(P, k, r2))
print({k: (T[k], F[k]) for k in sorted(T)})
