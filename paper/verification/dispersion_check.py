"""Independent check of the paper's claims, exact integer arithmetic only."""
import itertools, random, math, sys

def d2(p, q): return (p[0]-q[0])**2 + (p[1]-q[1])**2
def orient(a, b, c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
def incircle(a, b, c, d):
    """>0 iff d strictly inside circle(a,b,c), with a,b,c counter-clockwise."""
    rows = []
    for p in (a, b, c):
        dx, dy = p[0]-d[0], p[1]-d[1]; rows.append((dx, dy, dx*dx+dy*dy))
    (a1,a2,a3),(b1,b2,b3),(c1,c2,c3) = rows
    return a1*(b2*c3-b3*c2) - a2*(b1*c3-b3*c1) + a3*(b1*c2-b2*c1)

def strict_hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3: return pts
    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and orient(h[-2], h[-1], p) <= 0: h.pop()
            h.append(p)
        return h
    lo, up = half(pts), half(pts[::-1])
    return lo[:-1] + up[:-1]          # ccw, collinear points removed

def is_strictly_convex_ccw(P):
    n = len(P)
    return all(orient(P[i], P[(i+1)%n], P[(i+2)%n]) > 0 for i in range(n))

# ---------- brute force ----------
def brute_opt(P, k):
    best, arg = -1, None
    for S in itertools.combinations(range(len(P)), k):
        m = min(d2(P[a], P[b]) for a, b in itertools.combinations(S, 2))
        if m > best: best, arg = m, S
    return best, arg

def far_sets(P, k, r2):
    return [S for S in itertools.combinations(range(len(P)), k)
            if all(d2(P[a], P[b]) >= r2 for a, b in itertools.combinations(S, 2))]

# ---------- Algorithm 1 (CappedDecision), implemented as written ----------
def capped_decision(P, k, r2, stats=None):
    n = len(P); K = k - 2
    assert k >= 3
    far = lambda a, b: d2(P[a], P[b]) >= r2
    L = {}
    def count(chord, i, l, j):
        """#c with L_c(chord) not strictly inside circle(p_i,p_l,p_j); also checks monotonicity (Lemma 7)."""
        Lc = L.get(chord)
        if Lc is None: return 0
        ok = [Lc[c] is not None and incircle(P[i], P[l], P[j], P[Lc[c]]) <= 0 for c in range(K)]
        q = sum(ok)
        assert ok == [True]*q + [False]*(K-q), "threshold monotonicity violated"
        return q
    for gap in range(2, n):
        for i in range(n):
            j = (i + gap) % n
            if not far(i, j): continue
            best = [None]*(K+1)
            for t in range(1, gap):
                l = (i + t) % n
                if not (far(i, l) and far(l, j)): continue
                if stats is not None: stats['triples'] = stats.get('triples', 0) + 1
                q1 = count((i, l), i, l, j); q2 = count((l, j), i, l, j)
                v = min(K, 1 + q1 + q2)
                b = best[v]
                # smaller apex angle at l  <=>  current best strictly inside circle(i,l,j)
                if b is None or incircle(P[i], P[l], P[j], P[b]) > 0: best[v] = l
            Lc = [None]*K; cur = None
            for c in range(K, 0, -1):
                b = best[c]
                if b is not None and (cur is None or incircle(P[i], P[b], P[j], P[cur]) > 0): cur = b
                Lc[c-1] = cur
            L[(i, j)] = Lc
    for (i, j), Lc in L.items():
        if Lc[K-1] is not None:
            return True, reconstruct(P, k, L, i, j, K)
    return False, None

def reconstruct(P, k, L, i, j, K):
    def go(i, j, c):
        if c == 0: return []
        l = L[(i, j)][c-1]
        def cnt(ch):
            Lc = L.get(ch)
            if Lc is None: return 0
            return sum(Lc[x] is not None and incircle(P[i], P[l], P[j], P[Lc[x]]) <= 0 for x in range(K))
        q1, q2 = cnt((i, l)), cnt((l, j))
        c1 = min(q1, c-1); c2 = c-1-c1
        assert c2 <= q2, "witness reconstruction failed"
        return go(i, l, c1) + [l] + go(l, j, c2)
    return sorted([i, j] + go(i, j, K))

def capped_opt(P, k):
    if k == 2:
        m = max(d2(P[a], P[b]) for a, b in itertools.combinations(range(len(P)), 2))
        return m, None
    W = sorted({d2(P[a], P[b]) for a, b in itertools.combinations(range(len(P)), 2)})
    lo, hi, best, wit = 0, len(W)-1, None, None
    while lo <= hi:
        mid = (lo+hi)//2
        ok, w = capped_decision(P, k, W[mid])
        if ok: best, wit, lo = W[mid], w, mid+1
        else: hi = mid-1
    return best, wit

# ---------- Lemma 10 window ----------
def greedy_window(P, k):
    n = len(P); G = [0]
    dist = [d2(P[0], P[t]) for t in range(n)]
    while len(G) < k:
        t = max(range(n), key=lambda x: dist[x]); G.append(t)
        dist = [min(dist[x], d2(P[t], P[x])) for x in range(n)]
    g2 = min(d2(P[a], P[b]) for a, b in itertools.combinations(G, 2))
    per = sum(math.dist(P[t], P[(t+1) % n]) for t in range(n))
    return g2, per

# ---------- greedy structure (Section 6) ----------
def a_of(P, q, r2):
    n = len(P)
    for t in range(1, n):
        if d2(P[q], P[(q+t) % n]) >= r2: return (q+t) % n
def b_of(P, q, r2):
    n = len(P)
    for t in range(1, n):
        if d2(P[q], P[(q-t) % n]) >= r2: return (q-t) % n

def step(P, r2, chosen, fplus, fminus, side):
    """first point beyond the frontier on `side` (+1 ccw / -1 cw), strictly before the other frontier."""
    n = len(P); f, stop = (fplus, fminus) if side > 0 else (fminus, fplus)
    t = (f + side) % n
    while t != stop:
        if all(d2(P[t], P[c]) >= r2 for c in chosen): return t
        t = (t + side) % n
    return None

def seeded_search(P, k, r2, seed, fplus, fminus, mode):
    """mode 'alt': alternate sides (both starting sides tried); 'any': every left/right pattern."""
    def rec(chosen, fp, fm, nxt):
        if len(chosen) >= k: return True
        sides = [nxt] if mode == 'alt' else [+1, -1]
        for s in sides:
            t = step(P, r2, chosen, fp, fm, s)
            if t is None: continue
            if rec(chosen + [t], t if s > 0 else fp, t if s < 0 else fm, -s): return True
        return False
    if len(set(seed)) < len(seed): return False
    return any(rec(list(seed), fplus, fminus, s) for s in (+1, -1))

def V(P, k, r2, seed3, mode):
    for q in range(len(P)):
        a, b = a_of(P, q, r2), b_of(P, q, r2)
        if a is None: continue
        if seed3:
            if a == b or not d2(P[a], P[b]) >= r2: continue
            if seeded_search(P, k, r2, [q, a, b], a, b, mode): return True
        else:
            if seeded_search(P, k, r2, [q, a], a, q, mode): return True
    return False

def C1(P, k, r2, triple=False):
    for S in far_sets(P, k, r2):
        s = set(S)
        for q in S:
            a, b = a_of(P, q, r2), b_of(P, q, r2)
            if a in s and (not triple or b in s): return True
    return False

def G(P, k, r2):
    """RHS of (G_k): some arc A(i)=P[a(p_i), b(p_i)] contains an r-far (k-1)-set."""
    n = len(P)
    for i in range(n):
        a, b = a_of(P, i, r2), b_of(P, i, r2)
        if a is None: continue
        arc = [a]; t = a
        while t != b: t = (t+1) % n; arc.append(t)
        sub = [P[t] for t in arc]
        if len(sub) >= k-1 and far_sets(sub, k-1, r2): return True
    return False

# ---------- random convex instances ----------
def rand_convex(rng, n, kind):
    for _ in range(1000):
        if kind == 'square': pts = [(rng.randint(-10**4, 10**4), rng.randint(-10**4, 10**4)) for _ in range(8*n)]
        elif kind == 'circle':
            R = 10**5; pts = []
            for _ in range(3*n):
                th = rng.random()*2*math.pi; pts.append((round(R*math.cos(th)), round(R*math.sin(th))))
        elif kind == 'thin':
            pts = []
            for _ in range(3*n):
                th = rng.random()*2*math.pi; pts.append((round(10**5*math.cos(th)), round(2500*math.sin(th))))
        elif kind == 'disk':
            pts = []
            while len(pts) < 6*n:
                x, y = rng.randint(-10**4, 10**4), rng.randint(-10**4, 10**4)
                if x*x+y*y <= 10**8: pts.append((x, y))
        H = strict_hull(pts)
        if len(H) >= n:
            idx = sorted(rng.sample(range(len(H)), n))
            P = [H[t] for t in idx]
            assert is_strictly_convex_ccw(P)
            return P
    raise RuntimeError

if __name__ == '__main__':
    rng = random.Random(20261116)
    out = []
    # 1. capped DP (Algorithm 1 + Thm 11) vs brute force, witness validity, Lemma 7 monotonicity, Lemma 10 window
    trials = mism = winfail = 0
    for kind in ('square', 'disk', 'circle', 'thin'):
        for _ in range(60):
            n = rng.randint(5, 12); k = rng.randint(3, min(7, n))
            P = rand_convex(rng, n, kind)
            bf, _ = brute_opt(P, k); cp, wit = capped_opt(P, k)
            trials += 1
            if bf != cp: mism += 1; print('MISMATCH', kind, n, k, P, bf, cp)
            assert len(wit) == k and min(d2(P[a], P[b]) for a, b in itertools.combinations(wit, 2)) >= cp
            g2, per = greedy_window(P, k)
            if not (g2 <= bf <= 4*g2 and math.sqrt(bf) <= per/k + 1e-9): winfail += 1
    out.append(f"capped DP vs brute force: {trials} instances, {mism} mismatches; window violations: {winfail}")
    print(out[-1], flush=True)

    # 2. Table 3 counterexamples
    G4 = [(-976,42),(4,-49),(971,-41),(979,4),(964,29),(398,49)]
    V1 = [(-90482,144),(-45601,-629),(13859,-779),(38501,-677),(97951,32),(57181,535),(4941,792),(-19180,765)]
    V5 = [(-94942,276),(-57395,-1876),(-27771,-2581),(-2908,-2795),(62044,-1720),(87508,-655),(48547,2138),(2056,2796),(-24473,2629),(-70081,1423)]
    for name, P in (('G4', G4), ('V1', V1), ('V5', V5)):
        print(name, 'strictly convex ccw:', is_strictly_convex_ccw(P))
    r2 = 320756
    out.append(f"G4 example: 4-sets r-far={len(far_sets(G4,4,r2))}, RHS of G4 holds={G(G4,4,r2)}, capped decision={capped_decision(G4,4,r2)[0]}")
    r2 = 2556513005
    out.append(f"V1 example: r-far 4-sets={far_sets(V1,4,r2)}, V1={V(V1,4,r2,True,'alt')}, V3={V(V1,4,r2,True,'any')}, "
               f"C1'={C1(V1,4,r2,True)}, C1={C1(V1,4,r2)}, V5={V(V1,4,r2,False,'alt')}, V4={V(V1,4,r2,False,'any')}")
    r2 = 878078401
    out.append(f"V5 example: r-far 6-sets={far_sets(V5,6,r2)}, V5={V(V5,6,r2,False,'alt')}, V4={V(V5,6,r2,False,'any')}, C1={C1(V5,6,r2)}")
    for s in out[-3:]: print(s, flush=True)

    # 3. random tests of the properties at every feasible threshold
    cnt = {'C1': [0, 0], 'V4': [0, 0], 'G': [0, 0]}
    for it in range(400):
        kind = ('square', 'disk', 'circle', 'thin')[it % 4]
        n = rng.randint(6, 11); k = rng.randint(3, min(6, n))
        P = rand_convex(rng, n, kind)
        W = sorted({d2(P[a], P[b]) for a, b in itertools.combinations(range(n), 2)})
        for r2 in W:
            feas = bool(far_sets(P, k, r2))
            if k >= 3:
                cnt['G'][0] += 1; cnt['G'][1] += (feas != G(P, k, r2))
            if not feas: continue
            cnt['C1'][0] += 1; cnt['C1'][1] += (not C1(P, k, r2))
            cnt['V4'][0] += 1; cnt['V4'][1] += (not V(P, k, r2, False, 'any'))
    out.append(f"random property tests (tests, failures): {cnt}")
    print(out[-1])
    open('results.txt', 'w').write('\n'.join(out) + '\n')
