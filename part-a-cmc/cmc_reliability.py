"""
Network Reliability – Project Part A
Crude Monte Carlo (CMC) estimation of terminal reliability
for the dodecahedron network (20 reliable nodes, 30 unreliable edges).

UP criterion: Terminal Connectivity – all terminals are in the same
connected component of the subgraph of UP edges.

Author: Kobi Loya
"""
import random
import math
from itertools import product

# ---------------------------------------------------------------------------
# 1. The network (Figure 4.2 – the dodecahedron network)
#    EDGES[i-1] = (u, v) for edge number i, exactly as numbered in the figure
# ---------------------------------------------------------------------------
EDGES = [
    (1, 2),   # 1
    (2, 3),   # 2
    (3, 4),   # 3
    (4, 5),   # 4
    (5, 1),   # 5
    (6, 7),   # 6
    (7, 8),   # 7
    (8, 9),   # 8
    (9, 10),  # 9
    (10, 11), # 10
    (11, 12), # 11
    (12, 13), # 12
    (13, 14), # 13
    (14, 15), # 14
    (15, 6),  # 15
    (16, 17), # 16
    (17, 18), # 17
    (18, 19), # 18
    (19, 20), # 19
    (20, 16), # 20
    (15, 16), # 21
    (7, 17),  # 22
    (9, 18),  # 23
    (11, 19), # 24
    (13, 20), # 25
    (1, 8),   # 26
    (2, 10),  # 27
    (3, 12),  # 28
    (4, 14),  # 29
    (5, 6),   # 30
]
N_NODES = 20
TERMINALS = (4, 10, 17)          # T1, T2, T3 – as given by the lecturer

P_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
M_VALUES = [1000, 10000]


# ---------------------------------------------------------------------------
# 2. Disjoint Set Structure (Union-Find) with path compression + union by rank
# ---------------------------------------------------------------------------
class DSS:
    def __init__(self, n):
        self.parent = list(range(n + 1))   # nodes are 1..n
        self.rank = [0] * (n + 1)

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:       # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def network_is_up(state, edges=EDGES, n=N_NODES, terminals=TERMINALS):
    """state[i] == True  <=> edge i+1 is UP. Returns True if terminals connected."""
    dss = DSS(n)
    for (u, v), up in zip(edges, state):
        if up:
            dss.union(u, v)
    root = dss.find(terminals[0])
    return all(dss.find(t) == root for t in terminals[1:])


# ---------------------------------------------------------------------------
# 3. Crude Monte Carlo
# ---------------------------------------------------------------------------
def cmc_reliability(p, M, rng, edges=EDGES, n=N_NODES, terminals=TERMINALS):
    """Steps 1-6 of the scheme. p may be a scalar (same for all edges) or a list p_i."""
    probs = p if isinstance(p, (list, tuple)) else [p] * len(edges)
    r = 0                                                   # step 2
    for _ in range(M):                                      # step 5
        state = [rng.random() <= pi for pi in probs]        # step 3: x <= p_i -> up
        if network_is_up(state, edges, n, terminals):       # step 4
            r += 1
    return r / M                                            # step 6: unbiased estimator


# ---------------------------------------------------------------------------
# 4. Exact reliability (for validation only) – frontier-based dynamic
#    programming.  Returns N[k] = number of edge subsets with exactly k UP edges
#    for which the terminals are connected.   R(p) = sum N[k] p^k (1-p)^(m-k)
# ---------------------------------------------------------------------------
def exact_counts(edges=EDGES, terminals=TERMINALS):
    m = len(edges)
    tset = set(terminals)
    full = (1 << len(terminals)) - 1
    tbit = {t: 1 << i for i, t in enumerate(terminals)}
    last_use = {}
    for i, (u, v) in enumerate(edges):
        last_use[u] = i
        last_use[v] = i

    # state: (frontier tuple of (vertex, label) sorted, masks tuple per label) or 'OK'
    def normalize(assign, masks):
        relabel, new_assign, new_masks = {}, [], []
        for vtx, lab in sorted(assign.items()):
            if lab not in relabel:
                relabel[lab] = len(relabel)
                new_masks.append(masks[lab])
            new_assign.append((vtx, relabel[lab]))
        return (tuple(new_assign), tuple(new_masks))

    def zero():
        return [0] * (m + 1)

    states = {((), ()): zero()}
    states[((), ())][0] = 1
    for i, (u, v) in enumerate(edges):
        new_states = {}

        def add(key, vec, shift):
            tgt = new_states.setdefault(key, zero())
            for k, c in enumerate(vec):
                if c:
                    tgt[k + shift] += c

        for key, vec in states.items():
            if key == 'OK':
                add('OK', vec, 0)
                add('OK', vec, 1)
                continue
            for up in (False, True):
                assign = dict(key[0])
                masks = list(key[1])
                for w in (u, v):                    # introduce new vertices
                    if w not in assign:
                        assign[w] = len(masks)
                        masks.append(tbit.get(w, 0))
                if up:
                    a, b = assign[u], assign[v]
                    if a != b:
                        masks[a] |= masks[b]
                        masks[b] = 0
                        for w in assign:
                            if assign[w] == b:
                                assign[w] = a
                # retire vertices whose last edge was this one
                dead, ok = False, False
                for w in (u, v):
                    if w in assign and last_use[w] == i:
                        lab = assign.pop(w)
                        if lab not in assign.values():   # component left frontier
                            if masks[lab] == full:
                                ok = True
                            elif masks[lab] != 0:
                                dead = True
                if dead:
                    continue
                if ok:
                    add('OK', vec, int(up))
                else:
                    add(normalize(assign, masks), vec, int(up))
        states = new_states
    return states.get('OK', zero())


def exact_reliability(p, counts, m=len(EDGES)):
    return sum(c * p ** k * (1 - p) ** (m - k) for k, c in enumerate(counts))


# ---------------------------------------------------------------------------
# 5. Experiments
# ---------------------------------------------------------------------------
def main(seed=2026):
    rng = random.Random(seed)
    counts = exact_counts()
    print(f"Terminals: T1={TERMINALS[0]}, T2={TERMINALS[1]}, T3={TERMINALS[2]}")
    print(f"{'p':>5} | {'M=1000':>8} | {'M=10000':>8} | {'Exact R':>10}")
    print("-" * 42)
    rows = []
    for p in P_VALUES:
        r1 = cmc_reliability(p, 1000, rng)
        r2 = cmc_reliability(p, 10000, rng)
        R = exact_reliability(p, counts)
        rows.append((p, r1, r2, R))
        print(f"{p:>5} | {r1:>8.4f} | {r2:>8.4f} | {R:>10.6f}")

    # Part E: p=0.9, M=1000, 5 repetitions
    print("\nPart E: p=0.9, M=1000, 5 runs")
    runs = [cmc_reliability(0.9, 1000, rng) for _ in range(5)]
    mean = sum(runs) / 5
    s = math.sqrt(sum((x - mean) ** 2 for x in runs) / 4)
    R = exact_reliability(0.9, counts)
    print("runs:", runs)
    print(f"mean={mean:.5f}  s={s:.5f}  RE=s/mean={s/mean:.5f}")
    print(f"theoretical RE = sqrt((1-R)/(M*R)) = {math.sqrt((1-R)/(1000*R)):.5f}")
    print(f"exact R(0.9) = {R:.6f}, |mean-R|/R = {abs(mean-R)/R:.5f}")
    return rows, runs, counts


if __name__ == "__main__":
    main()
