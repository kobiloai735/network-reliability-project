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
# 4. Experiments
# ---------------------------------------------------------------------------
def main(seed=2026):
    rng = random.Random(seed)
    print(f"Terminals: T1={TERMINALS[0]}, T2={TERMINALS[1]}, T3={TERMINALS[2]}")
    print(f"{'p':>5} | {'M=1000':>8} | {'M=10000':>8}")
    print("-" * 28)
    rows = []
    for p in P_VALUES:
        r1 = cmc_reliability(p, 1000, rng)
        r2 = cmc_reliability(p, 10000, rng)
        rows.append((p, r1, r2))
        print(f"{p:>5} | {r1:>8.4f} | {r2:>8.4f}")

    # Part E: p=0.9, M=1000, 5 repetitions
    M = 1000
    print(f"\nPart E: p=0.9, M={M}, 5 runs")
    runs = [cmc_reliability(0.9, M, rng) for _ in range(5)]
    mean = sum(runs) / 5
    s = math.sqrt(sum((x - mean) ** 2 for x in runs) / 4)    # sample std (n-1)
    print("runs:", runs)
    print(f"mean={mean:.5f}  s={s:.5f}  RE=s/mean={s/mean:.5f}")
    return rows, runs


if __name__ == "__main__":
    main()
