"""
Network Reliability – Project Part B
Reliability estimation by the Destruction Spectrum (D-spectrum),
compared with Crude Monte Carlo (CMC), for the dodecahedron network
(20 reliable nodes, 30 unreliable edges, undirected).

UP criterion: Terminal Connectivity – all terminals are in the same
connected component of the subgraph of UP edges.

Authors: Kobi Loya, Roman Vikhman
"""
import math
import random

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
N_EDGES = len(EDGES)                 # n = 30
TERMINALS = (4, 10, 17)              # T1, T2, T3 – as given by the lecturer

P_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
M_VALUES = [1000, 10000, 50000]


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


def terminals_connected(dss, terminals):
    root = dss.find(terminals[0])
    return all(dss.find(t) == root for t in terminals[1:])


# ---------------------------------------------------------------------------
# 3. Destruction Spectrum (Monte Carlo estimate)
#
#    Edges are destroyed one by one in a random order (a random permutation).
#    The "anchor" r of the permutation is the index of the edge whose failure
#    turns the network from UP to DOWN.  f_i = P(r = i).
#
#    Computing r with DSS: destroying e_1, e_2, ... is the same as building the
#    network backwards, adding e_n, e_(n-1), ... to an empty network.  If the
#    terminals become connected exactly when e_r is added, then with e_r up
#    (and e_(r+1)..e_n up) the network is UP, and without it it is DOWN, so the
#    failure of e_r is the one that makes the network DOWN.
# ---------------------------------------------------------------------------
def anchor(perm, edges=EDGES, n_nodes=N_NODES, terminals=TERMINALS):
    """Index r (1..n) of the edge whose failure turns the network DOWN."""
    dss = DSS(n_nodes)
    for r in range(len(perm), 0, -1):            # add e_n, e_(n-1), ..., e_1
        u, v = edges[perm[r - 1]]
        dss.union(u, v)
        if terminals_connected(dss, terminals):
            return r
    raise ValueError("network is DOWN even when all edges are up")


def d_spectrum(M, rng, edges=EDGES, n_nodes=N_NODES, terminals=TERMINALS):
    """Estimate f_1..f_n from M random permutations. Returns list f (index 0 = f_1)."""
    n = len(edges)
    counts = [0] * n
    perm = list(range(n))
    for _ in range(M):
        rng.shuffle(perm)                          # uniform random permutation
        counts[anchor(perm, edges, n_nodes, terminals) - 1] += 1
    return [c / M for c in counts]


def reliability_from_spectrum(f, p):
    """R(N;p) = 1 - sum_{x=1..n} F(x) * C(n,x) * q^x * p^(n-x),  F(x) = f_1+...+f_x.

    F(x) = fraction of the C(n,x) sets of x failed edges that make the network DOWN.
    """
    n, q = len(f), 1 - p
    F, unrel = 0.0, 0.0
    for x in range(1, n + 1):
        F += f[x - 1]
        unrel += F * math.comb(n, x) * q ** x * p ** (n - x)
    return 1 - unrel


# ---------------------------------------------------------------------------
# 4. Crude Monte Carlo (same scheme as Part A)
# ---------------------------------------------------------------------------
def network_is_up(state, edges=EDGES, n_nodes=N_NODES, terminals=TERMINALS):
    dss = DSS(n_nodes)
    for (u, v), up in zip(edges, state):
        if up:
            dss.union(u, v)
    return terminals_connected(dss, terminals)


def cmc_reliability(p, M, rng, edges=EDGES, n_nodes=N_NODES, terminals=TERMINALS):
    r = 0
    for _ in range(M):
        state = [rng.random() <= p for _ in edges]
        if network_is_up(state, edges, n_nodes, terminals):
            r += 1
    return r / M


def relative_error(sample):
    """r.e. = s / mean, s = sample standard deviation (n-1)."""
    k = len(sample)
    mean = sum(sample) / k
    s = math.sqrt(sum((x - mean) ** 2 for x in sample) / (k - 1))
    return mean, s, s / mean


# ---------------------------------------------------------------------------
# 5. Experiments – sections (1)-(6) of the assignment
# ---------------------------------------------------------------------------
def main(seed=2026):
    rng = random.Random(seed)
    print(f"Terminals: T1={TERMINALS[0]}, T2={TERMINALS[1]}, T3={TERMINALS[2]}")
    print(f"n = {N_EDGES} edges, seed = {seed}")

    # (1) three D-spectra
    spectra = {M: d_spectrum(M, rng) for M in M_VALUES}
    print("\nTable 1 – Destruction Spectra")
    print(f"{'i':>3} | " + " | ".join(f"{'M=' + str(M):>9}" for M in M_VALUES))
    print(f"{'':>3} | " + " | ".join(f"{'f_i':>9}" for _ in M_VALUES))
    print("-" * 39)
    for i in range(1, N_EDGES + 1):
        print(f"{i:>3} | " + " | ".join(f"{spectra[M][i - 1]:>9.5f}" for M in M_VALUES))

    # (2) reliability from each spectrum
    print("\nTable 2 – R(N;p) (Destruction Spectrum)")
    print(f"{'p':>5} | " + " | ".join(f"{'M=' + str(M):>9}" for M in M_VALUES))
    print("-" * 41)
    for p in P_VALUES:
        print(f"{p:>5} | " + " | ".join(f"{reliability_from_spectrum(spectra[M], p):>9.6f}"
                                        for M in M_VALUES))

    # (3) CMC for the same p values
    print("\nTable 3 – R(N;p) (CMC)")
    print(f"{'p':>5} | " + " | ".join(f"{'M=' + str(M):>9}" for M in M_VALUES))
    print("-" * 41)
    for p in P_VALUES:
        print(f"{p:>5} | " + " | ".join(f"{cmc_reliability(p, M, rng):>9.6f}" for M in M_VALUES))

    # (4) 10 D-spectra with M=1000 -> R(0.95);  (5) 10 CMC runs with M=1000 at p=0.95
    M, p = 1000, 0.95
    r_ds = [reliability_from_spectrum(d_spectrum(M, rng), p) for _ in range(10)]
    r_cmc = [cmc_reliability(p, M, rng) for _ in range(10)]
    print(f"\nTable 4 – R(N; p={p}),  M={M}")
    print(f"{'i':>3} | {'Destruction Spectrum':>21} | {'CMC':>9}")
    print("-" * 40)
    for i, (a, b) in enumerate(zip(r_ds, r_cmc), 1):
        print(f"{i:>3} | {a:>21.8f} | {b:>9.4f}")
    m_ds, s_ds, re_ds = relative_error(r_ds)
    m_cmc, s_cmc, re_cmc = relative_error(r_cmc)
    print(f"mean| {m_ds:>21.8f} | {m_cmc:>9.4f}")
    print(f"   s| {s_ds:>21.3e} | {s_cmc:>9.3e}")
    print(f"r.e.| {re_ds:>21.3e} | {re_cmc:>9.3e}")
    return spectra, r_ds, r_cmc


if __name__ == "__main__":
    main()
