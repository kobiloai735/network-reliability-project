"""
Network Reliability – Project Part A – validation (not part of the submitted solution)

Computes the EXACT terminal reliability R(p) of the dodecahedron network,
to check the CMC estimates produced by ../cmc_reliability.py.

Run from the part-a-cmc folder:
    python validation/exact_validation.py

Author: Kobi Loya
"""
import math
import os
import sys
from collections import deque

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from cmc_reliability import EDGES, TERMINALS, P_VALUES  # noqa: E402


# ---------------------------------------------------------------------------
# Exact reliability (for validation only) – frontier-based dynamic
#    programming.  Returns N[k] = number of edge subsets with exactly k UP edges
#    for which the terminals are connected.   R(p) = sum N[k] p^k (1-p)^(m-k)
#
#    Idea: process the edges one by one.  A vertex is on the "frontier" from its
#    first edge until its last edge.  A state remembers only how the frontier
#    vertices are grouped into components and which terminals each group holds.
#    When a group leaves the frontier:
#       holds all terminals  -> success (absorbing state)
#       holds some terminals -> failure (branch discarded)
#       holds no terminals   -> irrelevant, continue
# ---------------------------------------------------------------------------
def bfs_edge_order(edges, start):
    """Reorder edges in BFS order from `start`.  The result does not depend on
    the order, but a BFS order keeps the frontier small: for the dodecahedron
    the peak number of states drops from ~739,000 to ~800."""
    adj = {}
    for i, (u, v) in enumerate(edges):
        adj.setdefault(u, []).append(i)
        adj.setdefault(v, []).append(i)
    seen, order, queue, visited = set(), [], deque([start]), {start}
    while queue:
        x = queue.popleft()
        for i in adj.get(x, []):
            if i not in seen:
                seen.add(i)
                order.append(edges[i])
            u, v = edges[i]
            y = v if u == x else u
            if y not in visited:
                visited.add(y)
                queue.append(y)
    # edges not reachable from `start` (other connected components) keep their order
    order += [e for i, e in enumerate(edges) if i not in seen]
    return order


def exact_counts(edges=EDGES, terminals=TERMINALS):
    terminals = tuple(dict.fromkeys(terminals))      # ignore duplicate terminals
    m = len(edges)
    if len(terminals) == 1:                          # a single terminal is always "connected"
        return [math.comb(m, k) for k in range(m + 1)]
    edges = bfs_edge_order(edges, terminals[0])
    full = (1 << len(terminals)) - 1                 # bitmask: all terminals
    tbit = {t: 1 << i for i, t in enumerate(terminals)}
    last_use = {}
    for i, (u, v) in enumerate(edges):
        last_use[u] = i
        last_use[v] = i

    def zero():
        return [0] * (m + 1)

    def add(table, key, vec, shift):
        tgt = table.setdefault(key, zero())
        for k, c in enumerate(vec):
            if c:
                tgt[k + shift] += c

    def normalize(assign, masks):
        """Canonical state: frontier (vertex, label) pairs + terminal mask per label."""
        relabel, new_assign, new_masks = {}, [], []
        for vtx, lab in sorted(assign.items()):
            if lab not in relabel:
                relabel[lab] = len(relabel)
                new_masks.append(masks[lab])
            new_assign.append((vtx, relabel[lab]))
        return (tuple(new_assign), tuple(new_masks))

    states = {((), ()): zero()}
    states[((), ())][0] = 1
    success = zero()                 # counts of edge subsets already known to be UP
    for i, (u, v) in enumerate(edges):
        # every subset that already succeeded stays UP whatever edge i does
        success = [a + b for a, b in zip(success, [0] + success[:-1])]
        new_states = {}
        for key, vec in states.items():
            for up in (False, True):
                assign = dict(key[0])
                masks = list(key[1])
                for w in (u, v):                    # introduce new vertices
                    if w not in assign:
                        assign[w] = len(masks)
                        masks.append(tbit.get(w, 0))
                if up:                              # merge the two components
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
                    for k, c in enumerate(vec):
                        if c:
                            success[k + int(up)] += c
                else:
                    add(new_states, normalize(assign, masks), vec, int(up))
        states = new_states
    return success


def exact_reliability(p, counts):
    m = len(counts) - 1              # number of edges, inferred from the counts
    return sum(c * p ** k * (1 - p) ** (m - k) for k, c in enumerate(counts))


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def main():
    counts = exact_counts()
    print(f"Terminals: T1={TERMINALS[0]}, T2={TERMINALS[1]}, T3={TERMINALS[2]}")
    print(f"UP subsets: {sum(counts):,} of 2^{len(EDGES)}")
    print(f"{'p':>5} | {'Exact R':>10}")
    print("-" * 19)
    for p in P_VALUES:
        print(f"{p:>5} | {exact_reliability(p, counts):>10.6f}")
    R = exact_reliability(0.9, counts)
    print(f"\nPart E check: exact R(0.9) = {R:.6f}, "
          f"theoretical RE (M=1000) = sqrt((1-R)/(M*R)) = {math.sqrt((1 - R) / (1000 * R)):.5f}")


if __name__ == "__main__":
    main()
