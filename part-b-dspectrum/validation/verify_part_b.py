"""
Network Reliability – Project Part B – validation (not part of the submitted solution)

Independent checks of dspectrum_reliability.py:
  1. anchor() (backward DSS construction) == direct forward destruction, on random permutations
  2. On small random graphs: R from the exact spectrum (all permutations) == brute-force R (all states)
  3. Exact D-spectrum of the dodecahedron (from exact UP-subset counts, Part A validation)
     -> reliability formula identity, and chi-square fit of the three MC spectra (seed 2026)

Run from the part-b-dspectrum folder:
    python validation/verify_part_b.py

Authors: Kobi Loya, Roman Vikhman
"""
import importlib.util
import itertools
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import dspectrum_reliability as ds  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "exact_validation", os.path.join(HERE, "..", "..", "part-a-cmc", "validation", "exact_validation.py"))
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)


def anchor_forward(perm):
    up = [True] * ds.N_EDGES
    for r, e in enumerate(perm, 1):
        up[e] = False
        if not ds.network_is_up(up):
            return r


def check_anchor(trials=20000):
    rng, perm = random.Random(1), list(range(ds.N_EDGES))
    for _ in range(trials):
        rng.shuffle(perm)
        assert ds.anchor(perm) == anchor_forward(perm)
    print(f"1. anchor() == forward destruction on {trials:,} random permutations: OK")


def check_small_graphs(graphs=40):
    rng, tested = random.Random(3), 0
    while tested < graphs:
        nv, m = rng.randint(4, 6), rng.randint(4, 7)
        edges = [tuple(rng.sample(range(1, nv + 1), 2)) for _ in range(m)]
        T = tuple(rng.sample(range(1, nv + 1), rng.randint(2, 3)))
        if not ds.network_is_up([True] * m, edges, nv, T):
            continue
        cnt, total = [0] * m, 0
        for perm in itertools.permutations(range(m)):
            cnt[ds.anchor(list(perm), edges, nv, T) - 1] += 1
            total += 1
        f = [c / total for c in cnt]
        for p in (0.1, 0.37, 0.5, 0.8, 0.95):
            r_bf = sum(math.prod(p if s else 1 - p for s in st)
                       for st in itertools.product([0, 1], repeat=m) if ds.network_is_up(st, edges, nv, T))
            assert abs(ds.reliability_from_spectrum(f, p) - r_bf) < 1e-12
        tested += 1
    print(f"2. R from exact spectrum == brute-force R on {graphs} random graphs x 5 values of p: OK")


def check_dodecahedron():
    assert ev.EDGES == ds.EDGES and tuple(ev.TERMINALS) == ds.TERMINALS
    n, N = ds.N_EDGES, ev.exact_counts()
    F = [1 - N[n - x] / math.comb(n, x) for x in range(n + 1)]
    f_exact = [F[x] - F[x - 1] for x in range(1, n + 1)]
    support = [i for i, v in enumerate(f_exact, 1) if v > 0]
    print(f"3. exact D-spectrum: nonzero for i = {support[0]}..{support[-1]}")
    worst = max(abs(ds.reliability_from_spectrum(f_exact, p) - ev.exact_reliability(p, N)) for p in ds.P_VALUES)
    print(f"   R from exact spectrum vs exact R: max difference {worst:.1e}")
    rng = random.Random(2026)                    # same stream as main(): the spectra of Table 1
    for M in ds.M_VALUES:
        f = ds.d_spectrum(M, rng)
        obs, exp_, o, e = [], [], 0.0, 0.0
        for i in range(n):                       # pool cells with expected count < 5
            o += f[i] * M
            e += f_exact[i] * M
            if e >= 5:
                obs.append(o); exp_.append(e); o = e = 0.0
        obs[-1] += o; exp_[-1] += e
        chi2, dof = sum((a - b) ** 2 / b for a, b in zip(obs, exp_)), len(obs) - 1
        print(f"   M={M:>5}: chi2 = {chi2:5.1f}, dof = {dof} (95% critical value ~{dof + 1.645 * math.sqrt(2 * dof):.1f})")
    print("   exact R(p): " + "  ".join(f"{p}: {ev.exact_reliability(p, N):.6f}" for p in ds.P_VALUES))


if __name__ == "__main__":
    check_anchor()
    check_small_graphs()
    check_dodecahedron()
