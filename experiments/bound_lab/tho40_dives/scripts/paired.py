#!/usr/bin/env python3
"""One cut arm against another on the same dives (python3 scripts/paired.py [--all] [--pair A,B]; default hand,cut06;
arms cut06, hand, comb).

hand_s<c>_a ran with cut06_s<c>_a's seed and hand_s<c>_b with cut06_s<c>_b's, so dive i of one is dive i of the other as
long as both drew the same children (the draws depend on the seed and on the children's QPB scores, not on the DNN bound,
until one arm fathoms a node the other branches). A pair is "same path" when the sampled depth-2, depth-3 and depth-4 nodes
agree over the depths both reached. Per stratum: mean over the pairs in each arm; the ratio of the stratum sums; a bootstrap
over pairs within strata for its 80% band."""
import glob, os, random, re, statistics as st, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import parse, dive_tot, top_tot, HERE, STRATA, fmt
use_all = "--all" in sys.argv
A, B = (sys.argv[sys.argv.index("--pair") + 1].split(",") if "--pair" in sys.argv else ("hand", "cut06"))
NAME = dict(hand="hand", cut06="zone 0.06", comb="combined")
def run(name):
    d = f"{HERE}/runs/{name}"
    if not os.path.exists(d + "/trace.txt"): return None
    done = "exit 0" in open(d + "/cmd.txt").read()
    if not done and not use_all: return None
    top, dives = parse(d + "/trace.txt")
    return top, (dives if done else dives[:-1])
def key(dv):
    k = {}
    for n in dv:
        if n["depth"] <= 4 and n["depth"] not in k: k[n["depth"]] = n["fixes"]
    return k
P = {}
for c in STRATA:
    for b in "abc":
        h, k = run(f"{A}_s{c}_{b}"), run(f"{B}_s{c}_{b}")
        if not h or not k: continue
        e = P.setdefault(c, dict(pairs=[], th=h[0], tk=k[0]))
        for dh, dk in zip(h[1], k[1]):
            kh, kk = key(dh), key(dk)
            e["pairs"].append((dh, dk, all(kh[d] == kk[d] for d in kh if d in kk)))
def tot(sample, keyq):
    H = sum(top_tot(P[c]["th"], keyq) + st.fmean(dive_tot(p[0], keyq) for p in sample[c]) for c in sample)
    K = sum(top_tot(P[c]["tk"], keyq) + st.fmean(dive_tot(p[1], keyq) for p in sample[c]) for c in sample)
    return H, K
S = {c: e["pairs"] for c, e in P.items() if e["pairs"]}
print(f"paired dives: {sum(len(v) for v in S.values())} in {len(S)} strata {sorted(S)}; on the same path: {sum(p[2] for v in S.values() for p in v)}\n")
print(f"## {NAME[A]} against {NAME[B]}\n")
print(f"A = {NAME[A]}, B = {NAME[B]}\n")
print("| stratum | pairs | same path | A: CPU s (mean) | B: CPU s (mean) | CPU A / B | DNN iterations A / B | nodes A / B |")
print("|---|---|---|---|---|---|---|---|")
it = lambda s, i: sum(tot(s, k)[i] for k in ("it1", "itc"))
for c in sorted(S):
    s = {c: S[c]}
    H, K = tot(s, "cpu"); n = tot(s, "nodes")
    print(f"| 26:{c} | {len(S[c])} | {sum(p[2] for p in S[c])} | {fmt(H)} | {fmt(K)} | {H / K:.2f} | {it(s, 0) / it(s, 1):.2f} | {n[0] / n[1]:.2f} |")
if S:
    H, K = tot(S, "cpu"); n = tot(S, "nodes")
    print(f"| all paired strata | {sum(len(v) for v in S.values())} | {sum(p[2] for v in S.values() for p in v)} | {fmt(H)} | {fmt(K)} | {H / K:.2f} | {it(S, 0) / it(S, 1):.2f} | {n[0] / n[1]:.2f} |")
    for label, SS in (("all pairs", S), ("same-path pairs only", {c: [p for p in v if p[2]] for c, v in S.items() if any(p[2] for p in v)})):
        H, K = tot(SS, "cpu")
        rng = random.Random(7); bs = []
        for _ in range(4000):
            s = {c: rng.choices(v, k=len(v)) for c, v in SS.items()}
            h, k = tot(s, "cpu"); bs.append(h / k)
        bs.sort()
        print(f"\nCPU ratio A / B, {label} ({sum(len(v) for v in SS.values())} in {len(SS)} strata): {H / K:.2f}, 80% band (bootstrap over pairs within strata) {bs[400]:.2f} to {bs[3600]:.2f}")
    print("\n| stratum | dive | same path | last depth A / B | CPU s A | CPU s B | ratio |\n|---|---|---|---|---|---|---|")
    for c in sorted(S):
        for i, (dh, dk, same) in enumerate(S[c]):
            a, b = dive_tot(dh, "cpu"), dive_tot(dk, "cpu")
            print(f"| 26:{c} | {i} | {'yes' if same else 'no'} | {max(n['depth'] for n in dh)} / {max(n['depth'] for n in dk)} | {fmt(a)} | {fmt(b)} | {a / b:.2f} |")
