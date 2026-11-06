"""wil50 vs sko49: total laptop CPU under the production candidate (full), from the
measured depth-3 subtree costs as a function of converged shortfall and each
instance's depth-3 shortfall sample. Interpolation is piecewise linear in log(cost)."""
import math, random, sys
# measured subtree CPU (laptop s, full arm) by converged shortfall
wil = {158: 3295, 316: 24204, 405: 353678}          # node 10 (pipeline), node 5, node 1
sko = {304: 47754}
if len(sys.argv) > 1: sko[int(sys.argv[1])] = float(sys.argv[2])   # node 4 when done
near_miss, fathomed = 700.0, 40.0                     # <=120 short: children fathom at depth 4; fathomed at depth 3
def curve(pts):
    xs = sorted(pts); 
    def f(s):
        if s <= 0: return fathomed
        if s <= 120: return near_miss
        if s <= xs[0]:   # between the near-miss plateau and the first point
            t = (s - 120) / (xs[0] - 120); return math.exp(math.log(near_miss) * (1 - t) + math.log(pts[xs[0]]) * t)
        for a, b in zip(xs, xs[1:]):
            if s <= b:
                t = (s - a) / (b - a); return math.exp(math.log(pts[a]) * (1 - t) + math.log(pts[b]) * t)
        # beyond the last point: extrapolate with the last slope
        a, b = xs[-2], xs[-1]; k = (math.log(pts[b]) - math.log(pts[a])) / (b - a)
        return pts[b] * math.exp(k * (s - b))
    return f
fw = curve(wil)
# sko49 multiplier at equal shortfall, from its measured points against the wil50 curve
mult = [sko[s] / fw(s) for s in sko]; m_sko = math.exp(sum(math.log(x) for x in mult) / len(mult))
wil_sample = [-1, -2, 404.8, 202.2, 240.7, 228.0, 315.7, 56.0, 202.6, 216.0, 43.6, 158.2, 181.7, 84.1, 60.6, 55.6]
sko_sample = [134, 175, 210, 258, 259, 263, 278, 280, 300, 304, 314, 337]
N_wil, N_sko = 29580, 13923
def total(sample, N, mult=1.0):
    return N * sum(fw(s) * mult for s in sample) / len(sample)
def boot(sample, N, mult, B=4000, seed=1):
    rng = random.Random(seed); v = sorted(total([rng.choice(sample) for _ in sample], N, mult) for _ in range(B))
    return v[int(0.1 * B)], v[int(0.5 * B)], v[int(0.9 * B)]
print("sko49 multiplier at equal shortfall:", " ".join(f"{s}:{sko[s]/fw(s):.2f}x" for s in sorted(sko)), f"geometric mean {m_sko:.2f}x")
for name, sample, N, mult in (("wil50", wil_sample, N_wil, 1.0), ("sko49", sko_sample, N_sko, m_sko), ("sko49 (no multiplier)", sko_sample, N_sko, 1.0)):
    T = total(sample, N, mult); lo, med, hi = boot(sample, N, mult)
    print(f"{name:22s} total {T:.3g} laptop CPU s = {T/86400/1000:.0f} laptop-core days / 1,000 cores (~{T*2/86400/1000:.0f} Xeon); "
          f"bootstrap 10-90%: {lo:.2g}..{hi:.2g}; share from nodes >300 short: {N*sum(fw(s)*mult for s in sample if s>300)/len(sample)/T:.0%}")
print("per-node cost at shortfall:", " ".join(f"{s}:{fw(s):.3g}" for s in (60, 134, 175, 210, 260, 300, 316, 337, 405)))
