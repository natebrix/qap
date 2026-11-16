"""The depth-1 sample: python3 pick_d1.py ROOT.state OUTPREFIX K N SEED
keeps the child on the path of QAPSOLN/tai50b.sln (certainty stratum) and N-1 of the others drawn at random
without replacement (random.Random(SEED)), dealt round-robin in drawn order into K state files OUTPREFIX{i}.state
(so a stage stopped early has still bounded a random subset). Prints the sample."""
import sys, random
sys.path.insert(0, ".")
import splitstate as S, nodes as N
src, pre, k, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
head, nodes, tail = S.read(src)
def fix(ln):
    p = [int(x) for x in ln.split()[5:5 + 50]]
    f = [i for i, x in enumerate(p) if x >= 0]; return f[0], p[f[0]]
on = [ln for ln in nodes if N.BK[fix(ln)[0]] == fix(ln)[1]]; off = [ln for ln in nodes if ln not in on]
random.Random(seed).shuffle(off)
sample = on + off[:n - len(on)]
for i in range(k): S.write("%s%d.state" % (pre, i), head, sample[i::k], tail)
print("on path:", ["%d:%d" % fix(ln) for ln in on]); print("sample order:", " ".join("%d:%d" % fix(ln) for ln in sample))
print("not sampled:", " ".join("%d:%d" % fix(ln) for ln in off[n - len(on):]))
