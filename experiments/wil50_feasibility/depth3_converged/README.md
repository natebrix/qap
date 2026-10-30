# Depth-3 clear misses run to convergence (Oct 30 2026, Apple M1)

Experiment 7.1 of ../NOTES.md: nodes 1, 2, 3 of the depth-3 sample
(../depth3/nodes.txt), `dnntest -f -T 1 -H 0 -i 2500 -v -u 48816`, the three
run concurrently on an Apple M1 (4 P + 4 E cores, Accelerate,
VECLIB_MAXIMUM_THREADS=1). Cutoff 48,814.05.

| node | short at 300 iterations | short at 800 | short at 2,500 | primal - dual at 2,500 | pres at 2,500 |
|---|---|---|---|---|---|
| 1 `34:40,28:22,14:29` | 515 | 419 | 405 | 1.2 | 1.2e-4 |
| 2 `18:43,11:29,46:44` | 322 | 224 | 202 | 1.0 | 1.2e-4 |
| 3 `27:44,5:16,17:9` | 380 | 277 | 242 | 1.9 | 1.1e-4 |

All three hit the 2,500-iteration limit with primal and dual within 2 of
each other, so these are the relaxation values to within a few units. The
gain from 300 to 2,500 iterations is 111, 120 and 139; from 800 to 2,500 it is
14, 22 and 35. So "short at 300 minus about 120" predicts the converged value
within 20 here, and 800 iterations gets within 35.

Applied to the other six clear misses of the sample (short at 300: 364, 449,
318, 353, 302, 304), the nine clear misses converge to roughly 185 to 405
short, median about 240. At the notes' 130 per level, a 240-short depth-3
node fathoms at depth 5 and a 405-short node at depth 6. One of the three
measured (and probably one or two of the nine) is in the depth-6 class.

Timing: 1.53 s per iteration with three processes running, against 0.57 s
for one alone on this machine (the 40-iteration check before the run). That
is a 2.7x per-process penalty from the shared AMX units at m = 47, much
worse than the 1.5x the M2 Pro notes saw at m = 28: at this size an Apple
chip gets little DNN throughput from a second process. 64 minutes wall for
the three.
