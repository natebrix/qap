# Which cheap bound predicts the DNN outcome? (nug30, depth 2, Oct 29 2026)

All 225 depth-2 nodes of nug30 (`bnbdrv -p PARAM/nug_deep.par -u 6124 -d 1 -o`),
each bounded by the solver's console (`b`: QPB, parametric QPB, PB, GLB) and by
the certified DNN bound (`dnntest -f -T 1 -H 0 -i 1000`, cutoff 6122.05). 165 of
the 225 are fathomed by the DNN bound within 1,000 iterations, the same count as
in the full nug30 runs. `nug30_depth2.csv` has every value; `analyze.py` (run it
from a directory holding the raw outputs) computed:

| predictor | AUC (fathomed vs missed) | best-threshold accuracy | Spearman with DNN lb |
|---|---|---|---|
| QPB (parametric) | 0.994 | 96.9% | 0.916 |
| QPB | 0.990 | 96.4% | 0.907 |
| lookahead estimate (known before the node is bounded) | 0.988 | 96.4% | 0.897 |
| GLB | 0.969 | 93.8% | 0.884 |
| PB | 0.920 | 86.2% | 0.745 |
| inherited parent bound | 0.594 | 73.3% | 0.258 |

Overlap of fathomed and missed nodes: QPB (parametric) 17 misses above the
lowest fathom and 9 fathoms below the highest miss, all within 5579.6-5603.0;
GLB 50 and 46; PB 58 and 81.

The purely spectral PB, QPB's ancestor, predicts worst; the purely combinatorial
GLB, a far weaker bound, predicts well. QPB, which like the DNN certificate
combines a spectral part with assignment problems, predicts best. One instance,
one depth, sibling nodes: suggestive, not conclusive.
