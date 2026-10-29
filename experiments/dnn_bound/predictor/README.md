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

## Do QPB-based branching choices serve the DNN bound? (`lines/`)

Five depth-2 nodes the DNN bound misses (its value 6122.0, 6106.6, 6092.5,
6074.5, 6046.1 against the cutoff 6122.05), nine candidate lines each (4 rows,
4 columns, and the most eccentric free location), every child bounded by QPB
(console) and by the DNN bound with the solver's settings (`-H 300`, cold):
1,260 children. `results.txt` has the per-line table.

- At the two near-miss parents every child of every line is fathomed: the line
  choice does not matter for fathoming (DNN CPU per line varies about 1.5x).
- At the far-miss parent (6046.1) lines give 6 to 17 child misses, and QPB
  scores rank them well (Spearman 0.92 for the average child QPB, 0.98 for the
  number of children above the child-level QPB threshold 5597.6).
- In between, the ranking is noisy; the average child QPB, lookahead's score,
  is sometimes backwards. Summed over the five parents, the line each score
  picks has 7 (average), 10 (minimum) and 6 (count above threshold) more child
  misses than the best line.
- The most eccentric location's line had 0, 0, 1, 1 and 6 child misses against
  best possible 0, 0, 0, 0 and 6, and was usually among the cheapest in DNN CPU.
