# wil50 / sko49: tho40's phase policy on complete subtrees (bound_lab `wil50_policy`, Nov 16 2026)

INTERIM (runs in progress; `python3 summarize.py` prints the table from `runs/`).

- wil50 a220 (`21:22,31:23,40:15`, cold `-F`): base 232 DNN nodes, 17,191 iterations [m `runs/a220_base`];
  `hoe2` = base + `dnn.giveup_warm = 100` + `cut.handover = 32:0:25`: the same 232 nodes, 12,699 iterations,
  **1.35x** [m `runs/a220_hoe2`]; with tho40-like patience far from the cutoff (`32:0:250`, `hoe`) 14,487, 1.19x.
- sko49 s188 (`20:6,40:4,41:16`): base 47 nodes, 4,162 iterations (the documented tree); `hoe` 4,128 (1.01x);
  `wc=1.4` 3,947 (1.05x). No gain there: the subtree is one root phase and 46 children that fathom at once.
- The zone arm is the base arm on these subtrees: every unfathomed node already runs a phase under the base zones.
- `wc` sweep on the far roots: `python3 wcsweep.py`.
