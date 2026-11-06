# sko49 or wil50: which is cheaper under today's search? (laptop, Nov 6)

The Nov 2 note put sko49 at 2-14x cheaper than wil50 on two models that
predate every lever measured since. This settles it on measured subtrees:
the matched tail probe (`../wil50_feasibility/tail_probe/`, base against
`full`, the production candidate) on two sko49 depth-3 nodes, compared with
the same probe on wil50 nodes 1 and 5 and the pipeline A/B on wil50 node 10.
Runs in `../wil50_feasibility/tail_probe/runs_sko49_n2/` and `runs_sko49_n4/`
(`PROBE_RUNS=... PROBE_ARGS="--node SPEC --instance qaplib/sko49.dat --ub
23386 --grid 7:7" run_matched.sh base full`; cutoff 23,384.05, 3 slots per
arm, binary 7d1dfc1).

## The two sko49 subtrees

| node (converged short) | arm | est. DNN nodes | est. CPU s (± 1 SE) | per level: fathomed/bounded |
|---|---|---|---|---|
| 2, `7:48,26:41,2:12` (304) | base | 2,762 | 103,704 (± 26,701) | L4 15/46; L5 91/93; L6 6/6 |
| | full | 2,357 | **47,754** (± 10,066) | L4 24/46; L5 64/66; L6 6/6 |
| 4, `11:26,45:36,35:28` (263) | base | 11,112 | 192,981 (± 144,919) | L4 25/46; L5 62/63; L6 2/3; L7 3/3 |
| | full | 1,247 | **46,327** (± 30,137) | L4 34/46; L5 35/36; L6 3/3 |

Matched nodes, base against full: node 2 level 5, 63 shared, 61/61
fathomed, 2.07x in iterations (level 4 0.59x, level 6 2.84x); node 4 level
5, 36 shared, 35/35, 2.86x (level 4 0.64x). Subtree: full 2.17x and 4.17x
less CPU than base. **The levers transfer to sko49 whole**: the same shape
as wil50 (level 4 a wash or a loss, level 5 2-4x), and the phase fathoms
more at depth 4 here (9 and 9 extra fathoms of 46, against 2 on wil50 node
1). Node 4's error bars are wide in both arms: one branched level-5 node of
36 (full) and 63 (base) carries levels 6 and 7.

## Cost per depth-3 node at equal shortfall

Laptop CPU under `full`, by the node's converged shortfall:

| instance | node | short | est. CPU s |
|---|---|---|---|
| wil50 | 10 | 158 | 3,295 (pipeline `full`, `../wil50_feasibility/pipeline_ab/`) |
| wil50 | 5 | 316 | 24,204 |
| wil50 | 1 | 405 | 353,678 |
| sko49 | 4 | 263 | 46,327 (10-90% about 16,000-76,000) |
| sko49 | 2 | 304 | 47,754 |

Interpolating wil50's three points in log(cost), sko49's nodes cost **2.3x
(node 2) and 3.7x (node 4; 1.3x at its lower error bar) a wil50 node of the
same shortfall**, geometric mean 2.9x. A sko49 unit is worth more than a
wil50 unit (135 per level against 130, on a cutoff half the size), so equal
shortfall is about equal depth-to-fathom; the extra cost is in wider deep
levels (node 2 has 1,320 estimated level-6 nodes against wil50 node 5's
690 at a similar shortfall) and longer deep solves.

## Totals (`compare_wil50.py`)

Each instance's depth-3 population (29,580 wil50 orbits, 13,923 sko49)
times the mean subtree cost over its depth-3 sample (wil50: 16 nodes of
`../wil50_feasibility/NOTES.md` section 2 with the converged shortfalls of
`stronger_bounds/README.md` step 2; sko49: the 12 of `NOTES.md` section
2), cost from the wil50 curve, sko49 at the 2.9x multiplier; near misses
(up to 120 short) at 700 s, depth-3 fathoms at 40 s. Bootstrap over the
sample for the range.

| | wil50 | sko49 |
|---|---|---|
| depth-3 nodes after symmetry | 29,580 | 13,923 |
| depth-3 sample | 2 fathomed, 5 near misses, 9 clear misses to 405 | 12 clear misses, 134-337 |
| total laptop CPU, `full` | **7.8e8 s** | **6.4e8 s** (3.8e8 with node 4 at its lower bar) |
| on 1,000 Xeon cores (2.0x the laptop per core) | **about 18 days** | **about 15 days** (9-15) |
| bootstrap 10-90% over the sample | 1.1e8 to 1.5e9 (2.5 to 35 days) | 4.8e8 to 8.0e8 (11 to 19 days) |
| share of the total in nodes over 300 short | 89%, from 2 of 16 nodes | 48% |

The wil50 figure agrees with the roadmap's independent 21-30-day center
(`../wil50_feasibility/roadmap/README.md` section 3), which is some
comfort about the method; its node-5 probe argued that center is
conservative, and this model says the same (18).

## Verdict

- **At the center the two cost the same within the error of either:
  about 15 against 18 days on 1,000 cores.** sko49 is not the 2-14x cheaper
  instance the Nov 2 models suggested, because its nodes cost 2-4x more
  at equal shortfall and that cancels its 2.1x narrower tree and the
  absence of a 400-short class.
- **sko49's estimate is four times tighter** (11-19 days against 2.5-35):
  wil50's total hangs on 2 of 16 sampled nodes, sko49's on none in
  particular. If the question is which total we could promise, sko49; if
  it is which is cheaper, neither.
- **wil50 keeps the tooling advantage**: symmetry file in production, the
  node-10 pipeline, two probed tail nodes, every lever measured there. sko49
  would need its own pipeline subtree and sym-file check (`sym/`) before
  any of this could be run as a search; a day or two, not a blocker.
- **Recommendation: stay on wil50, keep sko49 as the fallback with a
  trigger.** The wil50 sampling plan (first hours on the cores) measures
  the 400-short class directly; if it comes in above about 1.5e9 laptop
  CPU s (the top of this model's range), sko49 at 6e8 is the cheaper
  instance and the switch costs a day of setup. Below that, wil50 wins or
  ties, with everything already built.
