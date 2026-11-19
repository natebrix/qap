# tai50b on the M2 Pro: `tai50b_gz0nm.par` (Nov 18-19 2026)

`infra/phase1/recipes/tai50b_gz0nm.par` is `recipes/tai50b.par` (rule 9 at the bottom) with two changes:
`dnn.gu_zone = 0` (was 0.003) and `cut.giveup = 2400:100000,4000:0` (was `2400:0`): a phase still
within 100,000 of the cutoff at 2,400 iterations runs on to 4,000. Started from Nathan's terminal:

    OUT=$PWD/infra/laptop/runs/tai50b-gz0nm PARAMS=infra/phase1/recipes/tai50b_gz0nm.par SLOTS=8 \
      caffeinate -i infra/laptop/solve.sh tai50b

Commit 458f1379, 8 slots, incumbent 458,821,518 (best known + 1). The tree closed at 458,821,517, found
by the search and verified by the coordinator. `results.txt` has the per-level table, `export/` the node
and DNN traces (5,119 and 4,999 lines, 0 units not merged).

## Against the earlier tai50b solves

| run | rule | where | nodes | FW iterations | DNN nodes (fathomed) | cut phases (fathomed) | CPU s | wall |
|---|---|---|---|---|---|---|---|---|
| this | 9 | M2 Pro, 8 slots | **5,085** | **513,074** | 4,999 (4,171) | 828 (679) | **182,500** | **7.6 h** |
| Nov 16 (`bound_lab/laptop/runs_tai50b`) | 2 | M2 Pro, 10 slots | 8,206 | 572,251 | 8,163 (7,592) | 571 (358) | 247,515 | 11.0 h |
| `tai50b-r9` (`infra/phase1/runs/tai50b-r9`) | 9 | c7a | 6,629 | 664,865 | 6,567 (6,085) | 482 (294) | 419,520 | |

- **Nodes and iterations do not depend on the machine.** Against `tai50b-r9`, the same recipe without
  the two changes, the tree is 23% smaller and takes 23% fewer FW iterations.
- **CPU seconds compare only with the Nov 16 run**, the same laptop: 26% less, with 2 fewer slots and a
  different commit (11094c69). c7a seconds are not M2 Pro seconds (the root alone is 4,297 s there,
  1,043 s here).
- **More nodes go to a cut phase and more of those fathom:** 828 phases, 82% fathomed, against 482 and
  61% for `tai50b-r9`, most likely from `gu_zone = 0`. The DNN bound spends 39,500 s, against 174,411 s there.

## The near-miss extension (phases below the root, from `export/dnn.txt.gz`)

| phase ended at | phases | fathomed | CPU s |
|---|---|---|---|
| before 2,400 | 641 | 641 | 67,212 |
| 2,400 (more than 100,000 from the cutoff) | 52 | 1 | 26,721 |
| 2,401-3,999 (extended, then fathomed) | 37 | 37 | 17,548 |
| 4,000 (extended, never fathomed) | 97 | 0 | 28,528 |

134 phases were extended. 37 of them fathomed, each closing a node that would have been branched. The 97
that failed were all still within 100,000 of the cutoff at 4,000 iterations. Their extra 1,600
iterations cost about 11,000 CPU s (6% of the run).

## Not separated

The two changes were made together, so this run does not say how much each one is worth. Telling them
apart needs a run with one of them only.
