# tai50b, level by level: a rehearsal (bound_lab session `tai50b_size`, started Nov 16 2026 08:30 UTC)

**What this is.** The whole tai50b search run breadth-first from the saved long root, as an exact measurement of
the levels that `../tai50b_size/README.md` samples thinly. It is a rehearsal: a research build (three patch
files, one of them written today to work round a solver fault found this morning, `../tai_solves` section 5), in
one shared 4-core container, on the request of the bound-lab coordinating session. It is not a proof of record:
that would be a clean run under the coordinator, which is the repository owner's to launch. If this run reaches
the end with no open node it is strong evidence that tai50b's optimum is 458,821,517, and nothing more.

**How.** `bfs.sh` (as `../tai_solves/bfs.sh`): at each depth the open nodes are cut into up to 48 state files, each
bounded by one `bnbdrv -r FILE -d D -o OUT` process, 4 at a time; the children left open are merged for the next
depth. Every node starts warm from its parent with its pool.

- Root: `../tai50b_size/runs/root_long/root.state`, 453,552,676 certified (6,470 iterations, 4,341 CPU-s; bounded
  by the earlier build without `wc`; its 50 children are this run's depth 1).
- Binary, one from depth 1 on: `latest/` + `../regions/regions.patch` + `../tai50b_size/est_resume.patch` +
  `pool_retry.patch` (both steps), built outside the repository (`/home/user/qap_retry2/latest/bnbdrv`).
- Settings: `par/A_wc.par` = `../tai50b_size/par/A.par` with `dnn.beta_sched = w=0.36,wc=2.88`.
- Incumbent 458,821,518 (best known + 1), cutoff 458,821,517.05.
- Traces: `runs/dD_cK/nodes.txt` has one line per node (fixes, QPB / DNN / cut bounds, cutoff, outcome, branching
  line, CPU), `dnn.txt` one line per DNN bound, `err.txt` the retries, `cmd.txt` the command and exit code. All
  committed, every level.

**Resuming.** The state files and the parents' DNN states (123 MB each) are not in git and die with the container.
- Container alive: `./bfs.sh runs/open_D.state D 48` for the depth D named under Status (`progress.txt` has the
  last line written). A level that was interrupted is simply run again: delete `runs/dD_c*` first.
- Container gone: `open/open_D.txt` lists every open node (depth, inherited bound, fixes). Each can be restarted
  cold with `bnbdrv -p par/A_wc.par -u 458821518 -F FIXES` (it loses its parent's state and pool, so its first
  solve and phase start from scratch: roughly the cost of a depth-1 node each, several hundred CPU-s), or the
  levels can be rebuilt warm from the root by rerunning `bfs.sh` from depth 1 (the root state itself must first
  be rebuilt: `../tai50b_size/run_root.sh root_long 0`, 72 minutes).

**Stop conditions** (checked by `levels.py` after every level): a node on the best known solution's path fathomed
or bounded above 458,821,517; an incumbent below 458,821,517; a level more than 3x over the high prediction; or
disk or memory running short. The brief first had a wall limit of 12 hours (about 20:30 UTC); the coordinating
session relayed at 17:35 UTC on 2026-11-16 that the repository owner wants the run kept going, so that limit is
dropped (no script enforced it) and the run goes on level by level until the tree is finished or a condition fires.

## Status

- Levels complete: depths 1 to 4; 1,401 nodes below the root, 341,938 DNN iterations, 139,097 CPU-s = 38.6 core-hours (the root: 6,470 iterations, 4,341 CPU-s, run in `../tai50b_size`).
- Open nodes now: 1058 at depth 5 (`open/open_5.txt`).
- Incumbents found: none yet (the search starts from 458,821,518 = best known + 1).
- Path of the best known solution: 4 nodes bounded on it, **0 violations** (fathomed, or a certified bound above 458,821,517). Largest certified bound on it: 458,119,422.9497.

## The levels [m]

| depth | nodes bounded | branched | fathomed: QPB / DNN first solve / cut phase | open nodes carried to the next depth | DNN bounds, phases | DNN iterations | CPU-s | retries: pool, cold | on-path node: short of 458,821,517 | E2 predicted: nodes, CPU-s | E3 (scale 1) predicted: nodes, branched, CPU-s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 50 | 13 | 0 / 7 / 30 | 637 | 50, 43 | 50,769 | 27,731 | 0, 0 | 2,511,495 (B) | 50, 36,321 | 50, 5, 29,416 |
| 2 | 637 | 9 | 0 / 572 / 56 | 432 | 637, 65 | 107,281 | 46,393 | 1, 1 | 1,897,243 (B) | 654, 74,635 | 833, 14, 38,084 |
| 3 | 432 | 6 | 0 / 425 / 1 | 282 | 432, 7 | 29,612 | 10,354 | 1, 1 | 1,537,515 (B) | 141, 7,510 | 672, 11, 34,732 |
| 4 | 282 | 23 | 0 / 205 / 54 | 1058 | 282, 77 | 154,276 | 54,618 | 0, 0 | 702,094 (B) | 135, 46,558 | 517, 43, 70,500 |

The predictions are `../tai50b_size`'s, made without `wc=2.88` and from a 2,400-iteration phase limit under the same root; E3's depth 2 includes the 588 children of the 12 depth-1 near misses it prices separately.

## progress.txt

```
2026-11-16T08:24:54 depth 1: 50 open nodes
2026-11-16T10:24:37 depth 1 done: 637 open at depth 2; store 1.5G; disk free 28G
2026-11-16T10:24:44 depth 2: 637 open nodes
2026-11-16T14:01:59 depth 2 done: 432 open at depth 3; store 973M; disk free 29G
2026-11-16T14:02:02 depth 3: 432 open nodes
2026-11-16T14:47:22 depth 3 done: 282 open at depth 4; store 597M; disk free 29G
2026-11-16T14:47:25 depth 4: 282 open nodes
2026-11-16T18:46:31 depth 4 done: 1058 open at depth 5; store 2.1G; disk free 28G
```
