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
12 hours of wall (about 20:30 UTC).

## Status

- Levels complete: depths 1 to 1; 0 nodes below the root, 0 DNN iterations, 0 CPU-s = 0.0 core-hours (the root: 6,470 iterations, 4,341 CPU-s, run in `../tai50b_size`).
- Open nodes now: ? at depth 2 (`open/open_2.txt`).
- Incumbents found: none yet (the search starts from 458,821,518 = best known + 1).
- Path of the best known solution: 0 nodes bounded on it, **0 violations** (fathomed, or a certified bound above 458,821,517).

## The levels [m]

| depth | nodes bounded | branched | fathomed: QPB / DNN first solve / cut phase | open nodes carried to the next depth | DNN bounds, phases | DNN iterations | CPU-s | retries: pool, cold | on-path node: short of 458,821,517 | E2 predicted: nodes, CPU-s | E3 (scale 1) predicted: nodes, branched, CPU-s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 0 / 0 / 0 | ? | 0, 0 | 0 | 0 | 0, 0 | - | 50, 36,321 | 50, 5, 29,416 |

The predictions are `../tai50b_size`'s, made without `wc=2.88` and from a 2,400-iteration phase limit under the same root; E3's depth 2 includes the 588 children of the 12 depth-1 near misses it prices separately.

## progress.txt

```
2026-11-16T08:24:54 depth 1: 50 open nodes
```
