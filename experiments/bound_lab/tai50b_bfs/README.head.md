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

**Container restart, 2026-11-16 about 23:28 UTC.** The container was restarted during depth 6; the disk survived,
the processes did not. 44 of that level's 48 chunks had finished and were kept; the 4 unfinished ones
(`d6_c42`, `c45`, `c46`, `c47`, 36 of their 68 nodes already bounded) were redone from their input state files
(`bfs.sh` now keeps a chunk that ended with exit 0 and has its `out.state`). The table counts the redone chunks
once: the work lost in the interrupted processes, about 1.5 hours on 4 cores at most, is not in it.

**Second container restart, 2026-11-17 about 01:12 UTC**, during depth 7: 40 of its 48 chunks had finished and were
kept, the others were redone the same way. Both restarts came 70 to 90 minutes after the session's last turn, so
they look like the container being recycled while the session is idle; from here the session checks in every
40 minutes.

## Stopped, 2026-11-17 01:17 UTC, at the repository owner's request

Not on a stop condition: the owner said the estimate is good enough and the run can stop. **This is a rehearsal
stopped part-way, not a solve of tai50b.**

- **Complete levels: depths 1 to 6** [m]: 3,269 nodes below the root, 206,252 CPU-s = 57.3 core-hours (58.5 with
  the root), against 101.7 (E2) and 137.5 (E3) core-hours predicted for the same six levels. Not counted: the work
  lost in the two container restarts (under 2 core-hours each at most).
- **Depth 7, interrupted** [m]: its row in the table below is partial. 630 of its 704 nodes were bounded when it
  stopped (606 fathomed on the first solve, 8 in the cut phase, 16 branched), 17,168 CPU-s; 41 of its 48 chunks
  had finished. Its on-path node was bounded and branched, 305,475 short of 458,821,517.
- **What is not known** [e]: depths 7 and below. Depth 7 had 704 nodes where the estimate had 129 (E2) and 308
  (E3), and levels 4 to 7 each branched 16 to 23 nodes, so the tree had not started to narrow when the run
  stopped; each such level cost 7 to 15 core-hours. The earlier estimate's range (60 to 250 core-hours) is
  consistent with what was measured, with the measured levels at or below its low side, but the run does not
  show where the tree ends.
- Path of the best known solution through depth 7: 7 nodes, all branched, no certified bound above 458,821,517.
  No incumbent below 458,821,518 found. The pooled-solve fault occurred twice (depths 2 and 3), recovered by
  `pool_retry.patch` both times.
- **To resume** (container alive: `runs/open_7.state`, `runs/in_7_*.state`, `runs/store`, the finished chunks'
  `out.state` are on disk): `rm stop; ./bfs.sh runs/open_7.state 7 48` redoes only the unfinished chunks of
  depth 7 and goes on. Container gone: `open/open_7.txt` lists the 704 open nodes with their fixes (see Resuming
  above).
