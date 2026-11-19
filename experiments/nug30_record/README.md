# nug30 record attempt on the laptop (Nov 19 2026): one change at a time against g5e

The record is the Nov 16 run of `g5e` (`../bound_lab/nug30_phase`, `../bound_lab/laptop/README.md`):
1,302 nodes, 722 CPU s, **103 s wall**, 8 jobs on the M2 Pro, from incumbent 6124. The record is
counted in **wall time** and from the **optimum given** (`--ub 6124`), as every nug30 record so far.

`run.sh ARM[@JOBS] ...` runs the whole search with traces on, as the record was run. Each arm is
`par/g5e.par` with one setting changed (checked by diff). The solver is `70eac8d7` with no tracked
change under `latest/`. Each `runs/*/env.txt` has the commit and the load at the run's start.

## Results [m]

All on one evening (18:32-19:03). Both `g5e` runs, at the start and the end, give the same tree.

| run | change from g5e | nodes | DNN iterations | CPU s | wall s |
|---|---|---|---|---|---|
| `g5e` | none (start of the round) | 1,302 | 52,019 | 716 | **101.5** |
| `g5e-2` | none (end of the round) | 1,302 | 52,018 | 716 | **102.5** |
| `g5e-j9` | 9 jobs | 1,302 | 51,996 | 743 | 102.5 |
| `g5e-j10` | 10 jobs | 1,302 | 51,999 | 753 | 100.9 |
| `b3r40` | `dnn.beta = 0.3`, `dnn.refresh = 40` (sko42's) | 1,384 | 55,160 | 747 | 106.9 |
| `gu83-2` | `cut.giveup = 100:8,300:3` | 1,163 | 50,888 | 720 | 104.8 |
| `wcg10` | `cut.cp_wc_gap = 10` | 1,302 | 51,653 | 724 | 103.5 |
| `two` | `dnn.beta_sched = two,f=1.6667,t=60,w=0.36,wc=1.44` | stopped | | | |

Every completed run ends at 6124 with no open node. Run-to-run noise is about 1% of wall time.
During the round, Santa (the security agent) used about 60% of a core and `mediaanalysisd` about 35%.

- **The record reproduces** on today's code: 101.5 and 102.5 s against 103.
- **More jobs do not help.** The 9th and 10th workers run on the efficiency cores. They add CPU and
  give back about 1 s of wall time, which is within the noise. Wall time is already about 88% of
  CPU / 8.
- **sko42's beta and refresh lose 5%.** The `g2s` gain in the bound lab came from sko42's phase
  penalty, which `g5e` already has.
- **A looser stop rule shrinks the tree 11% but costs more than it saves.** It spends 122 s in cut
  phases against 100.
- **`cut.cp_wc_gap = 10` changes almost nothing.** Few phases start within 10 of the cutoff.
- **The two-phase schedule breaks the DNN bound on nug30.** From depth 5 on, 373 first solves
  certify values near -170,000 to -180,000 (`runs/two/dnn.txt`). The values are valid but useless,
  so nothing fathoms there and the tree grows without end. It was stopped after 18 minutes at 2,331
  nodes. The same schedule, with sko42's own `w=0.25,wc=0.7` on beta 0.3, is in `sko42.par`, so it
  should be checked there.
- `runs/gu83` is void. Its solvers were killed by mistake while `two` was being stopped. `gu83-2` is
  the rerun.

## Where the CPU goes in g5e (Nov 16 run, `../bound_lab/laptop/runs_nug30/g5e/dnn.txt`)

| depth | nodes | fathomed | CPU s | of which cut phase |
|---|---|---|---|---|
| 1 | 9 | 0 | 170 | 144 |
| 2 | 225 | 199 | 202 | 54 |
| 3 | 692 | 681 | 244 | 30 |
| 4 | 297 | 294 | 78 | 9 |
| 5 | 78 | 78 | 19 | 6 |

## Conclusion

For the cheap settings tried here, g5e is at a local optimum: no single change gains more than the
noise. A new record needs a structural change that cuts CPU, since wall time is already close to
CPU / 8. The candidates are:

- fewer or cheaper depth-3 first solves, where 692 nodes take 244 s;
- the depth-1 phases, which never fathom but whose pools carry the tree;
- the cost per DNN iteration.

When there is one to try, the next record should be run from 6126, with the optimum found rather
than given.
