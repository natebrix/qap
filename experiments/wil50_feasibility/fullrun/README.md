# The full run's settings

`fullrun.par` is the one place that says what the full wil50 run runs: the
params file `latest/PARAM/nug_deep_xbar9_all.par` (rule 9 through depth 6,
rule 2 below) with a `[settings]` section holding every DNN, cut-phase, pool
and state setting (`RUNNING.md`, "Settings"). `settings.sh` holds what is not a
setting: the instance, its symmetry file, the upper bound, the thread
variables, and the one diagnostic variable. The launcher (not built yet), the
sampling plan (`../roadmap/tree_size.md`; step 6 of `../roadmap/README.md`) and
any check meant to run what the full run will run source it rather than copy
anything:

```
. experiments/wil50_feasibility/fullrun/settings.sh
fullrun_check || exit 1      # the files exist and the binary's dump under fullrun.par is the golden one
latest/bnbdrv -p $FULLRUN_PAR -u $FULLRUN_UB -s $FULLRUN_SYM ... $FULLRUN_INSTANCE
```

The section's values are exactly the `cpcutt` arm of `../pipeline_ab/run_pair.sh`
(its 24 policy exports, `../pipeline_ab/runs_cpt/cpcutt/env.txt`, Nov 7, each
under its key), plus the sampling plan's depth spec. Before the settings
refactor (`docs/settings_plan.md`, Phase 2, Nov 9) they were environment
variables exported here; `bench/envpar.py` converted them, and the new golden
dump equals the old one value for value through the name map (the only
additions are the keys that were flags, `dnn.depths/giveup/warm`, and the
three `branch.*` keys that were compile-time constants).

| part | setting | why | measured |
|---|---|---|---|
| instance | `qaplib/wil50t.dat`, `latest/SYM/wil50t.sym`, `cut.side = 1` | the symmetry acts on the second matrix; a `-F i:j` spec on wil50.dat is `j:i` here | 15 depth-1 orbits against 50 (`../sym/README.md`) |
| depths | `dnn.depths = 1:4`, `dnn.giveup = 300`, `dnn.warm = 1`, `dnn.ext = 20` | the DNN bound from depth 1 (`../roadmap/tree_size.md` section 4) through depth 4, then 20 levels below a miss: from depth 3 down exactly the absolute depths of the measured node-10 arms (`-x 0:4:1:300:1` from a depth-3 node) | stage A ran it (`../sampling/`) |
| cut phase | `cut.cp = 100:1000` on the warm-Y, damped phase; `cut.cp_sweeps = 3:1` (since Nov 13) | one CP-ADMM solve with in-solve separation instead of Polyak rounds; decided Nov 7 | node 10 on wil50t: 1.760x in DNN CPU on the matched nodes, 2.20x end to end (`../pipeline_ab/runs_cpt/`); node 5 of the tail probe: 1.15x estimated, inside one SE (`../tail_probe/runs_node5_cpt/`); three sweeps: node 5's subtree 1.033x in DNN CPU and 1.098x in iterations on the 47 matched children, one more fathom, none lost (`../pipeline_ab/runs_sw2/`; `../review_nov10/README.md` section 2.9) |
| cut zone | `cut.zone = 0.0075` (since Nov 10; 0.005 before) | a node up to 366 short after its first solve runs a phase before it branches (none that far out fathoms by it), so its children inherit a pool and fathom more often in their first solve | the complete subtree of the 325-short depth-3 node `13:1,21:0,45:2`: 1.110x in DNN CPU, 1.149x in iterations, no fathom lost on 2,747 matched nodes (laptop, `../pipeline_ab/runs_zone005/`, `runs_zone0075/`; `../review_nov10/README.md` section 3.5) |
| state | `ship.bits = 32` | the parent's DNN state in state files | `../pipeline_ab/README.md` |
| speed | `dnn.beta_sched = two,...,w=0.25,wc=0.7`, `dnn.cross = 2`, `dnn.lazycert = 1`, `dnn.amx = 1`, `dnn.cross_first = 5` (since Nov 12) | the `bestcore` stack (the other core switches are the solver's only behaviour since Nov 9); `dnn.cross_first` adds one cheap certificate at iteration 5 of a warm call, 1.016x in DNN CPU on the node-10 subtree (`../review_nov10/README.md`, section 2.1, item 5) | `../pipeline_ab/README.md`, `../speed/` |
| pinned | `dnn.beta = 0.3`, `dnn.gu_zone = 0` | 34e347c (Nov 8) changed these defaults; the gated arm and stages A and B ran the old values (`docs/settings_plan.md`, step 0) | |

Not settled here, and the launcher's to set:

- **The frontier and the cutoff ladder** of the sampling plan.
- **The two CP add-ons** stay off: `pool.cp = 1` gave a further 1.05x
  in one same-session pair and `cut.cp_lazycert = 1` about 2%
  (`../roadmap/admm/cpadmm/README.md`, "In the search").

Changing anything in `fullrun.par` is a change of the full run: its own
commit, `settings.golden` regenerated in it
(`(. fullrun/settings.sh; fullrun_dump latest/bnbdrv > fullrun/settings.golden)`),
and stage A's frontier and stage B's class fractions re-examined
(`docs/settings_plan.md`, "Re-gate the full run").

### Nov 10: `cut.zone = 0.0075`, and the golden dump brought up to date

- **The change.** `cut.zone` from 0.005 (244 units) to 0.0075 (366). The
  golden dump changes in four lines: `cut.zone`, the settings hash, the
  params file's hash, and the count, 64 to 65. The fourth is
  `QAP_DNN_WARM_DIR`, a run setting 1d40615 added on Nov 9 and the stale
  golden lacked, so `fullrun_check` failed at HEAD from then until this
  commit. It passes again.
- **Re-gating, not done here.** Stages A and B ran 0.005. Under 0.0075, 2 of
  stage A's 15 depth-1 misses and 92 of its 578 depth-2 misses (244-366
  short) would run a phase. None can fathom by it (the farthest stage A
  phase that fathomed started 170.5 short), and the branching line comes
  from the first solve, so those nodes branch as before. What changes is
  below them: their children inherit a pool. Under the two depth-1 nodes
  that can shrink stage A's frontier, and under the 92 depth-2 nodes it
  moves stage B's children up. On the one subtree measured, nothing was
  lost and the work fell 1.11x; until stages A and B are re-run, their
  estimates overstate the work under these nodes.
- **State files written under 0.005** (stage A's, the depth-1 line test's)
  differ from `fullrun.par` in this setting now, so `-r` refuses them
  without `--allow-settings-change`.
