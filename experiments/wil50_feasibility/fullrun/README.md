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
| cut phase | `cut.cp = 100:1000` on the zone-0.005, warm-Y, damped phase | one CP-ADMM solve with in-solve separation instead of Polyak rounds; decided Nov 7 | node 10 on wil50t: 1.760x in DNN CPU on the matched nodes, 2.20x end to end (`../pipeline_ab/runs_cpt/`); node 5 of the tail probe: 1.15x estimated, inside one SE (`../tail_probe/runs_node5_cpt/`) |
| state | `ship.bits = 32` | the parent's DNN state in state files | `../pipeline_ab/README.md` |
| speed | `dnn.beta_sched = two,...,w=0.25,wc=0.7`, `dnn.cross = 2`, `dnn.lazycert = 1`, `dnn.amx = 1` | the `bestcore` stack (the other core switches are the solver's only behaviour since Nov 9) | `../pipeline_ab/README.md`, `../speed/` |
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
