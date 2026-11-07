# The full run's settings

`settings.sh` is the one place that says what the full wil50 run runs. The
launcher (not built yet), the sampling plan (`../roadmap/tree_size.md`; step 6
of `../roadmap/README.md`) and any check meant to run what the full run will
run should source it rather than copy the switches:

```
. experiments/wil50_feasibility/fullrun/settings.sh
fullrun_check || exit 1      # the files exist and latest/bnbdrv knows every switch
latest/bnbdrv ... -u $FULLRUN_UB -s $FULLRUN_SYM $FULLRUN_INSTANCE
```

Its environment is exactly the `cpcutt` arm of `../pipeline_ab/run_pair.sh`:
the same 28 `QAP_DNN_*`/`DNN_*` variables as `../pipeline_ab/runs_cpt/cpcutt/env.txt`,
checked with `diff` on Nov 7.

| part | setting | why | measured |
|---|---|---|---|
| instance | `qaplib/wil50t.dat`, `latest/SYM/wil50t.sym`, `QAP_DNN_CUT_SIDE=1` | the symmetry acts on the second matrix; a `-F i:j` spec on wil50.dat is `j:i` here | 15 depth-1 orbits against 50 (`../sym/README.md`) |
| cut phase | `QAP_DNN_CUT_CP=100:1000` on the zone-0.005, warm-Y, damped phase | one CP-ADMM solve with in-solve separation instead of Polyak rounds; decided Nov 7 | node 10 on wil50t: 1.760x in DNN CPU on the matched nodes, 2.20x end to end (`../pipeline_ab/runs_cpt/`); node 5 of the tail probe: 1.15x estimated, inside one SE (`../tail_probe/runs_node5_cpt/`) |
| state | `QAP_DNN_SHIP=32` | the parent's DNN state in state files | `../pipeline_ab/README.md` |
| speed | `DNN_R32_NONCERT=1`, `DNN_BETA_SCHED=two,...,w=0.25,wc=0.7`, `DNN_ADMM_CROSS=2`, the four `DNN_CORE_*`/`DNN_NOSYMW` switches | the `bestcore` stack | `../pipeline_ab/README.md`, `../speed/` |

Not settled here, and the launcher's to set:

- **The depth spec (`-x`).** The node-10 subtree A/Bs use `-x 0:4:1:300:1`
  from a depth-3 node, that is the DNN bound at depths 3-7. The sampling
  plan's stage A runs depths 0-2 with the DNN bound from depth 1.
- **The frontier and the cutoff ladder** of the sampling plan.
- **The two CP add-ons** stay off: `QAP_DNN_POOL_CP=1` gave a further 1.05x
  in one same-session pair and `DNN_CP_LAZYCERT=1` about 2%
  (`../roadmap/admm/cpadmm/README.md`, "In the search").
