# Regional equalities in CP-ADMM

Implementation `0fc7508b`, based on `ba47c211`. The warm replay also includes
`f734af31`, which obtains the capacity from the region's membership signature
to handle a repeated index in an imported pool. Enable with
`DNN_CP_REGION_EQ=1`. The default is zero.
The new setting changes the recorded settings hash. The replay harness uses
`--allow-settings-change` to reuse stock checkpoints for both arms.

For a region of size `s` in the reduced problem, the affine block enforces
`sum(p_i) = s` and `sum_{j!=i}(W_ij) = (s-1)*p_i` for each free element `i`.
It follows each half-space pass inside the Y-step's Dykstra sweeps. Each block
projection uses the lifted Frobenius metric. Intersecting the block with the
box and the half-spaces still takes multiple sweeps. The equalities already
hold at full DNN feasibility, so the comparison measures convergence to a
certified target.

The regional overlap Gram splits into a `2R` facility-mean system and an `R`
facility-difference system. Both pseudoinverses are cached for the node.
Moment updates reuse the region overlap kernels. The signed equality
multipliers are included in the normal, RHS and absolute magnitude passed to
the certificate. They are local to the projector. When a child builds its
own projector, region capacities are taken from its remapped free sets.

## Validation

`make test` passed all checks, including 4,536 brute-force cases and seven
exact node-count regressions with the new default-off setting. Brute-force
seeds 2 and 3 passed 4,536 cases each. The AWS coordinator scenario was skipped
because its optional dependencies were absent.

Build `latest/bnbdrv` and `experiments/dnn_cuts/dnncuts` with
`make -j2 -C latest bnb dnncuts` and the same LAPACK configuration used for the
benchmark. Then run:

```sh
python3 experiments/cp_equalities_261117/verify_small.py
```

The dense projection reference uses dependent and zero normals, tests both
region sides, checks the lifted correction and regional moments, and finds a
maximum error of `3.11e-15` in any projected entry. The maximum affine residual
is `1.13e-14`. Both signs of multipliers are tested, with a negative RHS, all
120 permutations and a lazy-certificate snapshot. The 18 reduced instance and
region-side combinations all returned CP and frozen-multiplier bounds below
their enumerated exact optima, with zero failed CP certificates.

## Cold-node diagnostic

The two arms used the same committed build and different values of
`DNN_CP_REGION_EQ`, with one OpenBLAS thread. Both used `tho40c.par` with
`cut.cp_iters = 400`. `bnbdrv` ran `tho40t.dat` at incumbent `240516`, with
fixes `26:10,8:0,2:32,6:39,33:24`, depth limit five, and `SYM/tho40t.sym`.
The first solve completed 200 cold iterations at the same certified bound
`237392.2812` for both arms. The phase used `wc=2.88` in both runs.

| Phase check | Off | On | On minus off |
| ---: | ---: | ---: | ---: |
| 100 | 238461.4102 | 238483.9101 | 22.4999 |
| 200 | 239203.4971 | 239221.0556 | 17.5585 |
| 300 | 239642.0406 | 239666.4357 | 24.3951 |
| 400 | 239890.3934 | 239917.5659 | 27.1725 |
| CP CPU seconds | 51.15 | 49.76 | -1.39 |

The first check at or above `239800` was 360 off and 350 on. Five earlier
targets from `238000` through `239500` were first reached at identical
checks. Both phases stopped at the cap below the cutoff `240514.05`. This
development pair overlapped shared workspace activity, and first-solve CPU
also varied by about 2.4% between arms.

## Small search control

Both nug12 arms found the exact objective `578`. The enabled arm's root bound
at check 400 was `571.3551`, up `0.2309` from `571.1242`. The searches took
24 and 28 DNN nodes off and on respectively, with nine CP phases each. CPU
seconds were 1.42 and 1.91 respectively including search and startup. This
control used `tho40c.par` with a `3:4` grid, depth range `0:11`, initial DNN
give-up set to 40, and the phase cap of 400.

## Warm replays with production stopping

Input state, fp32 Y/Z, recipe and binary hashes are in [corpus.json](corpus.json).
Each pair used the same one-node checkpoint and binary, ran back to back on
a shared host, and kept the production stopping policy.

The tho40t depth-three node `26:10,8:0,2:32` was replayed from the same
production parent checkpoint in both arms, including fp32 Y/Z and the cut
pool. Both retained 1,482 inherited cuts and reached `237871.8951` after
110 first-solve iterations. `tho40c.par` was used with production separation,
penalty and handover, with one OpenBLAS thread. No phase cap was added. Both
warm replays qualified for `wc=2.88`.

| Phase check | Off | On | On minus off |
| ---: | ---: | ---: | ---: |
| 100 | 237926.8980 | 237932.9462 | 6.0482 |
| 200 | 238277.7695 | 238286.4853 | 8.7158 |
| 300 | 238429.5377 | 238431.0371 | 1.4994 |
| 400 | 238493.9883 | 238495.4900 | 1.5017 |
| 600 | 238561.3236 | 238562.8141 | 1.4905 |
| 690, handover | 238577.8109 | 238579.4667 | 1.6558 |
| CP CPU seconds | 109.99 | 114.56 | 4.57 |

Both reached `238000`, `238300`, `238400` and `238500` on identical checks.
The enabled arm reached `238200` and `238550` one check earlier at 170 and
550 iterations. Both handed over at check 690, far below cutoff `240514.05`.

The depth-five node `26:10,8:0,2:32,6:39,33:24` was also replayed with its
production parent state. Both arms retained 1,718 inherited cuts and started
the cut phase at the same certified bound `240121.4553` after 110 first-solve
iterations. The first solve was `392.5947` below the cutoff.

| Phase check | Off | On | On minus off |
| ---: | ---: | ---: | ---: |
| 100 | 240121.4553 | 240121.4553 | 0 |
| 200 | 240225.4824 | 240217.1912 | -8.2912 |
| 300 | 240358.9179 | 240358.7682 | -0.1497 |
| 400 | 240426.3881 | 240430.2830 | 3.8949 |
| 600 | 240495.9725 | 240501.2407 | 5.2682 |
| Final certificate | 240514.1271 | 240514.2926 | 0.1655 |
| CP CPU seconds | 80.07 | 72.91 | -7.16 |

| Certified target | First check off | First check on |
| ---: | ---: | ---: |
| 240200 | 190 | 200 |
| 240300 | 260 | 260 |
| 240400 | 360 | 360 |
| 240450 | 450 | 440 |
| 240475 | 530 | 510 |
| 240500 | 630 | 600 |
| 240510 | 670 | 650 |
| 240514.05 (fathom) | 687 | 663 |

The check histories record the common checkpoints. The fathoming iteration
and final certificate are taken from the phase summary. The enabled arm
fathomed 24 iterations earlier, a 3.5% reduction. The depth-three arm reached
the same handover iteration with 1.66 more bound units. The CPU samples are
indicative under shared load. The tested configuration stays off by default,
and testing stopped at the node gate. The replays cover two warm tho40t
nodes at `wc=2.88`. The effect on other node types, penalties or problem
families is open.

```sh
python3 experiment_run_261117/replay.py run equality-warm-off tho40-d3-hard qap-equalities/latest/bnbdrv
python3 experiment_run_261117/replay.py run equality-warm-on tho40-d3-hard qap-equalities/latest/bnbdrv --env DNN_CP_REGION_EQ=1
python3 experiment_run_261117/replay.py run equality-warm-off tho40-d5-hard qap-equalities/latest/bnbdrv
python3 experiment_run_261117/replay.py run equality-warm-on tho40-d5-hard qap-equalities/latest/bnbdrv --env DNN_CP_REGION_EQ=1
```
