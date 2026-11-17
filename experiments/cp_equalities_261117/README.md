# Regional equalities in CP-ADMM

Implementation `0fc7508b`, based on `ba47c211`. Enable with
`DNN_CP_REGION_EQ=1`. The default is zero.

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
