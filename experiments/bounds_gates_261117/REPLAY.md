# Production warm-state replay

The replay uses `bnbdrv` for both arms. The baseline generates a shared parent state once. `-r` restores the child's assignments, inherited bounds, the parent's Frank–Wolfe starting point, fp32 DNN state, and cut pool. The solver constructs the conditioned child state and runs the first solve and cut phase via its ordinary entry point. The run manifest records the solver commit, executable hash, recipe and data hashes, state and warm blob hashes, environment, process CPU, and wall time. The trace separates first-solve and cut-phase CPU.

The seed history starts at a cold `-F` anchor. It recreates a local path using the committed production recipe and documented branch choices. The resulting warm state is shared across arms and can differ from the corresponding checkpoint in the live search, which starts above that anchor.

## Build

Use the same single-thread BLAS and LAPACK build for seed, baseline, and candidate. In this environment the SciPy wheel supplies the prefixed LP64 library.

```sh
export QAP_EXPERIMENT_SOURCE="$PWD/qap"
export QAP_EXPERIMENT_SEED="$PWD/qap-experiment-seed"
export QAP_EXPERIMENT_WORK="$PWD/experiment_run_261117"
seed_patch="$QAP_EXPERIMENT_SOURCE/experiments/tai60b_lines/force_line.patch"
git -C "$QAP_EXPERIMENT_SOURCE" worktree add -b research/experiment-seed-261117 "$QAP_EXPERIMENT_SEED" 7fdac382
git -C "$QAP_EXPERIMENT_SEED" apply "$seed_patch"
git -C "$QAP_EXPERIMENT_SEED" add latest/bnb.C
git -C "$QAP_EXPERIMENT_SEED" commit -m 'Use historical forced lines to seed shared experiment states'
qap_blas=$(python3 -c 'from pathlib import Path; import scipy; print(next((Path(scipy.__file__).parents[1] / "scipy.libs").glob("libscipy_openblas-*.so")))')
make -C "$QAP_EXPERIMENT_SEED/latest" -j3 DNN_LAPACK_DEFS='-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS' DNN_LAPACK="$qap_blas -Wl,-rpath,${qap_blas%/*}"
```

Seed patch SHA-256 `b3c0e372e2cbbfebf4941b3f2e55f6cba8245aa8f3362dd02890a91d3f6d0efa`. The seed commit for this run is `7fafc560e436f12c56781832de8c039c686707a6`; it changes `latest/bnb.C` only. Rebuild `dnn.o` when replacing an earlier non-OpenBLAS build, since `make` does not track changes to `DNN_LAPACK_DEFS`. Use an unpatched executable for the stock replay.

## Capture chains

Run these commands from the workspace, using the path to `replay.py` in this directory as `qap_replay`.

```sh
python3 "$qap_replay" capture tho40_d2 tho40 26:10,8:0 c32
python3 "$qap_replay" extract tho40_d2 tho40-d3-hard
python3 "$qap_replay" extract tho40_d2 tho40-d3-mid
python3 "$qap_replay" capture tho40_d3_to39 tho40 26:10,8:0,2:32 c39 --parent tho40-d3-hard
python3 "$qap_replay" extract tho40_d3_to39 tho40-d4-hard-parent
python3 "$qap_replay" capture tho40_d4_to24 tho40 26:10,8:0,2:32,6:39 c24 --parent tho40-d4-hard-parent
python3 "$qap_replay" extract tho40_d4_to24 tho40-d5-hard
python3 "$qap_replay" capture tho40_d3_to24 tho40 26:10,8:0,2:32 c24 --parent tho40-d3-hard
python3 "$qap_replay" extract tho40_d3_to24 tho40-d4-subtree
python3 "$qap_replay" capture tho40_d4_to8 tho40 26:10,8:0,2:32,33:24 c8 --parent tho40-d4-subtree
python3 "$qap_replay" extract tho40_d4_to8 tho40-d5-tight
python3 "$qap_replay" capture tai60b_d1 tai60b 44:45 r43
python3 "$qap_replay" extract tai60b_d1 tai60b-d2
python3 "$qap_replay" capture tai60b_d2 tai60b 44:45,43:58 r1 --parent tai60b-d2
python3 "$qap_replay" extract tai60b_d2 tai60b-d3-tight
```

The node trace records when a forced parent fathoms. If a parent fathoms, that line has no children to extract. Keep each `*.state` file and its `warm/` folder together. The comparator checks both hashes before solving.

## Replay

```sh
python3 "$qap_replay" run stock-400 tho40-d3-hard "$QAP_EXPERIMENT_SOURCE/latest/bnbdrv" --cap 400
python3 "$qap_replay" run scaling-400 tho40-d3-hard "$PWD/qap-scaling/latest/bnbdrv" --cap 400 --env DNN_CP_DSCALE=0.25 --env DNN_CP_DSCALE_UNTIL=200
python3 "$qap_replay" run equalities-400 tho40-d3-hard "$PWD/qap-equalities/latest/bnbdrv" --cap 400 --env DNN_CP_REGION_EQ=1
python3 "$qap_replay" run stock-control nug12 "$QAP_EXPERIMENT_SOURCE/latest/bnbdrv" --complete
python3 "${qap_replay%/*}/summarize.py" --json "$QAP_EXPERIMENT_WORK"/runs/{stock-400,scaling-400,equalities-400}/tho40-d3-hard
```

The only policy change for `--cap` is a hard phase iteration cap. Real cutoffs and all earlier stopping rules still apply. The extra terminal certificate can change the last reported bound. Use equally capped repetitions to compare actual CPU at fixed iteration budgets and uncapped repetitions to compare actual CPU to fathom. Omit `--cap` for the complete production phase. Use `--complete` for a subtree gate after the node screen. Node totals come from the run trace so the seed's counters never enter the subtree totals.

The stock recipes are `infra/phase1/recipes/tho40c.par` for the 5×8 grid on side 1 and `tai60b.par` for the 86 overlapping cluster regions on side 0. Both use `w=0.36,wc=2.88`; tho40 retains `w` when the first solve ends within 300 objective units of the cutoff. The cut phase inherits its production launch logic, including tho40 handover and tai60b's 1,000-iteration descendant give-up. The nug12 control applies a 3×4 grid on side 1, starts with incumbent 580, and must close at the known optimum 578.
