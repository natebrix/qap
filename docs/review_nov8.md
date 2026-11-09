# Code review, November 8 2026

Five independent read-only passes over `latest/` and `bench/` at commit 3297267,
one per file group (dnn.C; dnncut.C; bnb.C + branch.C; linalg.C, qpbnd.C, fw.C,
lap.C, dlap.C, assign.C; driver.C, state.C, util.C, stats.C, estimate.C,
dnntest.C, bench/parallel.py, bench/bench.py), each asked for out-and-out bugs,
tidying, and linear-algebra replacements. The items below are consolidated and
ranked. **[v]** = verified by hand against the code after the review; **[r]** =
reviewer's finding, read but not independently re-derived. Nothing has been
changed; every item here is a proposal.

The headline: **no path was found in the production configuration where the
certified DNN bound claims a value its dual point does not justify, and no
fathoming or symmetry error on the instances run so far.** The wrong-proof
risks found are latent: a malformed symmetry file for nug22, a checkpoint path
that drops `-d` frontier nodes, and user-supplied upper bounds off the parity
grid.

## A. Could make a proof wrong (none affects nug30, tho30, kra30, tho40 runs to date)

| # | where | finding | status |
|---|---|---|---|
| A1 | `latest/SYM/nug22.sym`, `branch.C:179-229` | The file lacks its `|J2| |J3|` size line, so `read_symmetry` reads |J2|=6, |J3|=17 and runs 11 tokens short; the stream fails silently and J3 holds garbage. Below every root child only locations 12-17 are branched. A nug22 run with this file explores a fraction of the tree and can report a wrong optimum. `read_symmetry` validates nothing: no stream check, no range check, no sortedness check, no record of n or of which matrix the file applies to (hence the `wil50t`/`sko42t`/`tho40t` transposed copies). | [v] file is short by 11 tokens; all 21 other files parse with 0 leftover tokens. Fix: add `2 12` to nug22.sym; make the reader fail hard on stream failure, out-of-range entries, non-ascending lists; record n in the file. |
| A2 | `bnb.C:419-433` | The stop path (SIGUSR1 / time limit) re-inserts the `too_deep` nodes set aside by `-d` before `write_state`; the periodic `-k` checkpoint does not. A resume from a periodic checkpoint of a `-d` run loses the frontier. `parallel.py`'s frontier pass does not use `-k`, so dormant today. | [v]. Fix: write `too_deep` at every checkpoint. |
| A3 | `driver.C:218-221`, `bnb.C:478, 1224`, `qap.C:181` | The fathoming tolerance DELTA (1.95 on even-valued instances, 0.95 on integer ones) assumes the incumbent is an attained value. A hand-given `-u` or `-U` value off the grid (e.g. `-u 6125` on nug30) fathoms nodes containing strictly better solutions and reports none. All bench runs pass the known optimum, so unexposed. | [v] logic; [r] exposure. Fix: round a user `-u` up to the grid in the driver, or reject off-grid values. |
| A4 | `bench/bench.py:75-85` | `read_optimum` reads token 0 as n and token 1 as the optimum; ten `QAPSOLN/*.sln` files hold only the optimum on line 1 (nug27, nug28, nug30, tho30, kra30a/b and -old, esc32d, tai25a), so `--ub opt` passes the first permutation entry as the bound, the root is fathomed, and the result is recorded "ok". Not a solver bug, but a false pass in the test harness. | [v] format of the files. Fix: fix the ten files (kra32.sln was fixed Nov 8) and make the reader require two tokens on line 1. |

## B. Bugs that weaken bounds or results without invalidating them

| # | where | finding | status |
|---|---|---|---|
| B1 | `qpbnd.C:339-347` | In `update_STG`, h is computed as the row norms of UᵀX1 but the trace term tr(X1 T1 X1ᵀ) with T1 = Uᵀdiag(t)U (eigenvectors as rows of U) needs the row norms of U·X1ᵀ. The g term for S is right. Every parametric S,T refresh after the first optimises the wrong functional; the bound stays valid (any feasible s,t is valid). A fix changes node counts, so it needs an A/B against the nug30/tho30 trees. | [v] by derivation; reviewer checked numerically. |
| B2 | `driver.C:118, 478`, `assign.C:173` | `-F` is documented to give the subtree root depth 0; `fix()` increments depth per fix and nothing resets it, so depth = number of fixes. `-x`/`-d` count from the true root (this is what the tho40 depth-2 sweep relied on) and `stats->root_bound` is never set in a `-F` run, so the rule-selection gap and the DNN gate are skewed for every node of the subtree. | [v]. Decide the intended semantics, then fix the code or the help text, and set `root_bound` from the first node when it is 0. |
| B3 | `dnn.C:2619, 2637-2639, 2840-2848` | Opt-in only: with `DNN_CERTX_USE=1` and `par.child_bounds`, `child_lb` adds the standard certificate's reduced costs to a polished certificate's value, which the outer-LAP duals do not justify. No caller in the search sets `child_bounds` (only `dnntest -C`). | [v] unreachable from bnb. Fix: exclude child_bounds under CERTX or take both parts from the standard certificate. |
| B4 | `dnn.C:2398-2412, 2441` | Opt-in only (lean CERTX modes): a non-improving fp64 re-projection leaves Zc modified in place while `ci` describes the unmodified multiplier; the pairing is rigorous only to fp64 rounding, which the lambda shift covers in practice. | [r]. Fix: restore Zc when the refinement does not improve. |
| B5 | `qpbnd.C:749-755`, `lap.C:62, 356-360, 920` | The `pbnd` path (console only) sorts n entries of an (n-1)-vector; the int JV `lap()` has an uninitialised `j2`, garbage duals at dim 1, truncation toward zero of negative costs and BIG=1e8. All dead in the branch and bound (`dlap.C` fixed each of these). | [r]. Delete or route through dlap. |
| B6 | `bnb.C:976-978`, `bnb.h:26` | `qs_int` is `int`; `factor*loc_in_row` exceeds 2^31 in dive estimates on nug28-30. Only the `l_level`/`e_level` estimate columns; all node counters are 64-bit. | [v] type. Fix: use `ctr_type`. |
| B7 | `bnb.C:1079, 1224, 478` | Fathom tests use the raw QPB `bound`, not `node.bound` (which includes the inherited and predicted bound); a node whose predicted bound already exceeds the cutoff still runs the DNN and is branched. And FW's early exit uses `inc` where `inc - DELTA` would do. | [r]. Both are lost fathoms / wasted iterations, not wrong answers. |

## C. Robustness and operations

| # | where | finding | status |
|---|---|---|---|
| C1 | `bench/parallel.py:259-273` | The checkpoint barrier: at each interval the coordinator signals every worker, launches nothing until all have exited, then writes. With long tasks (tho40) the running count drains 7, 6, ..., 0 before every checkpoint: 29% CPU utilisation over 5 h. The fix is a barrier-free checkpoint: the pool plus the input node lists of running tasks is a valid state at any instant (a resume redoes at most the in-flight tasks). The same should back the SIGTERM path, which today SIGKILLs workers and points at a stale checkpoint. | [v]. Highest-value operational fix before any multi-day run. |
| C2 | `bench/parallel.py` | The coordinator is single-threaded and does all state I/O inline: every task exit re-parses a task file, every launch re-emits parent DNN states (10-14 MB each at m=40), a checkpoint builds the whole file in memory about three times (9 GB peak for the 3 GB tho40 checkpoint). Dedup of parent states is by object identity, so a state split across task files is written twice. A content-addressed store for parent states (one file per parent, referenced by id) would cut I/O and memory by an order of magnitude. | [r] (consistent with the 3.0 GB checkpoint observed). |
| C3 | `bench/parallel.py:243-251`, `bench/bench.py:254-260` | `bench.py --jobs N` with a timeout SIGKILLs only the coordinator and orphans the workers. Reviewer could not explain from the code why `timeout 13200` did not stop the coordinator on Nov 8 (a SIGTERM to the process does stop it, as the later TaskStop showed): the signal most likely never reached it. Recommendation: a `--wall-limit` inside parallel.py rather than an external `timeout`. | [r]. |
| C4 | `bench/parallel.py:261-264` | The checkpoint SIGUSR1 can hit a worker before `driver.C:404` installs the handler (the steal path waits `--steal-after` seconds; the checkpoint path does not), killing it and aborting the run with "worker failed". Low probability per event, once per checkpoint per fresh worker. | [r]. |
| C5 | `driver.C:146-331` | Option parsing: `i++` without `i < argc` (a trailing `-u` dereferences NULL), `atoi`/`atof` accept garbage (`-u 3O00` becomes 0 and prunes the whole tree with exit 0), unknown options silently ignored, `-F` silently ignored with `-r`, `-t` silently ignored without `-o`. | [r]. |
| C6 | `bench/parallel.py:114`, `state.C:342-344` | A checkpoint is written with the version of the file read, but worker-returned nodes carry version-6 warm lines; resuming a v≤5 checkpoint produces a file bnbdrv rejects. Only across format upgrades. | [r]. |
| C7 | `dnncut.C:1363, 1730, 1372-1376` | Three unreachable-today crashes: NULL `warm` with `first_iters == 0`; `export_pool` with `pool_out` set and `other_orig` NULL; `trim`'s `keep > size` at kcuts > 500,000. One guard each. | [r]. |
| C8 | `bnb.C:39-55` | `read_incumbent_file` adopts the objective even when the permutation is missing or invalid, so the printed permutation can belong to an older objective. Also running workers never learn a better incumbent (the `-U`/SIGUSR2 mechanism exists but parallel.py does not use it); irrelevant when `--ub` is the optimum. | [r]. |

## D. Tidying (highest value first)

- **The DNN_CERTX experiment** (~650 lines in `dnn.C:1462-1790, 2417-2622`) and the superseded switches whose measurements are in (`DNN_NOSYMW`, `DNN_CORE_NOSYMV`, `DNN_CHOLQR_PASSES`, `DNN_R32_NONCERT`, `DNN_GUARD_DIV`, `DNN_CORE_REFRESH_PARTIAL`, the ~35 `getenv` reads in `bnb.C:1080-1180`): fold the record recipe into defaults, retire the losers, move the CERTX code to its own translation unit or keep only mode 32. The recipe today is 20 environment variables that every launcher must reproduce.
- **Env switches read inconsistently** in `dnn.C`: a dozen are latched once per process in function-local statics, others re-read every call; `dnn.h` says "no global state". One per-call struct (as `CertXEnv` and `BetaSched` already are) would settle it.
- **Duplicated blocks**: the Rump constant (`dnn.C:1180-1186` vs `1718-1722`), the Gilmore-Lawler loop three times (`cert_poly`, `gl_full`, `hg_polish`); in `dnncut.C` the 18-line export block (`2031-2051` vs `2177-2199`), the scale loop (`2002-2010` vs `cp_lifted_scale`), regions construction three times, four cut-key encodings, and two parameter structs for the same projector (`DNNCPParams` vs `DNNCutParams`); `lap()` vs `lap0()` in `lap.C`; `_glb0`/`leader_matrix` in `qpbnd.C`; two reductions of differing rigour in `assign.C` (`reduceQAP` with LowerSum vs `reduceC`/`computeShift` with plain doubles).
- **Dead code**: `#ifndef BNB_STANDALONE` paths in `bnb.C` (always defined), the lgb/boost branching block and ~10 unused selection rules in `branch.C` (two reference undeclared variables and no longer compile if enabled), `bnb_trace`/`print_stack`/`pop_alt`, `linalg.C` `warmstart`/`mtvm_mlt0/1`/`permrows0`/`permcols_old`/`rowdot` and the undefined `mvmt_mlt`, `fw.C` `search_dir`/`deltaG`/`deltaG_alt`/`qp_obj` and the never-read `best_X` copies, the UNROLL blocks in both LAPs, `dsyevr` declared and unused, `print_deep_nodes`, `-l` parsed and never read.
- **Comments that contradict the code**: `-F` "depth is 0"; `dnn.C:1384-1391` "only Z is carried over" vs the `ymode` path; `dnn.h:28` "100 m^4 bytes" (about 132 m^4 in subspace mode, 825 MB at m=50); `bnb.C:1070-1072` and `assign.h:90` on `dnn_ext` not being in state files (it is, version 4); `state.C:282` "version 1-3" (version is 6); `dnncut.h:130` `prefix` "(validation)" while it is the default; `RUNNING.md:53` "initial upper bound 100000" (1e100); the flag table lacks `-U` and `-F`; `branch.h`'s symmetry-format block differs from the reader.
- **Silent no-ops**: the CP cut form ignores `kcuts5/kcuts7/kbox/delta/theta/freeze/pgate` without warning; `-x` ignores a sixth field and accepts `dmax < dmin`; `QAP_DNN_CUT_SIDE` is read twice.
- Small leaks and hazards: `util.C:283/295/307` `new char[50]` per call; `util.C:126-130` params-file field overwriting the `-s` pointer; `DNNResult` without a constructor set field by field in several places; `int w` shadowing `Work w` at `dnn.C:2715`; `strtri`'s hidden-length prototype as `int` where every other wrapper uses `size_t`.

## E. Linear algebra and efficiency

Per DNN iteration at K = (m-1)²+1 ≤ 2,402 the cost is three K×K×kb `sgemm` (Chebyshev filter and Rayleigh-Ritz), the `dsyrk` forming R (about K²·nneg), and a full `ssyevd` refresh amortised over 20 iterations; the O(m^5) transforms and the O(m^4) Y/Z steps are an order of magnitude below. Against that:

| # | where | finding | expected gain |
|---|---|---|---|
| E1 | `dnn.C:1104-1117` | `R = W + FFᵀ` is formed by `dsyrk` in fp64 from fp32 eigenvectors; `DNN_R32_NONCERT=1` does it in fp32 on the 19 of 20 non-certificate iterations and is already in the record recipe. Make it the default. | already taken in the recipe; a default, not a new win |
| E2 | `dnn.C:924-959` | The refresh runs `ssyevd` on all K eigenpairs though only nneg + guard low ones and the top one are used; `DNN_CORE_REFRESH_PARTIAL` hand-chains `ssytrd` + `sstedc` (all K tridiagonal vectors) + `sormtr`. `ssyevr`/`sstemr` with RANGE='I' does strictly less. Trade: bitwise reproducibility against `ssyevd`. | 5-10% per iteration [r] |
| E3 | `dnn.C:869-874` | With `tight_lam = 0` the Weyl update materialises W − Wprev, takes its norm and copies, three K² memory passes per iteration; `QAP_DNN_TIGHT_LAM=2` (in the recipe) replaces them. Make it the default; if the Weyl path stays, fuse the three passes. | in the recipe already |
| E4 | `dnn.C:2146-2147, 723-731` | In the default configuration `forward` writes Wf next to W and then marks it stale because W is symmetrised afterwards, so `proj_psd_subspace` converts W to Wf again: one wasted K² fp32 store per iteration. | a few ms per iteration at K≈2,100 [v] |
| E5 | `dnn.C:796-801` | `DNN_EIG_SINGLE`: the fp32 eigenvector matrix is widened to K² doubles only to build F; build F from the fp32 matrix as the subspace path does. | one K² pass, SINGLE mode only [r] |
| E6 | `dnn.C:1164-1188` | Each Rump attempt copies −S into A and factorises; starting a decade below the last successful shift fails first whenever the shift is at the noise floor, costing an extra partial `dpotrf` per check. Try the last successful shift first. | one partial Cholesky per check [r] |
| E7 | `dnn.C:1241-1244` | `cert_poly` scans all n1² entries of L + scale·Zm for max|M| at every check just to set a penalty any valid P would do; use quantities already at hand. | one O(m^4) pass per check [r] |
| E8 | `dnncut.C:1015-1032, 1384-1392, 1545` | `separate()` collects up to ~1M violated triples (45 MB) and fully sorts them to take sep_k; `partial_sort`/`nth_element` and a reused buffer. `rviolation` is called four times per triple on the same six loads. | the visible knob in the CP form [r] |
| E9 | `dnncut.C:209-222, 405-433` | `aggregate()` is O(|S|²·m²) per region and `assemble()` writes triangles cut by cut; the sigma grouping already in `CutProjector` does both in one pass over L. | O(cuts·|S|²) → O(n1²) [r] |
| E10 | `dnncut.C:1068, 1089-1091, 750-761` | Per-check allocations in `dense_term`, `cert_term_at`; `pstep` rebuilds Csig/Dsig every P-step though `hpass` maintains the analogous sums incrementally. | allocation and one O(ns·R·m²) pass per P-step [r] |
| E11 | `fw.C`, `qpbnd2.C` | `fw()` and `qpbnd2()` allocate their workspaces per call; rules 4, 7, 9 call `qpbnd2` once per child (m² per node). `docs/meschach_removal.md` measured `m_get`/`m_free` at 7.3% on nug20. A per-size workspace. Also `best_X` is copied every improving iteration and never read. | several % of QPB time [r] |
| E12 | `fw.C:378-417` | `deltaG_alt2` (dG = A P B) is hand-written scalar loops; a `dgemm` would be 2-4x faster at n = 30-50 but the makefile pins the arithmetic (`-ffp-contract=off`) for identical node counts across machines. Note only. | blocked by reproducibility |
| E13 | `qpbnd.C:219-282` | `init_ST` recomputes the two symmetric eigendecompositions at every node though they depend only on which rows and columns are free; `need_WU` exists for this and is never used. About 1% of a node. | low |

Checked and found sound by the reviewers: the Householder/rotation transform order; all `dsyrk`/`ssyrk`/`strsm`/`sgemm` layout arguments under the row-major-as-column-major convention; LAPACK workspaces and every `info` test; `(size_t)` casts on every index product; `proj_psd` using the exact negative part on certificate iterations so W − R is NSD by construction; the lazy-certificate proxy being an upper bound on the full certificate; the warm-child Z mapping; the lifted L matching the codebase's objective convention; the Hungarian repaired duals being valid under early exit; the side-1 CP projector being side 0's to the bit; every certified cut value being of the form GL(L + Gᵀμ) − hᵀμ with μ ≥ 0; rule 8/9 index layout and loops; warm-start propagation to children; the J1/J2/J3 logic for every well-formed symmetry file; steal and resume paths neither dropping nor duplicating nodes; atomic state writes; 64-bit node counters.

## Suggested order

1. A1 (nug22.sym and a validating reader), A2, A3: an afternoon, and they close every latent wrong-proof path found.
2. C1 (barrier-free checkpoint) and C3 (`--wall-limit`): before any multi-day run.
3. B1 (the QPB h term), as an A/B on nug30 and tho30 trees.
4. Fold the record recipe into defaults and retire the measured-out switches (D, E1, E3); then E4, E6, E7, E8 as cheap per-iteration wins.
5. C2 (parent-state store) before n ≥ 42 is attempted in earnest.
