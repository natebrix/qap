# ucheap/suggest: parameters for the retries (Nov 23)

Evidence and its limits: `../README.md`, "sko42 and tai35b".

- **nug30** (the laptop record): `g5e_r11x2_gu83_rt` on cognac, 82.8 s against 103 (f3666c27,
  `experiments/nug30_record/README.md` there). `g5e_r11x2` failed on the laptop for want of `dnn.retry`, not for
  its lines: one pooled solve failed there and its subtree inherited the state.
- **sko42**: **`sko42_d2r10.par`**.
  `experiments/estimates/sko42/sko42.par` with rule 10 on every line at depth 2 only (five blocks: the second
  and third are depth 2, 25 lookahead iterations as measured), `branch.la_exp = 32`. The first block is the root
  and the 12 depth-1 nodes on sko42, left at rule 9: the run has sko42-261120b's 438 depth-2 nodes and differs
  from the depth-2 lines down. Over five depth-2 subtrees of that run: 0.51 of its nodes, 0.70 of its CPU; from
  depth 3 down rule 9 was the best of the three rules (rules 10 and 11 at two lines: 1.04 of its CPU).
  **`sko42_d2r10_sb.par`** (Nathan's choice, Nov 23; on cognac since 91cf6dbe as
  `infra/phase1/recipes/sko42_d2r10_sb.par`, sko42's default in `recipes/defaults.txt`): the same with DNN strong branching at the root and depth
  1 (`sb.depths = 0:1`, `sb.all = 0`, `dnn.depths = 0:40`, as `infra/phase1/recipes/tho40c_sb.par`). Unmeasured as
  a whole. At the 12 depth-1 nodes of sko42-261120b strong branching took rule 9's line at 11, for 18 CPU h of
  firings (`experiments/sb_depth1` in cognac); the root's firing run here
  (`../sko42root/`, Nov 23) **chose r22, not the run's r36**: score 1.057, then r35 1.233, r1 1.250, r5 1.302;
  r36 1.779, tenth of 42 rows; 504 child solves, none retried, 12,688 CPU s. So a run with this file shares no
  node below the root with sko42-261120b: the depth-2 evidence above is from another tree, its 18 CPU h of
  depth-1 firings are at other nodes, and its total cannot be split between the root's row and rule 10. No tree
  under r22 has been measured.
  `sko42_r10.par`: rule 10 in every block but the first (0.61 / 0.78 over eight subtrees; its depth-2 lookahead
  is 12 iterations at most nodes, measured with 25). `sko42_r10_sb.par`: that with `sb.depths = 0:1`,
  `sb.all = 0`, `dnn.depths = 0:40`, as `infra/phase1/recipes/tho40c_sb.par`. Not measured on sko42.
- **tai50b**: no change. `infra/phase1/recipes/tai50b.par` as it is in cognac (rule 9 at every depth). On tai35b
  rule 11 gave 1.06 of rule 9's nodes, rule 10 1.10 (1.33 of the CPU), DNN strong branching at the root and
  depth 1 0.96 (1.81 of the CPU).

On cognac at 91cf6dbe or later `infra/phase1/start.sh sko42` and `caffeinate -i infra/laptop/solve.sh sko42` take
`sko42_d2r10_sb.par` by default. Another of these files, copied in:
`PARAMS=<path from the repository root> caffeinate -i infra/laptop/solve.sh sko42`.

**wil50** (Nov 24): its next-run file is cognac's `infra/phase1/recipes/wil50_d2r10.par` (another session's, Nov 23;
rule 10 at depth 2 as here, no strong branching). `wil50-est.spec` here is a draft of an overnight estimate under
it, for `infra/phase1/specs/` on cognac once Nathan says so: its keys pass `start.sh --spec`'s check; the dry run
itself needs the AWS CLI, which this container lacks.
On cognac since 1a7c43bf (`infra/phase1/specs/wil50-est.spec`), with `sko42_d2r10_sb.par`'s first block through
depth 1 (it went through depth 2, so a depth-2 node more than 0.20 of the root's gap short would have kept
rule 9). `experiments/wil50_pair` (in progress here, for cognac): rule 10 against rule 9 at wil50 depth-2 nodes
on the laptop.
