# ucheap/suggest: parameters for the retries (Nov 23)

Evidence and its limits: `../README.md`, "sko42 and tai35b".

- **nug30** (the laptop record): in cognac, `experiments/nug30_record/README.md`, "Nov 23":
  `caffeinate -i experiments/nug30_record/run.sh g5e g5e_r11x2 g5e_r11x2_gu83`.
- **sko42**: `sko42_r10.par`. `experiments/estimates/sko42/sko42.par` with rule 10 on every line in every block
  but the first, `branch.la_exp = 32`. The first block is the root and the 12 depth-1 nodes on sko42, left at
  rule 9: the run has sko42-261120b's 438 depth-2 nodes and differs from there down. Over eight subtrees of that
  run: 0.61 of its nodes, 0.78 of its CPU.
  `sko42_r10_sb.par`: the same with `sb.depths = 0:1`, `sb.all = 0`, `dnn.depths = 0:40`, as
  `infra/phase1/recipes/tho40c_sb.par`. Not measured on sko42.
- **tai50b**: no change. `infra/phase1/recipes/tai50b.par` as it is in cognac (rule 9 at every depth). On tai35b
  rule 11 gave 1.06 of rule 9's nodes, rule 10 1.10 (1.33 of the CPU), DNN strong branching at the root and
  depth 1 0.96 (1.81 of the CPU).

To run one on the laptop from a checkout of cognac at 626230c8 or later, with the file copied in:
`PARAMS=<path from the repository root> caffeinate -i infra/laptop/solve.sh sko42`.
