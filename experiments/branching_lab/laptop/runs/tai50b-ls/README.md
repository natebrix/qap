# tai50b with the landscape signal choosing lines (mode a), on the M2 Pro: stopped, verdict reached

`../tai50b_ls.sh` as written: worktree at cognac 0e6c2ecf with `experiment_build.patch` and `lsla.patch`,
`QAP_LS=a:2:0:10`, `QAP_LS_SEARCH=30:3000:2.5:300`, `QAP_ROOT_CACHE=off`, 8 slots, incumbent 458,821,518
(best known + 1), everything else `recipes/tai50b.par` (rule 9 at every depth). The tai25b check passed
(5 firings, 26 nodes, root line c7). Started 22:23 Nov 19; stopped by Nathan at 10:01 Nov 20, 11 h 15 min
in, once the decision rule could no longer come out "has legs".

## Verdict

**Fails the decision rule, as far as it was run: the tree is larger than tai50b-r9's whole tree would allow,
and the cause is the signal's root line.**

- **At the stop:** 5,606 nodes (5,022 fathomed), 254,577 CPU s; 47 units queued and 8 leased, holding about
  610 open nodes at depths 3-6. Even if every open node fathomed with no children, the tree would be about
  6,200 nodes. With the children that depth-5 nodes have on tai50b, it would end well above 6,960, the
  "harmful" bound. "Has legs" (at most 6,300) was out of reach once the open nodes were counted. The final
  count is not measured.
- **The root.** The signal chose **c26** (its #1: predicted 12 survivors, score 96.33). Rule 9's line at
  this root is **r31**, as in `cognac infra/laptop/runs/tai50b-gz0nm` (the same root state: root bound
  455,449,196.18 in both). The signal ranked r31 **#23** (predicted 15 survivors).
- **Measured at depth 1:** under c26, **28 of 50** children branched. Under r31, 10 of 50 did in tai50b-r9
  (40 fathomed) and 9 in gz0nm. The signal predicted fewer survivors for c26 than for r31 and got the order
  wrong by nearly 3x.
- **That fixed the tree's width:** depth 2 has 1,372 nodes (28 x 49) against tai50b-r9's 490 (10 x 49). Depths
  3 and 4 are 2,778 and 1,264 against 576 and 282.

| depth | nodes (branched), this run | tai50b-r9 nodes |
|---|---|---|
| 0 | 1 (1) | 1 |
| 1 | 50 (28) | 50 |
| 2 | 1,372 (58) | 490 |
| 3 | 2,778 (27) | 576 |
| 4 | 1,264 (13) | 282 |
| 5 | 199 (4) | 782 |
| 6 | 2 (0) | 630 |

Counts are unique nodes in the merged unit traces and the frontier (`nodes.txt.gz`, 5,666 lines: 131
branched, 5,077 fathomed by the DNN bound, 458 by a cut phase). They include a few nodes from units still
leased at the stop, so they run slightly above the coordinator's 5,606.

## What this does and does not answer

- **Answered:** on tai50b, mode (a) at the root is harmful. Its top line made a tree 2.8x wider at depth 2
  than rule 9's line, which it ranked 23rd. This matches the tai60b depth-9 miss at 24:2 (`../lsla/README.md`:
  c30 predicted 2 survivors, had 12 branching children). That makes two large misses at chain nodes, one of
  them the root.
- **Not answered:** whether the signal's choices at depths 3-9 help. They were made below the wrong root line,
  so this tree cannot isolate them. The landscape README's prediction (the signal agrees with rule 9 at depths
  0-2) did not hold at depth 0 in this build.
- **Cheapest next check (not run):** the same script with the signal kept off the root and depth 1-2 (rule 9
  there), so its first firing is at depth 3. Same decision rule, about 8-12 h on the laptop. Or, cheaper,
  check the c26/r31 inversion against the landscape tables: why the signal's predicted survivors at this root
  (12 vs 15) disagree with the measured 28 vs 10.

## Signal cost

131 firings (depth 0: 1, 1: 28, 2: 58, 3: 27, 4: 13, 5: 4), 34 past the gate, 4,675 CPU s in all (1.8% of
the run; the root's firing 267 s).

## Files

`ls.log` (one line per firing), `run.log`, `coord.log`, `cmd.txt` (the environment), `nodes.txt.gz` (as above).
No `results.txt` or export: the run was stopped before it finished. Its store is on the laptop in
`../qap-lsla/infra/laptop/runs/tai50b-ls/`, so the same command resumes it.
