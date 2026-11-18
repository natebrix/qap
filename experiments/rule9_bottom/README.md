# Rule 9 at the bottom of the tai-b trees: a whole-tree A/B on tai50b (Nov 17 2026)

Question (Nathan): every production recipe branches by rule 9 through depth 6 and by rule 2 below.
Should rule 2 go entirely? `../bound_lab/rule2_ab/` said it made no difference on tho40, but there
depth 8 and below was 1.9% of the subtree's CPU (7.5% of tho40-b's tree at 25 h).

## What the two tai-b trees show [m]

From the node traces: tai60b-261117's live tree on the website at 16:03 UTC Nov 17 (42,305 nodes
merged, run still going) and tai50b-261116's (`infra/phase1/runs/tai50b-261116/nodes.txt.gz`, final).

- **The bottom block decides most of the tree.** Depth 7 and below: tai60b 81% of the merged nodes
  and 75% of their CPU; tai50b 66% of the nodes and 35% of the CPU.
- **Rule 2 picks one location per depth.** It counts the children whose QPB bound plus U_ij reaches
  the cutoff (`select_ind_max_fathoms`, `latest/branch.C`). QPB is 18-27% below the cutoff at every
  branched tai60b node, so no child counts and the tie-break, the largest U column sum, decides:
  - tai60b: location 54 at all 97 depth-7 branchings, location 2 at 243 of 315 at depth 8;
  - tai50b: location 20 at 9 of 10 at depth 7, 10 at 10 of 12 at depth 8, 1 at 29 of 32 at depth
    10, 33 at 16 of 16 at depth 11.

  Rule 9 at depths 4-6 used 6-18 distinct lines a depth. A line fixed by the locations' geometry is
  what the eccentric rule (6) does; on wil50 it was 2-4.7x worse than the marginal line at depth 3
  (`../wil50_feasibility/tuning/lines/`).
- **The children of rule-2 branchings branch more often** (children of a depth's branchings that
  branch in turn):

  | parent depth | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|
  | rule | 9 | 9 | 9 | 2 | 2 | 2 | 2 |
  | tai60b (complete parents only) | 3.6% | 1.9% | 2.2% | **5.9%** | 1.4% | 1.6% | 0.1% |
  | tai50b | 1.8% | 1.7% | 2.1% | **2.8%** | **2.8%** | **5.6%** | 1.2% |

  Suggestive, not proof: the counts are small and depth is confounded with the rule.
- **Rule 9 costs about nothing here**: a QPB lookahead over two lines (tho40's A/B: 745 against 748
  CPU s).

## The A/B

- **Base:** tai50b-261116, the run of record (rule 2 below depth 6): 8,216 nodes, 133 c7a
  core-hours, 4.4 h on one c7a.16xlarge. The laptop's solve under the same settings agreed to 0.1% in
  nodes (`compare.py` below: 8,225 against 8,214 trace lines; up to 10% apart at single depths 6-9,
  where the two BLAS libraries round apart). That is the noise.
- **Arm:** `infra/phase1/recipes/tai50b_r9.par`, `tai50b.par` with the bottom block's rule 2 changed to
  9 (one line), the same regions, incumbent and unit sizing (`recipes/tai50b.sh`). `latest/` has
  changed since the base only in the estimator (8bcd08b), which a search does not run.
- **Launch** (the launcher, once tai60b has ended and freed the Spot quota):

  ```
  infra/phase1/start.sh tai50b
    Run name:          tai50b-r9
    Params file:       infra/phase1/recipes/tai50b_r9.par
    Incumbent:         Enter (458,821,518, the best known + 1)
    Instances:         1, c7a.16xlarge, Spot, 64 slots, time limit 72
    Change the other settings?   y; then Enter until "Public or private status": private
  ```

  Private keeps an A/B arm off the website's catalog; public works too if you want to watch it.
  About 133 core-hours: 4.5 hours, a few dollars at Spot. Then `infra/phase1/collect.sh tai50b-r9`
  and commit `infra/phase1/runs/tai50b-r9/`.
- **Read it:**

  ```
  python3 experiments/rule9_bottom/compare.py base=infra/phase1/runs/tai50b-261116/nodes.txt.gz \
      r9=infra/phase1/runs/tai50b-r9/nodes.txt.gz
  ```

- **Decision:** the arm must close the tree at 458,821,517. Fewer nodes and less CPU by more than 2%:
  rule 9 at every depth in every tai-b recipe, and in tho40c.par, where it is harmless. Within 2%:
  no effect on tai50b, where the bottom is 35% of the CPU; tai60b's bottom is 75%, so then a
  subtree A/B on tai60b's depth-6 nodes before deciding. Worse: keep rule 2.
