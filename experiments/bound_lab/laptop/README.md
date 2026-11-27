# Laptop runs for the bound lab (Nov 16)

Scripts for the M2 Pro, run from the repository's root on this branch
(`claude/cool-babbage-ilnogj`, whose `latest/` is `cognac`'s plus commit b86352e, the bound lab's
nine default-off keys; `integration/README.md`). Both build the solver first.

| script | what | time [e] |
|---|---|---|
| `tho40ab_cut.sh` | the tho40ab subtree (`-F 26:10,8:0`) through `qapctl local`, as `tho40b_analysis/giveup_ab/run.sh` does it, with two arms: `cut06` (tho40-b's settings + the grid cut phase at zone 0.06; settings only) and `comb` (+ early entry, hand-over, phase penalty). The reference arm g200 is read from giveup_ab's results (`REF=...`, or `ARMS="g200 cut06 comb"` to run it here) | 10-20 minutes an arm if the subtree responds as the warm subtrees did; g200 about an hour |
| `nug30.sh` | nug30 whole, record recipe (`base`) against `nug30_phase`'s best arm (`g5e`), 8 jobs | 3 minutes + under 2 |
| `tai50b.sh` | tai50b whole through `qapctl local` with the cluster-region cut phase, from incumbent best known + 1 (`INST=tai35b` or `tai40b` for the tests). One settings file for the root and the rest: the root's phase runs to 6,000 iterations, phases below stop at 2,400 (no stall rule). It watches the disk and stops the run under 10 GB free | about 100 cloud core-hours by the dives, so 5-10 hours here; unknown beyond that [e] |
| `overnight.sh` | the first two, then `tai50b.sh` from the branch's newest commit | |

Tested on the cloud box (Nov 16): `tho40ab_cut.sh` on the 37-node subtree `26:10,31:0,34:32,2:7`
(g200 37 nodes, 109 CPU s; cut06 1 node, 65; comb 1 node, 57); `nug30_summary.py` on the cloud
runs of `nug30_phase`. Neither has run on macOS.

What to look for in `tho40ab_cut.sh`'s table: the dives put `cut06` at 5.8x less CPU than tho40-b's
settings over the whole tree and `comb` at about twice that again; the warm depth-3 subtrees of
this subtree gave 7.8x for `comb`. Every arm must end with best objective 240516 and no open node.

`tai50b.sh` was tested on the cloud box with `INST=tai35b SLOTS=4` (Nov 16, commit 5f644a13's
solver): 981 nodes, 188,121 DNN iterations, 6,558 CPU s, 37 minutes of wall, 112 units, no retry;
it found 283,315,445 from 283,315,446 and closed the tree. The level-by-level cloud run with the
stall rule took 922 nodes, 125,232 iterations and 5,275 CPU s: without the stall rule 47 unfathomed
phases run to 2,400 iterations, half again the iterations and a quarter more CPU. tai50b itself has
not been run this way anywhere. It is a rehearsal: `latest/` here carries the pooled-solve retry
(`../tai_solves/pool_retry.patch`), which has no key and is not covered by the integration checks.

## Results on the M2 Pro (Nov 16; `runs_tho40ab/`, `runs_nug30/`, pushed from the laptop)

tho40ab subtree `-F 26:10,8:0` through `qapctl local`, 10 slots, solver at 5f644a13 [m]:

| arm | nodes | CPU s | wall s | g200 CPU / arm |
|---|---|---|---|---|
| g200 (tho40-b's settings; `giveup_ab`) | 22,895 | 31,594 | 4,173 | 1 |
| cut06 (cut phase at zone 0.06, settings only) | 1,312 | 9,796 | 1,497 | 3.23 |
| comb (+ early entry, hand-over, phase penalty) | 877 | 3,894 | 650 | 8.11 |

Every arm ends at 240516 with no open node. The warm depth-3 subtrees of this subtree had given
7.8x for comb (`../tho40_search`); the dives give cut06 5.8x over the whole tree (`../tho40_dives`),
and 3.2x on this subtree.

nug30 whole, 8 jobs [m]: base 3,716 nodes, 111,182 iterations, 1,261 CPU s, 179 s; g5e 1,302 nodes,
52,018 iterations, 722 CPU s, **103 s** (1.74x). Both end at 6124.

`tho40ab_cut.sh`'s table step failed on the laptop's Python 3.11 (a backslash inside an f-string);
fixed after that run, and the table above was built by hand from the arms' results.

## tai50b solved on the M2 Pro (Nov 16; `runs_tai50b/`, pushed from the laptop)

`tai50b.sh` (overnight.sh's third step), 10 slots, solver at 11094c69 (5f644a13's `latest/`), from
incumbent 458,821,518 (the best known + 1) [m]:

- **The tree closed at 458,821,517.** 8,206 nodes (8,000 fathomed, 1 enumerated), 1,082 units all
  merged, nothing open. 11.0 hours of wall (02:36 to 13:56 laptop time), 247,515 CPU s (68.8 core-hours,
  slots busy 97.9%). The DNN bound at 8,163 nodes (7,592 fathomed, 89,890 s); the cut phase at 571
  (358 fathomed, 154,893 s).
- **The solution was found by the search** at 13:56:01, ten seconds before the end, by unit
  `u2.1.4.1.4.3.0.3.2.3.2.2.3.0.0.2.0.1.8.14.11.3.8.0.0.2`, and verified by the coordinator
  (`coord/common.py verify_incumbent`: the exact integer objective recomputed in Python from
  `qaplib/tai50b.dat` and the reported permutation). With an incumbent of the best known + 1, nodes on
  its path could not be fathomed (cutoff 458,821,517.05), so the search had to reach its leaf.
- **The levels:** 1 / 50 / 490 / 576 / 282 / 782 / 585 ... nodes at depths 0-6, the deepest node at depth
  47. 40 of the 50 depth-1 nodes fathomed. Against the level-by-level rehearsal (`tai50b_bfs`, its
  phases stopped by the stall rule): 1,398 nodes at depths 1-4 here, 1,401 there; 37 of 50 depth-1 nodes
  fathomed there, the longer phases here fathoming three more.
- **Against the estimates:** the dives' 5,000-13,000 nodes and 60-250 cloud-box core-hours
  (`../tai50b_size`). 68.8 M2 Pro core-hours is about 125 on the cloud box and 170 on c7a.

What it is: a complete branch and bound with certified bounds, on one machine, through the
coordinator's local mode, with the region family `line:1,top:5,linegap:4:0.03`. Not yet a run of
record: the next one is `infra/phase1/start.sh tai50b` on cognac (`docs/next_runs.md`), which carries
the same settings with the region file's hash in its record; and a second route, another region family,
as for tai35b and tai40b.
