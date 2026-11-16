# Laptop runs for the bound lab (Nov 16)

Two scripts for the M2 Pro, run from the repository's root on this branch
(`claude/cool-babbage-ilnogj`, whose `latest/` is `cognac`'s plus commit b86352e, the bound lab's
nine default-off keys; `integration/README.md`). Both build the solver first.

| script | what | time [e] |
|---|---|---|
| `tho40ab_cut.sh` | the tho40ab subtree (`-F 26:10,8:0`) through `qapctl local`, as `tho40b_analysis/giveup_ab/run.sh` does it, with two arms: `cut06` (tho40-b's settings + the grid cut phase at zone 0.06; settings only) and `comb` (+ early entry, hand-over, phase penalty). The reference arm g200 is read from giveup_ab's results (`REF=...`, or `ARMS="g200 cut06 comb"` to run it here) | 10-20 minutes an arm if the subtree responds as the warm subtrees did; g200 about an hour |
| `nug30.sh` | nug30 whole, record recipe (`base`) against `nug30_phase`'s best arm (`g5e`), 8 jobs | 3 minutes + under 2 |

Tested on the cloud box (Nov 16): `tho40ab_cut.sh` on the 37-node subtree `26:10,31:0,34:32,2:7`
(g200 37 nodes, 109 CPU s; cut06 1 node, 65; comb 1 node, 57); `nug30_summary.py` on the cloud
runs of `nug30_phase`. Neither has run on macOS.

What to look for in `tho40ab_cut.sh`'s table: the dives put `cut06` at 5.8x less CPU than tho40-b's
settings over the whole tree and `comb` at about twice that again; the warm depth-3 subtrees of
this subtree gave 7.8x for `comb`. Every arm must end with best objective 240516 and no open node.
