# Cut-phase trims (Nov 3, cloud): give-up deferral, primal gate, zone

Item 3 of `../roadmap/README.md` (tree memo items 4 and 5). Both rules were
replayed on logs we already had before anything ran; one survived.

## 1. Give-up deferral: killed by the replay

The idea (`QAP_DNN_GU_ZONE`, Goux's D6): a plain-DNN give-up decided near the
cutoff is deferred so that a slow uncut fathom does not pay for a cut phase.
`replay_giveup.py` re-runs the search's extrapolation rule (`latest/dnn.C`,
`hopeless_after` 300 as in `-x 0:4:1:300:1`, a check every 10 iterations)
on the certified traces of the 16 depth-4 children
(`stronger_bounds/depth4/uncut_*.log`, 1,500 iterations with no give-up),
with and without the deferral, charging a give-up the phase that actually
ran on that child (`cuts_*.log`).

| policy | CPU over the 16 children |
|---|---|
| the rule as it runs today (no zone) | 22,572 s |
| zone 0.003, deferral 100 iterations (the default) | 22,967 s |
| deferral 200 | 23,250 s |
| deferral 300 | 23,585 s |
| deferral 500 | 24,219 s |

The three slow uncut fathoms (nodes 3, 5, 13 at 510, 600 and 630
iterations) **fathom under the rule as it stands**: the extrapolation never
fires on them, because their bound keeps rising geometrically. The rule
fires on the nine near misses between 300 and 780 iterations, and every one
of those needs the phase (two fathom with it, seven do not). So the
deferral has nothing to save and costs 100-500 iterations on each miss: 2-7%
on the sample. `QAP_DNN_GU_ZONE` stays off. (The memo's premise came from
the depth-4 table, where the phase was started from the give-up point at
400 iterations by construction, not by the rule.)

## 2. Primal gate: passes the replay, in the library and the search

The idea: the certified bound never passes the Lagrangian primal value, so
a phase whose best primal after a few rounds is well below the cutoff
cannot fathom soon. `replay_pgate.py R M` applies "after round R, stop when
the best primal of rounds 1..R is more than M below the cutoff" to every
logged phase with a `FINAL` line: 88 phases on wil50 nodes from
`stronger_bounds/{laptop/dual,dual,local7,dualtest,depth4}` and the
30-run budget study, 25 of which fathomed (in any round or the quiet round).

| R | M | false stops (would have fathomed) | misses stopped (of 63) | closest fathom / closest caught miss (primal short by) |
|---|---|---|---|---|
| 2 | 15 | 8 | 51 | |
| 2 | 30 | 1 | 46 | |
| 3 | 15 | 1 (`dual/five_14`, a 5-point variant) | 46 | 23.0 / 25.3 |
| 3 | 25 | 0 | 46 | 23.0 / 25.3 |
| **4** | **20** | **0** | **44** | **12.4 / 26.6** |
| 4 | 30 | 0 | 42 | |

Round 4 with margin 20 is the choice: no false stop, 44 of 63 misses
stopped, and 14 units between the closest fathom (budget study `n14_Y12`,
12.4 short after round 4) and the closest caught miss (`laptop/dual/stepA_14`,
26.6). Node 14 is the borderline case in every variant: a near miss that
fathoms only through the quiet round from the best multipliers. Round 3
with margin 25 catches two more misses but sits 2 units from a false stop.

Implementation: `DNNCutParams::pgate_rounds / pgate_margin`
(`latest/dnncut.h`), checked after the stall rule in `dnn_cut_rounds`; a
gated phase skips the quiet round and reports `gated = 1`. Switches:
`QAP_DNN_CUT_PGATE=R:M` in the search (`latest/bnb.C`), `-pgate R:M` in
`dnncuts`. Default off. `smoke_pgate_n13.log`: the gate firing on node 13
with short rounds.

What it is worth: a stopped phase costs 4 rounds instead of 12 (plus no
quiet round): about 1,200 of a failed phase's roughly 2,000 iterations. The
rule-9 arm ran 27 phases, 7 of which failed; stopping 5 of them saves about
5 x 1,200 x 0.55 s = 3,300 s of 54,621, about 6%; the wider zone admits more
nodes to the phase, and with the gate the extra failures are cheap, so the
model's optimum moves to 0.0035 (tree memo: +4%).

## 3. The search A/B

`../stronger_bounds/search/run_ab5.sh`: the run_ab4 rule-9 arm (322 nodes,
54,621 CPU-s) with `QAP_DNN_CUT_PGATE=4:20` and `QAP_DNN_CUT_ZONE=0.0035`,
everything else unchanged, 4 workers.

| arm | nodes | fathomed | by the DNN bound | by the cut phase | branched | phases run | CPU-s | wall | vs rule 9 |
|---|---|---|---|---|---|---|---|---|---|
| rule 9 (run_ab4) | 322 | 315 | 295 | 20 | 7 | 27 | 54,621 | 4 h 22 min | - |
| rule 9 + gate 4:20 + zone 0.0035 | stopped after 7 nodes (`search/partial_ab5/`) | | | | | | | | |

**Abandoned (Nov 3).** With the deferral dead, the two surviving trims are
worth about 1.1x together (6% from the gate on the rule-9 arm's 7 failed
phases, a modelled 4% from the zone), and a four-hour run to confirm a 1.1x
is not worth the machine time while larger levers (pool inheritance, the
Y12s phase budget, speed) are unmeasured. The gate stays in the code,
default off, replay-validated; switch it on with `QAP_DNN_CUT_PGATE=4:20`
in any later production configuration, where it is free.
