# The next tho40 and tai50b runs, after the bound lab (Nov 16)

`experiments/bound_lab/README.md` (Nov 15-16) puts tho40's tree at 226 c7a CPU-days with the grid
cut phase (no solver change), against 1,308 as tho40-b runs, and tai50b at about 100 core-hours
with cluster-region cuts. This page decides what that means for tho40-b and plans a tai50b run of
record. The subtree A/B in `experiments/tho40b_analysis/giveup_ab/` settles tho40's settings first.

## tho40-b: let it end, then a new run

**What is left under the current settings.** At its 72-hour limit (Nov 17 16:46 UTC) tho40-b will have
done about 190 c7a CPU-days: 103.7 at 39.4 h, at 63 CPU-days a day. The Knuth estimates of the
whole tree under its settings are 1,308 (laptop, 3,322 dives, read as low) and 1,036 (se 463; the bound
lab's 48 dives). That leaves about 850-1,120 CPU-days: 14-18 more days on its one c7a.16xlarge.

**Under the cut phase** the whole tree is 226 CPU-days (se 39, 51 dives; bootstrap 179-273, too narrow
by the study's own account), and about 110-215 with the hand-over rule and the phase penalty. The part
tho40-b finished would cost those settings about a fifth of what it cost tho40-b: perhaps 35 CPU-days
of the 226.

**A switch of tho40-b's remaining frontier would take:**
- **New settings in a run in flight.** The coordinator fixes the params file at `qapctl init`. Changing
  it means rewriting the run's record in the store, which no command does today.
- **`--allow-settings-change` on every resume.** The policy hash covers every setting's value, so a
  unit's state file written under the old settings is refused without it. The no-patch block changes
  values only, so it needs no new binary but still changes the hash.
- **A new binary on every worker at once, for the full block.** The integration build has 80 settings
  where cognac has 71 (`integration/README.md`, finding 3).
- **Pools in the state files.** `pool.on` makes the state files carry cut pools. The units written
  before the switch carry none, so their nodes start pools from scratch. Pools pass through the
  coordinator intact: `giveup_ab/README.md`, the same 38-node tree as one process over 38 units.
- **A proof that rests on two settings regimes**, recorded per unit.

It saves about 35 CPU-days at most, about half a day on one instance, or $5-10 at Spot rates
($0.0077 a core-hour on a c7a.4xlarge, Nov 15). **Not worth the change in flight.**

**Recommendation.**
1. tho40-b ends at its limit, as decided. The Nov 17 collect exports its traces.
2. Settle the settings on the subtree: `giveup_ab`'s cut06, np10 and full12 against g200, g100 and
   g100z0 (the laptop).
3. A fresh run, tho40-c, with the chosen block, from the root. It needs the path check first
   (`experiments/tho40b_analysis/path_check/dive.py`, its params given with `--params`). With the
   full block, it also needs the integration build reviewed and merged (`integration/REVIEW_cognac.md`).
4. **Size:** 226 CPU-days is 3.6 days on one c7a.16xlarge, or 1.8 on two (`infra/phase1/start.sh
   tho40`, two Spot instances: the quota is 128 vCPUs). The combined arm, if the
   subtree confirms it, would take 1.7-3.4 days on one. That is the "couple of days" target, on the
   quota there is.
5. **Read the dives as likely low.** The estimator that produced 1,308 read about half the size of a
   subtree whose size was known (`knuth_tho40/README.md`). The bound lab's dives used its weighting
   too, so tho40-c's estimate should be redone with `QAP_EST_CHILD_QPB=1` before the quota is sized.

## tai50b: what a run of record needs

The bound lab's evidence: 275 Knuth dives and a model put the tree at 5,000-13,000 nodes and 60-250
core-hours on the cloud box, about 1.4 times that on c7a. tai35b and tai40b were solved end to end in
922 and 733 nodes by two region families each. A level-by-level rehearsal runs on
`claude/cool-babbage-ilnogj-tai50b` (`experiments/bound_lab/tai50b_bfs`). It is not a proof of record:
three patch files, outside the coordinator.

**The build.** One commit, built on every worker from the run's source archive:
- cognac with b86352e's `latest/`: `cut.regions`, and the `-e -r` fix for estimates.
- `pool_retry.patch` with the review's changes:
  - a settings key, so the record shows it;
  - each retry recorded in the DNN trace;
  - the trigger tied to the cut-phase gate.
- The review's two small fixes before a run that depends on them: the `-e -r` root bound set before
  step 1, and `cut.regions` validated at startup.

**The regions file.**
- **On every worker:** commit it in the repository, next to the run's params file. `launch.sh` ships
  `git archive HEAD`, so every worker gets the same file at the same path. The path is relative to
  `latest/`.
- **In the settings record:** today only the path is hashed (finding 4). Read the file at startup, add
  its hash to the record (a `cut.regions.hash` entry beside `par.hash`), and fail at startup, not at the
  first phase, if it is missing or has out-of-range indices (finding 14). A resume with a different file
  then fails the `-r` check.
- **The run record:** the commit and the file's sha256.

**The root.** The bound lab's root is a long one (`cut.stall = 0`, `cut.cp_iters = 6000`; 6,470
iterations, 72 minutes on one core), separate from the settings below it. The coordinator takes one
params file per run. Either:
- seed the run from that root's state file, which `qapctl init` cannot do today; or
- give the root its own settings by depth in the params file.

The first is simpler and keeps the root measured in `tai50b_size` exact.

**Traces and storage.**
- `--trace`, both traces, exported at the end (the collect does it since 3c2516c).
- A tai50b parent's DNN state is about 123 MB in the rehearsal (n = 50), and the store keeps one copy per
  parent until its children are done. Size the disk for the open frontier, a few dozen nodes wide at
  depths 4-5 (`tai50b_size`, section 7): 100 GB per instance is ample.

**The incumbent and the checks.**
- **Incumbent:** the best known + 2 (458,821,519), the rule for new runs (`infra/phase1/README.md`,
  Starting a run): the run must find the best known itself and prove nothing better exists. tai50b's
  values are integers but not proven even, so + 1 (the rehearsal's) is the tightest; + 2 adds the
  nodes whose bound lies within one more unit.
- **Path check first:** `path_check/dive.py` toward `QAPSOLN/tai50b.sln` with the run's settings.
  `dive.py` handles grid symmetries only, and tai50b has none: it needs a no-symmetry mode, a few lines.
- **A rehearsal under the coordinator:** tai35b and tai40b end to end (922 and 733 nodes, about 1.4
  CPU-hours each): the node counts must match the bound lab's.
- **A second route,** another region family, as for tai35b and tai40b: a cheap independent
  confirmation once the first run ends.

**The size of the run.** About 140 c7a core-hours (60-350): a few hours on one c7a.16xlarge, a few
dollars on Spot. The wall is set by the depth structure, not the cores: the root is one 72-minute solve,
and depths 4-5 hold only a few dozen open nodes, so 64 slots will not all be busy. A c7a.8xlarge would
do.
