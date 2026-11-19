# coord: the coordinator, worker supervisor and run controls

Phase 1 of `docs/distributed/plan.md`: the work model of `bench/parallel.py`
(open nodes in the solver's own state files, resumed with `-r`, returned
with `-o` at a CPU budget or on SIGUSR1) moved onto a queue, a ledger and an
object store, so that nothing waits on a global barrier and any process can
die at any time.

```
python3 -m coord.qapctl init RUN --backend local:DIR --instance qaplib/wil50.dat \
    --params PARAM/nug_deep_xbar9_all.par --ub 48816 --slots 4 \
    --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" [--env DNN_R32=1 ...] [--trace]
python3 -m coord.qapctl coordinator RUN --backend local:DIR      # one per run
python3 -m coord.qapctl worker RUN --backend local:DIR --slots 4 # one per machine
python3 -m coord.qapctl status RUN --backend local:DIR
python3 -m coord.qapctl pause|resume|requeue-leased RUN --backend local:DIR
python3 -m coord.qapctl local RUN --backend local:DIR ...        # init + both, one command
```

Backends: `local:DIR` (one machine: a directory for the store, one SQLite
database in WAL mode for the queue and the ledger) and
`aws:BUCKET[/PREFIX]?region=R&name=N` (S3, two DynamoDB tables, which
`qapctl aws-setup` creates, and SQS queues; `&endpoint=URL` points it at a
local `moto` server for tests). Both implement the same small interface
(`backend.py`). Each run has its own two lanes, large and small (SQS: the queues
`N-units-RUN-large` and `-small`; local: its rows of the messages table), created
with the run and deleted when it finishes or is stopped (a paused run keeps them),
so runs on one backend never see each other's messages. Until Nov 19 every run
shared two queues, and another run's large units kept a run's small ones unread
(`docs/review_nov18_runs.md`, item 1). A run from before that, or one whose queues
are gone, gets them when its coordinator starts, which sends its queued units to
them (`ensure_queues`; a unit sent twice is harmless).

The readout (`coord/readout.py`, Nov 19): one Markdown page per run, `runs/RUN/readout.md`, which
the coordinator writes when the run finishes or is stopped (a dives run's finish writes its
source's: an estimate run's readout is its estimate) and `qapctl readout RUN [--save] [--out F]`
makes at any time. Its first line is the verdict (COMPLETE, PART-WAY, RUNNING, PAUSED, STOPPED);
the warnings follow: a status not written for 15 minutes, dead units, and an estimate's flags
(`coord/dives.py`, `flags`): PART-WAY (pending dives units imputed), GAP, HEAVY TAIL (the 10 largest
dives over half the rest), FEW DIVES (most strata one dive), SETUP-BOUND (setup over the dives'
CPU), DEEP TAIL, WIDE (a jackknife standard error over 25%), each from the tho40-est diagnosis
(`experiments/tho40_est_diag`). `infra/phase1/collect.sh` fetches it; the GitHub workflow
`phase1-readouts` commits every run's to the branch `phase1-readouts` hourly.

An estimate read part-way (Nov 19): a dives run's units go in a random order in one lane, and the
pending units' open nodes are imputed, depth by depth, from the merged units' estimates per open
node at the same depth (a ratio estimator); a delete-a-group jackknife over units gives the
standard error. On tho40-est-dives' record (28 units pending): 983 nodes and 1.9 CPU h imputed (the
diagnosis's 950 and 1.7 by hand), a jackknife standard error of 22% of the nodes and 12% of the CPU.
The estimator's settings (`latest/estimate.C`): `QAP_EST_BUDGET` counts from the end of the
strata's setup; `QAP_EST_MIN_DIVES` dives per branched stratum whatever the budget;
`QAP_EST_SAMPLE_DEPTHS=D:F,...` a sampling rate by depth; `qapctl dives`/`estimate-flow
--dive-env QAP_EST_X=V` passes them to every dives unit (a spec's `ESTIMATE_ENV`).
`QAP_EST_CHILD_DNN=CAP` (Knuth's dives, `kQ`): each child a dive may take weighted by its own
short DNN solve, CAP iterations warm from its parent's state, in place of its QPB bound (the
thesis's design with a predictor fit for a DNN tree; `experiments/estimates_nov19`).

The decision point (Nov 19; `coord/dives.py`, `decision()`): an estimate is worth its CPU only
while it is a small share of the tree it estimates. A dives run marks a pilot, a random tenth of
its units (10 at least, half of a smaller run's: the first of its shuffle, so they are sent
first). When every pilot unit has merged, `estimate-flow` compares the CPU the whole dives run
will take (the pilot's mean per unit times the units) with the rest of the tree as the pilot
alone estimates it. Above `--decide` (default 0.15), the dives stop and the solve goes on to the
end: RUN resumed (and marked `estimate_decision`), RUN-dives stopped with the decision, its
running units cancelled (their workers kill them at the next heartbeat) and its queues deleted;
`worker --follow` goes back to RUN; the readout's verdict is SOLVING, then COMPLETE, with the
decision's line. Otherwise the dives go on to the end, and the estimate ends with its cost (the
estimator's CPU against the rest it found). The pilot is waited for whole because the first
units to merge are the small ones: on tho40-est2's record its first 55 read a rest of 42 CPU h
against the 5,610 its dives found. Replayed on the collected runs (400 random pilots each,
`experiments/estimates_nov19/decide_sim.py`): tai50b-est decides to solve every time, after a
median 3.1 estimator CPU h of the 37.5 its dives took; tho40-est2 never (4.8-8.4%, 10-90%).

The root cache: `init` stores its frontier pass (the state file, the parent states it names, its
traces and report) under `runs/.cache/frontier/KEY/`, KEY a hash of the solver binary, the CPU's
model name, the data file, the files the flags name, the flags, the incumbent, the frontier depth,
the traces and `bnbdrv --print-settings` (the run class aside, as in the settings hash). A later
`init` with the same KEY reads it back instead of running the frontier pass: the root bound is paid
once per problem and settings (72 minutes for tai50b, hours for tai60b). `QAP_ROOT_CACHE=off`
runs it regardless. `tests/coord_test.py rootcache`: nug18 with the DNN and `--frontier-depth 1`,
the second run from the cache, both runs equal to the serial run, DNN iterations by level included.

## How a run works

- **Units.** A unit is a state file holding some open nodes (consecutive in
  stack order, so siblings, which share a parent's warm start and cut pool,
  stay together) plus a ledger row. `init` runs the frontier pass (the root
  or the `-F` node, to `--frontier-depth`) and splits its open nodes into the
  first units.
- **Workers** lease a unit from the queue (the message's visibility is the
  lease; the ledger lease is conditional: queued, or a lease that ran out),
  run `bnbdrv -r in.state -o out.state -t BUDGET -k CKPT -u INC -U inc.txt`,
  and heartbeat while it runs: extend the lease, upload the solver's latest
  checkpoint every `--ckpt-upload` seconds, pass a better published
  incumbent to the solver (inc.txt + SIGUSR2, the solver's new `-U`), and
  send SIGUSR1 when the coordinator asks the unit to return. Results go to the
  store before the ledger; the completion is conditional on still holding
  the lease, so a late or duplicate attempt is discarded.
- **The coordinator** merges each finished unit exactly once: it marks the
  unit merged with the next sequence number, adds its statistics to the run
  totals, and splits the open nodes it returned into child units (the split is
  recorded first and child ids derive from the parent's, so a replay creates
  the same units). The totals are saved to a new file named by the last
  sequence they include, and then `totals.json` points at it, so one write
  switches both. A restarted coordinator replays everything merged after
  that sequence. It keeps
  about `--reserve` x slots units queued, asks the longest-running unit to
  return when slots sit idle with nothing queued (the steal, without a
  barrier), re-queues failed units, verifies an incumbent a unit found (the
  objective recomputed from the original data) before publishing it, and
  writes `status.json`, then `results.txt` and `final.state`.
- **Faults.** A dead worker's units come back when their leases run out (or
  at once with `requeue-leased` after a known machine restart) and resume
  from their last uploaded checkpoint. SIGTERM to a worker drains it: every
  solver is asked to stop at its next node, one still inside a node after
  the grace period is killed and its checkpoint is returned. Solvers run in
  their own sessions, and a solver that dies of a signal while draining also
  returns its checkpoint.

- **Traces.** With `init --trace` every attempt uploads `trace.txt`
  (`QAP_DNN_TRACE`, a line per DNN node) and `nodes.txt` (`QAP_NODE_TRACE`, a
  line per node: format in `latest/bnb.C`) next to its `out.state`, and the
  frontier pass its own under `runs/RUN/frontier/`. The whole tree is
  `frontier/nodes.txt` plus, for each merged unit, the `nodes.txt` of the
  attempt the ledger records: an attempt that resumes a checkpoint starts from
  the lines that checkpoint accounts for (`ckpt.nodes.txt`) and a killed
  solver's trace is cut to its checkpoint, so each node is in exactly one of
  them. Other attempts' files overlap them. `qapctl nodes RUN --out FILE[.gz]` writes that whole tree, and
  `infra/phase1/collect.sh` keeps it with a finished run's records as `nodes.txt.gz`
  when it is under 25 MB. The DNN trace is carried across checkpoints the same way
  (`ckpt.trace.txt`), uncut, and `qapctl export RUN --out DIR` (`coord/export.py`)
  writes both whole: `nodes.txt.gz` and `dnn.txt.gz`, a DNN line for each node whose DNN
  bound ran, matched by its fixes, in the same order. The collect exports every
  `--trace` run, finished or not (part-way: the units merged so far), to
  `infra/phase1/exports/RUN/`, out of git; the GitHub workflow keeps it as the job's
  artifact for 90 days. A public `--trace` run also publishes that tree as it
  grows (`coord/tree.py`, with each status publish): every merged node, in
  shards under `public/runs/RUN/tree/` that brixius-web's tree viewer reads;
  only the shards a new node changed are written again, and a restarted
  coordinator, or `qapctl publish`, rebuilds the same tree with the same ids
  (`tests/node_trace_test.py coord coordkill`). `python3 -m coord.tree TRACE DIR
  --run NAME` writes a trace's tree to a directory the same way.
- **Parent states and disk.** A node's DNN bound starts from its parent's final
  state when it ships (`ship.bits`), and a unit's state files would carry a copy
  per open node: `coord/warm.py` keeps one per parent in `runs/RUN/warm/`
  instead (`docs/parent_state_store.md`). What goes when: a unit attempt's
  directory on the worker when the attempt is over; a merged unit's state files
  once the totals that hold it are saved; a parent when no unmerged unit can
  start from it; the run's parents when it finishes. `qapctl stop` gives a run
  up, `qapctl purge-states` deletes a stopped or finished run's state files.

- **Estimate runs** (`coord/dives.py`; `docs/estimates_plan.md`). An estimate is the start
  of a solve. The solve runs for a wall-clock budget, then `qapctl pause` (or `estimate-flow`)
  pauses it, and its units return their open nodes. Then `qapctl dives SRC --as NEW` makes a
  dives run: one unit per queued unit of SRC, which reads SRC's state file and parent states
  where SRC keeps them, never copies or deletes them. Its worker runs the estimator's survivor
  dives (`bnbdrv -e -r`, `QAP_EST_SURV`, `latest/estimate.C`) for the unit's budget of wall
  seconds and keeps their JSON lines (`dives.jsonl`). Its coordinator merges units without
  splitting them. Its `results.txt` is `qapctl estimate NEW`'s report:
  - SRC's merged totals plus every open node's exact part and its dives' mean, by depth;
  - a bootstrap interval;
  - the deepest level the dives reached, and the open nodes no dive reached.

  `qapctl estimate-flow RUN --solve-wall S --dive-budget T [--decide F]` does all of it in place
  of the coordinator, restartably; `worker --follow` goes from RUN to RUN-dives (and back, when the
  decision point stops the dives for the solve: above). SRC stays paused:
  `qapctl resume` sends its queued units again and the solve goes on. The dives run has
  queues of its own, so SRC's queued messages wait in SRC's. `tests/dives_test.py [--flow]` (nug17): the estimate within 15% (1.03
  and 1.06 on the first runs), and the resumed solve's totals equal a serial run's; `--decide`
  (in `--flow` too): the decision point at 1% stops the dives after the pilot, the solve finishes
  with the serial run's nodes on the worker that came back, and the readout says so. On EC2:
  `QAP_ESTIMATE=SOLVE_S:DIVE_S[:B]` on every instance (`infra/phase1/node.sh`; start.sh asks
  for it under the other settings). On the laptop: `ESTIMATE=... infra/laptop/solve.sh PROBLEM`.
  `node.sh` in its local mode against moto, from a clean archive, gave three outcomes (Nov 18):
  - **nug18, `QAP_ESTIMATE=8:4`:** paused at 8 s, 15 dives units holding 409 open nodes; the
    worker followed. The estimate: 131,122 nodes (5-95%: 125,933 to 137,343) against 136,207,
    and 40.2 CPU s against 40.1.
  - **nug18 at 25 s:** the tree closed inside the budget. That is the solve's answer, and the
    check passed.
  - **nug17 at 3 s:** the tree closed while the units returned. Nothing was open, so the run
    finished as a solve.

  Workers see a pause at their next status look, up to 10 s later, and finish the unit they hold.

With a fixed incumbent the totals equal a serial run's exactly, whatever
fails (`tests/coord_test.py`, `make coord`: nug17 with 0.3-second units,
about 90 units a run):

| scenario | what happens | totals |
|---|---|---|
| normal | units return at the budget, steals | = serial |
| kill9 | the worker and its solvers SIGKILLed mid-run, a new worker | = serial |
| kill9ckpt | the same with longer units: the killed units resume from their uploaded checkpoints | = serial |
| coordkill | the coordinator SIGKILLed mid-run and restarted | = serial |
| crash_marked, crash_applied, crash_saving | the coordinator stops (`QAP_COORD_CRASH`) in its 20th merge: after marking a unit merged, after applying it, between writing the totals and `totals.json`; a new one replays | = serial |
| drain | SIGTERM to the worker's process group mid-run | = serial |
| noub | no incumbent: the run finds the optimum and publishes it verified | optimum |
| aws | normal, on S3/SQS/DynamoDB (a local moto server, region eu-north-1, where a new bucket must name its region) | = serial |
| tworuns | two runs on one backend at once, each on its own queues: run a paused partway with units queued (its worker drained), run b finishes meanwhile, then a resumed with a new worker | both = serial |
| awstworuns | tworuns on S3/SQS/DynamoDB (moto); no queue left once both finish | both = serial |
| stopped | `qapctl stop` partway; a run stopped before any worker came has its queues deleted; another run then finishes | = serial |

On wil50 (`experiments/wil50_feasibility/coord_search/README.md`): the
node-10 subtree with the production settings, one 4-slot worker: the same
322-node tree as `bench/parallel.py`, 1.09x fewer DNN iterations (23 units
and no barrier keep more parents' warm starts: 228 warm-started nodes against
176), the same 81% busy (one subtree's tail, which no driver can split). Its
CPU and wall ratios against the parallel.py run (1.29x, 1.28x) are mostly a
faster host after a container restart (1.23x per iteration on identical
nodes).

The public status (`publish.py`): every `--publish-every` seconds (300) and
at the end, the coordinator writes `public/runs/RUN/status.json`, `totals.state`
(the merged statistics, every per-level table of the solver's report),
`report.txt` (that report) and `history.json` (one sample per publish), for
brixius-web's live run page; `--visibility private` publishes nothing, and
`qapctl publish RUN` publishes at once. `infra/phase1/publish_setup.sh` serves
`public/` through CloudFront with the bucket kept private.

On EC2, `infra/phase1/` (README.md there) runs a run on one instance: the
launcher, the instance script with a Spot-notice watcher that drains the
worker, and the acceptance check against `bench/parallel.py`.

Not built yet (later phases of the plan): the Terraform and the
auto-scaling group, the dashboard, pruning queued units against a new incumbent (the solver fathoms
them at pop time anyway), per-depth ETA, ledger snapshots, the AWS
`count_by_status` without a scan.
