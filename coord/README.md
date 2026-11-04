# coord: the coordinator, worker supervisor and run controls

Phase 1 of `docs/distributed/plan.md`: the work model of `bench/parallel.py`
(open nodes in the solver's own state files, resumed with `-r`, returned
with `-o` at a CPU budget or on SIGUSR1) moved onto a queue, a ledger and an
object store, so that nothing waits on a global barrier and any process can
die at any time.

```
python3 -m coord.qapctl init RUN --backend local:DIR --instance qaplib/wil50.dat \
    --params PARAM/nug_deep_xbar9_all.par --ub 48816 --slots 4 \
    --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --env QAP_DNN_POOL=1 ... [--trace]
python3 -m coord.qapctl coordinator RUN --backend local:DIR      # one per run
python3 -m coord.qapctl worker RUN --backend local:DIR --slots 4 # one per machine
python3 -m coord.qapctl status RUN --backend local:DIR
python3 -m coord.qapctl pause|resume|requeue-leased RUN --backend local:DIR
python3 -m coord.qapctl local RUN --backend local:DIR ...        # init + both, one command
```

Backends: `local:DIR` (one machine: a directory for the store, one SQLite
database in WAL mode for the queue and the ledger) and
`aws:BUCKET[/PREFIX]?region=R&name=N` (S3, two SQS queues, two DynamoDB
tables; `qapctl aws-setup` creates them; `&endpoint=URL` points it at a
local `moto` server for tests). Both implement the same small interface
(`backend.py`).

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
  the same units). The totals are saved with the last sequence they include;
  a restarted coordinator replays everything merged after that. It keeps
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

With a fixed incumbent the totals equal a serial run's exactly, whatever
fails (`tests/coord_test.py`, `make coord`: nug17 with 0.3-second units,
about 90 units a run):

| scenario | what happens | totals |
|---|---|---|
| normal | units return at the budget, steals | = serial |
| kill9 | the worker and its solvers SIGKILLed mid-run, a new worker | = serial |
| kill9ckpt | the same with longer units: the killed units resume from their uploaded checkpoints | = serial |
| coordkill | the coordinator SIGKILLed mid-run and restarted | = serial |
| drain | SIGTERM to the worker's process group mid-run | = serial |
| noub | no incumbent: the run finds the optimum and publishes it verified | optimum |
| aws | normal, on S3/SQS/DynamoDB (a local moto server) | = serial |

On wil50 (`experiments/wil50_feasibility/coord_search/README.md`): the
node-10 subtree with the production settings, one 4-slot worker: the same
322-node tree as `bench/parallel.py`, 1.09x fewer DNN iterations (23 units
and no barrier keep more parents' warm starts: 228 warm-started nodes against
176), the same 81% busy (one subtree's tail, which no driver can split). Its
CPU and wall ratios against the parallel.py run (1.29x, 1.28x) are mostly a
faster host after a container restart (1.23x per iteration on identical
nodes).

Not built yet (later phases of the plan): the Terraform and the Spot
instance-metadata handler (the drain path is what it will call), the
dashboard, pruning queued units against a new incumbent (the solver fathoms
them at pop time anyway), per-depth ETA, ledger snapshots, the AWS
`count_by_status` without a scan.
