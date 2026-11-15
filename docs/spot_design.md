# Spot: running a run on a fleet that grows and shrinks

Nov 15 2026. The design note for item (c) 1 of `docs/tho40b_results_plan.md`:
workers run on Spot from now on, so the launcher and the coordinator must fit a
fleet whose size changes during a run. Nothing here is built yet; a resumed
tho40-b runs cognac's head, so each change below keeps today's behaviour when
the new parts are absent.

## What works today

- **Spot launch:** `launch.sh` launches Spot by default (`--on-demand` opts out),
  one-time requests that terminate on interruption, and checks the Spot vCPU
  quota (`L-34B43A08`) before it launches.
- **More workers:** `launch.sh RUN --worker-only` adds a worker instance to a run
  that has a coordinating one (role `all`).
- **Interruption notice:** `node.sh` polls the instance metadata for a Spot
  notice every 5 s and drains the worker: each running unit returns its open
  nodes, or its last checkpoint if a node does not end within 90 s, before the
  instance goes. The coordinator test `drain` covers this path by SIGTERM; no real
  Spot notice has exercised it.
- **Restarts:** a coordinator killed at any point replays what it had merged
  (`coord_test.py` coordkill and the crash scenarios), and `launch.sh RUN`
  resumes a run that exists.
- **A dead worker's units** come back when their leases run out (`--lease`,
  30 min), or at once with `qapctl requeue-leased`.

## What assumes a fixed fleet

`--slots` is fixed when the run is created (`qapctl init`, the recipe's slot count)
and read from the run's config everywhere:

| where | what it does with `slots` | on a changing fleet |
|---|---|---|
| `Coordinator.balance` | idle slots = slots - leased units; steals only then | more workers than slots: never idle, so no steals and workers wait on an empty queue; fewer: steals that are not needed |
| `common.plan_split` | keeps reserve x slots units queued | the queue is sized for the wrong fleet: too short for a larger one, too long for a smaller one |
| `publish.py`, `results.txt` | busy fraction = busy wall / (slots x wall) | over or under 100% as the fleet changes |
| `launch.sh` | one coordinating instance (role `all`), coordinator and worker together | a Spot interruption of that instance stops the coordinator too: the run waits for a relaunch |

## The changes, in order

1. **A worker registry.** Each worker writes `runs/RUN/workers/<host>.json` (its
   slots, instance type, market, start time, last beat) with every heartbeat
   (60 s), and deletes it when it drains.
2. **Live slots in the coordinator.** Live slots = the sum over workers that beat
   within the last 3 minutes, or the config's `slots` when no worker writes the
   registry (runs and workers from before the change). `balance`, `plan_split`
   and the busy fraction use it; the busy fraction integrates slot-seconds over
   the run instead of dividing by one count. The public status shows the live
   slots and workers.
3. **A coordinator role.** `launch.sh --coordinator-only`: the coordinator alone,
   on a small on-demand instance (it needs little CPU), so a Spot interruption
   never stops the run. Role `all` stays for single-instance runs.
4. **Faster requeue on interruption.** After the drain, the Spot watcher marks
   the worker `leaving` in the registry; the coordinator requeues any unit still
   leased to a worker that is leaving and has stopped beating, instead of waiting
   out the 30-minute lease.
5. **A fleet that refills itself.** An EC2 Auto Scaling group (or EC2 Fleet) of
   workers, from a launch template built from `node.sh`'s user data, over several
   instance types of one architecture and every Availability Zone in the region,
   capacity-optimized. Replacements join the run on their own. Until then,
   `--worker-only` launches by hand (or from the workflow) do the same job.

One architecture per run: x86 (c7a, c7i, m7a, m7i and their c6 kin) or Graviton
(c7g, m7g), not both. Floating point and the BLAS differ between them, so a
node's bound, and so the tree, can differ slightly with where it ran. Either
way the proof is valid, but one architecture keeps a run's numbers reproducible.

## The test, before a large run

- **Locally (`coord_test.py`):** workers with different slot counts joining and
  leaving mid-run; the totals equal the serial run's (fixed incumbent), no unit
  is lost, the live slots follow the workers, and steals happen when a larger
  fleet runs out of queued work.
- **On EC2:** a subtree run with a known tree on a coordinator-only on-demand
  instance and 2 to 4 Spot workers. The tho40ab subtree (`-F 26:10,8:0`, 22,650
  nodes, 78,804 CPU s with ship.bits = 32) takes about 30 minutes on one c7a.16xlarge.
  Inject interruptions with AWS Fault Injection Service
  (`aws:ec2:send-spot-instance-interruptions`) on one worker, then another, mid-run.
  Pass when the totals match a clean run's, no unit is lost or merged twice, every
  interrupted unit came back within a minute of its notice, and the replacement
  workers joined on their own.

## Decisions for the group

- Is the coordinator on on-demand while every worker is on Spot the right
  split? It costs a small instance for the run's whole length, and saves a stall
  at every interruption of the coordinating instance.
- Which architecture for the next large run? x86 c7a measured 2x faster per node
  than Graviton c7g on nug30 (`infra/phase1/runs/nug30-c7a`, `nug30-a`); Graviton
  Spot is usually cheaper per vCPU. Price per node decides.
