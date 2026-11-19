# Design review: running experiments on phase 1, November 18 2026

Nov 18 put four runs on the phase 1 fleet: the `tho40-est` estimate, its dives run, and two sko42 arms
that differ only in `cut.zone`. All four have stopped. Each hit a different weakness in how runs are
set up, shared and collected. None of them makes a proof wrong. Together they mean an A/B on the fleet is
not yet a controlled experiment.

Everything is parked for now. These are proposals, not changes. Items are ranked by how much each one
undermines a comparison.

## State at writing

| Run | State | Notes |
|---|---|---|
| `tho40-est` | paused (by design, for its estimate); 559 units merged, 550 queued | instances terminated; resumable from its store |
| `tho40-est-dives` | 524 of 552 units merged; 28 never leased | estimate read part-way (below) |
| `sko42-z10-261118` | Spot reclaimed it at 22:03 UTC, 30 min in; 187 units merged | `cut.zone` 0.010; launched without `--trace` (`runs/sko42-z10-261118/NOTE.md`) |
| `sko42-z25-261118` | Spot reclaimed it at 22:01 UTC; 147 units merged | `cut.zone` 0.025, traced |

The `tho40-est-dives` estimate, read part-way with `qapctl estimate`:

- **Whole tree:** 1.87M nodes and 4,500 CPU h, with a 5–95% range of 3,500 to 5,550.
- **It undercounts.** The deepest dive reached depth 7, yet the estimate puts about 1M nodes at depth 6.
  Whatever lies below depth 7 is not represented, so read 4,500 CPU h as a floor.

The laptop sko42 run at `cut.zone` 0.025 (`infra/laptop/runs/sko42-zone25/`) is unaffected and still
running.

## 1. Runs that share the queues starve each other's small lane

**What we saw.** The last 28 `tho40-est-dives` units were never leased. All 28 were in the small lane,
at depth 5 or 6, and their messages had been sent. At the same time, the 128 tho40 slots sat idle for
more than 20 minutes. The large queue held 557 visible messages: the two sko42 arms' units plus
`tho40-est`'s 519 paused ones.

**Why.** Every run on a backend shares `qap-units-large` and `qap-units-small`.

- `SqsBackend.receive` (`coord/backend.py:467`) asks the large queue first and only tries the small
  one when large comes back empty.
- A worker that receives another *running* run's message hides it for `FOREIGN_HIDE_S = 1` second
  (`coord/worker.py:46`, `:506`), so the message is visible again almost at once.
- A paused run's message is hidden for 300 s, which is long enough.

So while any other live run has units in the large lane, a worker almost never reaches the small one.

**What it does to an A/B.** Two arms on the same backend are not independent: each one's large-lane
units delay the other's small-lane units. Total CPU and the tree should be unaffected. Wall time, cost
and the order in which units run are affected, and the order can change when the incumbent is found.

I first read the sko42 arms' "34 of 64 slots leased" as starvation. That was wrong: it was the arms'
last status before Spot reclaimed them, written while they were still ramping up. Starvation is proven
only for the dives run's tail.

**Options:**

- **(a) Per-run queues.** Name them `{name}-{run}-units-{lane}`, created at `qapctl init` and deleted
  when the run finishes. This removes the problem and the foreign-message path entirely. It costs some
  queue lifecycle work, and a dives run would read its source's queue.
- **(b) Don't let the large lane block the small one.** Either keep looking within the same receive
  when a message is foreign, or alternate which lane a slot asks first.
- **(c) Hide another live run's message for longer,** for example 30 s, as the paused case already
  does.

(b) and (c) are small changes but leave the coupling in place. For experiments, (a) is the one that
makes arms independent.

**Status, Nov 19: done, (a).** Each run has the queues `qap-units-RUN-large` and `-small`, created
with it and deleted when it finishes or is stopped; the foreign-message code is gone. A run from
before gets its queues when its coordinator next starts (its queued units are sent to them).
`tests/coord_test.py` `tworuns` and `awstworuns`: one run paused with units queued while the other
finishes, then resumed (`coord/README.md`).

## 2. An A/B is launched by typing at prompts

**What we saw.** The baseline arm's recipe has
`COORD=--budget 600 ... --checkpoint-every 120`, without the `--ckpt-upload 300 --trace` that was
meant to follow. The treatment arm's recipe has them. So the two arms differ in tracing, and tracing
costs a few % of CPU, as well as in `cut.zone`. The cause of the truncation is unknown: either the
input was cut off, or it was typed short between meetings. `start.sh` reads it with `read -r -p`
(`infra/phase1/start.sh:56`).

**Proposal.** Make an A/B a single command from files:

```
start.sh sko42 --ab experiments/estimates/sko42/sko42.par experiments/estimates/sko42/sko42_zone25.par
```

or a recipe file per arm that is committed before launch. The launcher would refuse arms whose recipes
differ in anything other than `PARAMS`, and its summary would show the differences between the two
recipes. The prompts would stay for one-off runs.

**Status, Nov 19: done.** `start.sh --spec F`: every answer from a committed spec; `ARMS` in it makes
the A/B, whose arms differ in their params files only by construction, the difference shown before
the launch (`infra/phase1/specs/`). Item 5 with it: the run's record is committed and pushed at launch.

## 3. Spot for comparisons

**What we saw.** Both sko42 arms were in eu-north-1c and were reclaimed two minutes apart, 30 minutes in
(`instance-terminated-no-capacity`). The runs can resume. But an arm that keeps getting interrupted runs
for hours longer, and resumes add requeues and lost unit attempts that an uninterrupted arm doesn't have.

**Proposal.**

- An A/B defaults to on-demand. That's about $2.50 to $3 an hour per c7a.16xlarge, against $0.63 on
  Spot.
- Failing that, it defaults to Spot in different zones, with the interruptions and the CPU of lost
  attempts reported next to the totals.
- `collect.sh` reports the instances' `StateReason` and the Spot request status, so an interruption is
  never mistaken for a hang.

Here, the stopped status files looked like a stalled run until the instance state was checked.

**Status, Nov 19: a coordinator on demand and a self-replacing Spot worker group** (agreed with the
proposal's author in place of all on demand): `infra/phase1/README.md`, Starting a run. A reclaimed
worker is replaced by the group, the root bound is on the on-demand coordinator and in the root
cache, and `collect.sh` shows each instance's `StateReason` and the group's activities. For an A/B,
`--workers-on-demand` puts the group on demand.

## 4. Collect can come back without the traces, and say so only in passing

**What we saw.**

- **tai60b, first collect:** no traces, and nothing in the run folder showed they were missing.
- **tai60b, second collect:** `qapctl export` succeeded. Then `qapctl nodes` hit an S3 read timeout and
  left `nodes.txt.gz.tmp` behind. That step re-reads the 41,784 units the export had just read, and took
  74 minutes before failing.
- **tho40-est collect:** the SSO token expired partway, so `traces:` printed a `TokenRetrievalError` as
  one line among many.

**Proposals:**

- Take the run's `nodes.txt.gz` from the export, which is the same lines from the same source, instead
  of reading the store a second time. That is what was done by hand for `tai60b-261117` (2984bb36).
- Check the session at the start of a collect and stop if fewer than about 30 minutes are left on it.
- End with a one-line verdict, such as `complete`, `INCOMPLETE: no traces`, or `INCOMPLETE: no node
  trace`, and write the same verdict to `runs/RUN/collect.txt`, so a missing piece is visible in git.

**Status, Nov 19: done,** with the run's readout (`coord/readout.py`) printed before the verdict, and
the readouts of every run committed hourly to the branch `phase1-readouts` without approval
(`.github/workflows/phase1-readouts.yml`).

## 5. Records with no owner

**What we saw.**

- `runs/tho40-est/` on the laptop has `hosts/`, `region.txt` and `status.json`, but no `recipe.sh` or
  `launches.tsv`. The coordinator's log, however, says it was launched with
  `runs/tho40-est/recipe.sh on the laptop`.
- None of this run is committed.

**Proposal.**

- `launch.sh` commits, or at least stages, `runs/RUN/recipe.sh` and `launches.tsv` the moment it
  launches.
- A run started from another session commits its record before the instances come up, so either
  session can collect it.
- Whoever launched `tho40-est` should commit its recipe and the estimate's output.

## Before the next A/B

1. Per-run queues, or (b) at the least.
2. An A/B launched from two committed params files, with the launcher refusing any other difference.
3. On-demand instances, or interruptions reported with the results.
4. A collect that states whether it is complete.

Then rerun the sko42 `cut.zone` pair from scratch rather than resuming the current arms. They differ in
tracing, and both have been interrupted once.
