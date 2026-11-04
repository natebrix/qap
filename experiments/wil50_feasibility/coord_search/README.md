# The node-10 subtree through the coordinator

`run_coord_search.sh`: the pool-inheritance configuration of
`../pool_search/` (stacked arm + `QAP_DNN_POOL=1`), run through the new
coordinator (`coord/`, plan.md phase 1) with one 4-slot worker instead of
`bench/parallel.py`. Units: CPU budget 3,600 s, new units of about 1,800
CPU-s, 2 x slots queued, steal after 10 s, lease 600 s, heartbeat 60 s,
checkpoint uploads every 300 s. Nov 4, 18:21-19:57 UTC, no restart.

## Result

| driver | nodes / fathomed / DNN-fathomed | CPU s | wall s | slots busy | tasks or units | warm-started nodes | first-solve iterations | cut-phase iterations |
|---|---|---|---|---|---|---|---|---|
| `bench/parallel.py` (pool_search, 15:36-17:39) | 322 / 315 / 307 | 23,905 | 7,342 | 81.4% | 80 tasks, 3 checkpoints | 176 | 46,550 | 17,600 |
| coordinator | **322 / 315 / 307** | 18,496 | 5,736 | 81.6% | 23 units, no barrier | **228** | **39,770** | 19,100 |

- **The same tree, exactly** (node, fathom and DNN-fathom counts; FW
  iterations 30,816 in both).
- **The CPU and wall ratios (1.29x, 1.28x) are mostly the machine, not the
  driver.** The container restarted at 17:49 onto a host with the same CPU
  model name that is faster: on the 90 nodes that ran with the same number of
  iterations, the same warm start and the same pool in both runs, an
  iteration took 1.23x longer in the parallel.py run (median; quartiles 1.16
  and 1.29), and node 10 itself (cold, alone, 300 iterations) took 147 s
  there against 114 s here. CPU comparisons across container restarts need
  a benchmark on each host, not just the model name (`hosts.txt` logs the
  name only).
- **What the coordinator itself changed: 1.09x fewer DNN iterations**
  (64,150 to 58,870), host-independent. With 23 units and no checkpoint
  barrier instead of 80 tasks and three barriers, 73 more nodes started from
  their parent's in-memory Y and Z (21 fewer the other way), so the first
  solves needed 1.17x fewer iterations (the pooled depth-5 children: 68
  iterations on average against 128); the cut phases ran 1,500 more.
- **The slots were not busier.** Idle slot time 18.4%: 89% of it in the
  last 25 minutes, when the queue was empty and one to three units were
  finishing (the last one ran alone for 13 minutes: a node in its cut phase,
  which no driver can split), and 10% in the first five minutes, while node 10
  was bounded alone. In between, all four slots were busy. On one subtree
  with four slots, the critical path decides the busy fraction; the
  barrier and steal costs the coordinator removes show at fleet scale, where
  thousands of subtrees overlap their tails.

## Reading

The coordinator is a correct drop-in for `parallel.py` on a real wil50 run:
the same tree, a 1.09x saving in iterations from keeping warm starts, and
everything the fleet needs (leases, checkpoint resume, drains, incumbent
publication, a restartable merge). Its 1.15-1.35x fleet-efficiency estimate
(`roadmap/node_speed.md` section 4) stays an estimate: a 4-slot subtree run
cannot show it. Phase 2 of plan.md (100 Spot cores with induced
interruptions) is the place to measure it.
