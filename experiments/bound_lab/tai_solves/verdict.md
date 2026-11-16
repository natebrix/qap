| instance | verdict | optimum found and proved | nodes | DNN iterations | CPU-s | wall | first published solve |
|---|---|---|---|---|---|---|---|
| **tai35b** | **solved, valid** | 283,315,445 (= `QAPSOLN/tai35b.sln`), no open node | **922** | 125,232 | 5,275 | 10 min root + 35 min on 3 cores | 2.6 million nodes |
| **tai40b** | **solved, valid** | 637,250,948 (= `QAPSOLN/tai40b.sln`), no open node | **733** | 75,582 | 4,862 | 21 min root + 49 min on 1-2 cores | 278,000 nodes |

- Both from an incumbent one above the optimum, with tai50b's settings plus `dnn.beta_sched = w=0.36,wc=2.88`.
- **No sign of an invalid cut.** On the optima's paths 32 and 37 nodes were bounded and none fathomed; the
  largest certified bounds there are 283,315,444.9994 and 637,250,947.9995. No search ended above the optimum.
@FILE:verdict_second.md@
- **A solver fault was found and patched on the way** (section 5): with `wc=2.88`, a child's first solve from
  its parent's state and pool sometimes fails outright (a certified value of -3e11, or none). The bound stays
  valid, but the node gets no phase and its children inherit the failure; tai35b's tree blew up at depth 8
  until the failed solve was redone without the pool, and cold if that fails too (`pool_retry.patch`). It did
  not occur on tai40b or in the tai50b session.
  A tai50b run with `wc` needs the patch.
- What it means for tai50b (section 6): the code and the mechanism hold on two instances first solved at
  great cost, and both trees are one path with a small fringe, as the tai50b report assumed. They are easier
  instances than tai50b (roots 0.89% and 0.18% under the optimum against 1.15%), so they do not test its number.
