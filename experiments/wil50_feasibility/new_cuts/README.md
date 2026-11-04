# New cut families for the wil50 bound (Nov 4): three research agents

Launched at Nathan's request: three research agents, up to eight hours each,
looking for valid inequalities (or stronger relaxations usable as Lagrangian
cuts) that add enough certified bound at depths 3-5 to make a real dent in
the wil50 tree. About 130 units of extra bound at depths 3-5 is one level, a
factor of about 47 in the tree (`../NOTES.md` section 5).

| directory | mandate | core |
|---|---|---|
| `literature/` | the strongest known QAP relaxations and inequality families beyond DNN, from the literature, and what makes them strong | 1 |
| `grid/` | the location side: wil50's 5 x 10 Manhattan grid, its row/column decomposition, packing and dispersion limits on the grid | 2 |
| `facility/` | the facility side: wil50's flow matrix, facility-subset inequalities, disjunctive and lift-and-project cuts | 3 |

Each agent keeps `NOTES.md` current while it works and ends with
`REPORT.md`. Solver runs are pinned one at a time to the agent's core.

## Coordinator's note: the host changes across container restarts (Nov 4, 09:50)

This box is not one machine. On Nov 3 (AMX work, 15:00-19:30 UTC) it was an
"Intel Xeon Processor @ 2.10GHz" with `amx_bf16 amx_tile amx_int8`; since the
restart at about 05:19 Nov 4 it is an "Intel Xeon Processor @ 2.80GHz" with
no AMX (kernel 6.18.44-fc-v233; v230 before). On the current host the
identical degree-4 fp32 benchmark on node 10 (`dnntest -T 1 -d 4 -g 1.0 -G
1.0 -H 0 -i 100 -u 48816 -x 39:27,30:48,25:49`, the same certified bound to
the last digit) runs at 506 ms per iteration against 409 on the AMX host;
bundle D without AMX at 464. `DNN_AMX=1` falls back to sgemm silently here.
Consequences: CPU seconds measured across restarts are not comparable unless
the host is known, and which host each earlier run used was not recorded
(the stacked node-10 A/B ran through three restarts). From here every long run
should log `grep -m1 'model name' /proc/cpuinfo` and the AMX flags at start
and at every resume.
