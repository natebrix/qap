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
