# Projected RLT3 (quadruple consistency): the violation sample, tabled (laptop, Nov 6)

Roadmap item 8: before building any RLT3 cut type, sample the violation of
projected RLT3 inequalities at converged DNN iterates, as step 8 of
`../README.md` did for the RLT2 triples (median about 0.01). Prompted by Jeff
Linderoth's RLT-n plus PDLP suggestion (Nov 5 thread). Budget 90 minutes.

For k facilities the k(k-1)/2 pairwise blocks of Y must be the marginals of
one joint injective placement; k = 3 is step 8's triple consistency (the
projection of RLT2 onto Y), k = 4 the projection of RLT3. `rlt3_sample.py`
separates the most violated inequality <G, Y> <= 0 with potentials G in
[-1, 1] by the same LP as Astra's `../triples/triplet_separate.py`, rounds G
to 1/1024 and takes the right-hand side as the exact maximum of sum G over
all placements, so every reported violation belongs to a valid inequality.
On the triple 4,8,36 of node 13 it reproduces step 8's 0.0886 to five
digits (and 0.0008 for 0,12,34).

## Why the quadruple LP is out of the laptop's reach

- **Merging locations loses the violation.** Mapping all but the s heaviest
  locations of the k facilities to one "elsewhere" cell gives a valid, much
  smaller LP (an inequality on the merged variables is a sum of Y entries),
  but the triple 4,8,36 then shows 0.0001 at s = 10, 14 and 20 against
  0.0886 at the full 47: the violation lives in how the mass is spread over
  many locations, not in the heavy cells.
- **Restricting the potentials loses it too.** Fixing G to 0 where the block
  entry is 0 (2,453 free entries of 6,627) gives violation 0: the inequality
  needs its negative potentials on the zero-mass entries to stay valid.
- So the quadruple LP has to be the full one: 47·46·45·44 = 4,280,760
  placement rows, 6·47² = 13,254 columns, 25.7 M nonzeros. Row generation
  (solve on a subset of placements, scan all 4.3 M for the most violated,
  add them, repeat) converges slowly because the LP is degenerate: with
  3,000 rows added per round, 8 LPs and 424 s without converging; with
  150,000 rows per round the second LP alone exceeded 420 s; the full LP in
  one shot did not finish inside its 480 s limit (HiGHS through scipy, one
  core). A triple (97,290 rows, 6,627 columns) takes 10 s.

A quadruple therefore costs at least 10 minutes on the laptop, and a sample
of the size that makes the question answerable (a few hundred, on three
nodes) is hours of cloud time, not a 90-minute laptop job. Tabled.

## By-product: node 10's triple violations (the first at this node)

Step 8 ran its 60 triples on nodes 13 and d2 only. `n10.Y.bin` is node 10's
(`39:27,30:48,25:49`) uncut 1,500-iteration iterate (`n10.uncut.log`,
certified 48,655.84, the step-2 table's 48,655.9); 30 random triples
(`n10.k3.log`, `n10.k3.json`):

| node | triples | violated (> 0.001) | median | p90 | max | sum |
|---|---|---|---|---|---|---|
| 10 | 30 | 24 | 0.0175 | 0.0412 | 0.0492 | 0.54 |
| 13 (step 8) | 60 | 60 | about 0.01 | | 0.162 | 1.45 |
| d2 (step 8) | 60 | 60 | about 0.01 | | 0.113 | 1.64 |

The same order as nodes 13 and d2 and below the 0.1-0.2 of the regional
triangles' violations: at fixed multipliers the gain available from all 30
is mu times 0.54. The RLT2 projection is satisfied at node 10 as at the
other two.

## If it is picked up again

- Run it where the hours are cheap: the cloud box, `rlt3_sample.py PREFIX
  --fixes SPEC --k 4 --count N --from-triples '../triples/n13.triplet_*.json'
  --seconds 3600`, which spends the first half of the sample extending the
  most violated known triples (RLT3's best case) and prints each
  quadruple's excess over its best sub-triple, the quantity that says what
  RLT3 adds beyond RLT2.
- Or a smarter separator: dual simplex warm-started between rows-generation
  rounds (not reachable through scipy's linprog), or the primal form, a
  distribution over the 4.3 M placements matching the six blocks, by column
  generation with a 4-facility assignment pricing problem.
- The decision rule stands: a median excess near 0.01 closes selective RLT3
  on wil50; only several times the regional triangles' violations would
  justify a cut type.
