# copositive: cuts from the completely positive cone (Nov 25 2026)

Nathan's question: the QAP is an exact linear program over the completely positive (CP) cone (Povh-Rendl 2009), and
the DNN relaxation is its first approximation. Are cuts from copositive matrices a new source of bound?

## Answer

**No, for cuts on small sets of assignments; nothing of the kind is violated at nug12's root or at a wil50
depth-3 node.** Two reasons, one structural and one measured.

- **Structure.** Every entry of Y is a product of 0/1 variables, so on any index set S the feasible Y_S lie in
  conv{yy'} over the feasible 0/1 vectors y, which is inside the CP cone. A copositive cut on Y_S is therefore
  valid for that hull and implied by its facets: it is an exact subgraph constraint (Gaar-Rendl), and at
  |S| = 5 the hull's facets are the triangle and pentagonal inequalities [memory, Deza-Laurent]. Horn-matrix cuts
  (Dong-Anstreicher 2013), Burer-Dong's full 5x5 separation and Hildebrand's other extreme rays are all inside
  this. The same holds for regional aggregates that are 0/1 (the families closed on Nov 4: 5-point, local-7).
- **Measured** (`runs/`): the DNN solution of nug12's root (567.99) and the converged uncut iterate of wil50
  node 13 (`21:28,8:41,28:19`, m = 47, 48,627.29; 1,500 iterations).

| | nug12 root | wil50 node 13 |
|---|---|---|
| Horn cuts, 5-cycles of the strongest scaled entries (`cop5.py`): cycles tested, least value (negative = violated) | 0.9M, +0.027 | 1.4M, +0.055 |
| Horn cuts, local search over all 5-sets (`cop5_ls.py`): starts, least value | 300, +0.011 | 60, +0.046 |
| largest atomic triangle violation, all triples (`hull5.py`) | 0.0208 | 0.0066 |
| largest distance from the exact hull, 5-sets grown from the 20 worst triangles | 0.0052 | 0.0016 |
| the same sets' triples alone | 0.0052 | 0.0016 |

  The hull distance of a 5-set equals that of its worst triangle (a quarter of the triangle's violation, in the
  largest-entry norm): on these sets the exact hull of five assignments holds nothing beyond triangles. Atomic
  triangles at wil50 are violated by at most 0.0066, against 0.073-0.080 for the regional triangles the cut phase
  uses (`../wil50_feasibility/new_cuts/grid/REPORT.md`).

## What this does not cover

- Only the uncut iterate was tested at wil50, not a triangle-phase one; the 5-sets for the hull test were grown
  from violated triangles (a pentagonal violation with no violated triangle inside it would be missed there; the
  Horn searches, which do not start from triangles, found none).
- Cuts with large support. The CP formulation is exact, so every valid inequality is a copositive cut in
  principle, the coarse-problem cut of the grid report included; the cone gives no cheap way to find those. Its
  own hierarchy (Parrilo's K^1) is, in dual form, a consistent family of one conditional DNN matrix per
  assignment: about 2,200 PSD blocks of order 2,200 at wil50. No QAP computation of it was found in the literature.
- Literature read for this (Nov 25, by a subagent, primary sources): Burer-Anstreicher-Dur 2009, Dong-Anstreicher
  2013, Burer-Dong 2013, Bomze-Frommlet-Locatelli 2010, Hildebrand 2012, Povh-Rendl 2009,
  Chen-Lasserre-Magron-Pauwels 2022. None applies CP cuts, K^1 on subsets or a sublevel hierarchy to the QAP.
  Motzkin-Straus cuts on the conflict graph reduce to q(alpha - q) >= 0 for the count q of assignments in the set,
  which the DNN implies.

## Files

`cop5.py`, `cop5_ls.py`, `hull5.py` (arguments: a `-dumpY` file); `runs/*.txt` their output. The iterates:
`latest/dnntest -u 580 -T 1 -H 0 -i 20000 -dumpY runs/nug12_root qaplib/nug12.dat` and
`latest/dnntest -x 21:28,8:41,28:19 -u 48816 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0 -v -dumpY runs/n13_uncut qaplib/wil50.dat`.
