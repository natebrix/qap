# Why the signal chose c26 at tai50b's root (Nov 20 2026)

The laptop's whole-tree run (`../runs/tai50b-ls/`) failed because the signal ranked c26 first at the root
(12 predicted survivors; 28 of its 50 children branched) and rule 9's r31 23rd (15 predicted; 10 branched).

**Answer [m]: the local search, not the threshold. Sixteen of c26's children have good completions that a search
started from the incumbent does not find.** With one incumbent-free search per cell added, and every local
optimum credited to all the cells it satisfies, the signal ranks r31 first of 100 lines and counts 26 for c26.

## The evidence

`tai50b_root_all.txt` (the C signal, `t0all.C` = lsla's `t0.C` printing every child; slack 3,372,321, 30 kicks
+ 3,000 on cells with 1 < x <= 2.5), joined with the depth-1 outcomes of the laptop run (c26) and tai50b-r9
(r31) in `root_children.tsv`:

| line | branched | of them x <= 1 | of them x > 1 | fathomed with x <= 1 |
|---|---|---|---|---|
| c26 | 28 | 12 (facilities 19-30, x 0.00-0.25) | **16 (facilities 34-49, x 7.6-9.0)** | 0 |
| r31 | 10 | 8 | 2 (x 1.07, 1.77) | 7 (x 0.31-0.81, all by the cut phase or first solve within 5k of the cutoff) |

The 16 are 418k-766k short of the cutoff. Read literally, x = 8 says their best completion is 27M above the best
known, which would make the bound's gap there 8 times the root's.

`hardsearch.py` (3 random starts, the fix forced, 3,000 kicks each, 4 s a child) → `hardsearch.txt`:

| child | outcome | signal's x | random-start ILS x |
|---|---|---|---|
| 39:26, 34:26, 37:26, 44:26 | branched | 7.6-9.0 | **0.61-0.66** |
| 31:40, 31:2 | branched | 1.77, 1.07 | 0.56, 0.70 |
| 0:26, 15:26, 31:1 | fathomed | 1.79, 6.99, 3.84 | 4.69, 2.12, 5.46 |

Taking the smaller of the two x, all nine fall on the right side of x = 1. The incumbent with one exchange,
a descent and 30 kicks of 6 exchanges stays in the incumbent's arrangement; these children need another one.

## A fix, tried in Python (`harvest.py`, `harvest_node.py`)

Per free cell: (a) the signal's base search (incumbent, the cell forced, descent, 30 kicks); (b) a random
completion with the cell forced, descent, 100 kicks. Every local optimum either search visits is credited to all
the cells it satisfies. Predicted survivors of a line = its cells whose best cost is within the slack.

| node | line: measured branching children | base search alone: count (rank) | with (b) and the harvest: count (rank) |
|---|---|---|---|
| tai50b root (`harvest_k100.txt`, 147 CPU s) | c26: 28 | 12 | **26** (38 of 100; tp 26, fp 0, fn 2) |
| | r31: 10 | 11 | 14 (**1** of 100; tp 9, fp 5, fn 1) |
| tai60b depth 7, 10:39 (`harvest_d7_10_39.txt`, 199 s) | c59: 3 | 3 (15) | 3 (9 of 106) |
| | r11: 4 | 4 (17) | 4 (11) |
| | c47: 5 | 2 (6) | 3 (7) |
| | c54: 9 | 9 (87) | 9 (84) |
| tai60b depth 8, 14:54 (`harvest_d8_14_54.txt`, 186 s) | c59: 3 | 3 (16) | 3 (8 of 104) |
| | r3: 4 | 4 (17) | 4 (9) |
| | c2: 7 | 7 (53) | 7 (50) |
| tai60b depth 9, 24:2 (`harvest_d9_24_2.txt`, 170 s) | c30: 12 (the signal's miss) | 2 (9) | 4 (25 of 102) |
| | r8: 4 | 3 (12) | 3 (3) |
| | r3: 4 | 4 (26) | 4 (15) |

Ranks are by count alone, ties in the order found (no chain-model tie-break here), so a rank inside a tie means
little; the counts are the result. CPU is Python around `landscape/ls.c`, one core of the cloud box; the refined
C signal the run used took 267 s at this root on the laptop.

- At tai50b's root the fix repairs the miss: the count for c26 goes from 12 to 26 against a measured 28.
- At tai60b's depth-7 and depth-8 chain nodes the counts equal the measured branching children on six of seven
  bounded lines (c47: 3 against 5), with or without the fix.
- At tai60b's depth-9 node the fix moves c30 off the top (count 4, rank 25) and r8, the best line measured, to
  3rd, but c30's count is still 4 against 12. Eight of c30's twelve branching children are less than 110k short
  of the cutoff, at a slack of 326k: near misses of the bound, which a count of completions within the slack
  does not see. That part of the miss is the threshold, not the search.
- At all three tai60b nodes the fewest predicted survivors (2) belong to r46 and c10, which nobody has bounded.

## Not done

The fix in the C signal (`lsla.C`), a whole tree with it, any tuning of the 100 kicks, and the chain-model
tie-break. r46 at the tai60b nodes is unbounded.
