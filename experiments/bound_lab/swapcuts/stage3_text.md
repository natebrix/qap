Driver `swapdrv.C` on dnncuts' ADMM defaults (PRSM, degree 4, `-T 1`), cold first solve of 300
iterations (wil50, tai50b: 400). All bounds below are certified values; "plain" and "T" arms
run in the same driver, same rounds, so the comparisons are at equal iterations.

### 3a. Violations at converged plain iterates [m]

| node | d | m | cuts | violated | sum of violations (cost x probability) | share on one fixed facility | log |
|---|---|---|---|---|---|---|---|
| tho40 gate `26:10,8:0,39:32,38:7` | 4 | 36 | 5,184 | 192 (3.7%) | 1,413 | 52% on 38@7, 30% on 8@0 | `logs/s3_tho40_d4_ex.log` |
| tho40 d4_d `26:10,8:0,39:32,30:7` | 4 | 36 | 5,184 | 181 (3.5%) | 1,346 | 88% on 8@0 | `logs/s3_tho40_d4d_cp1400.log` |
| the same, after 800 / 1,400 CP iterations of the triangle phase | | | | 314 / 299 | 1,505 / 1,433 | 85% on 8@0 | `logs/s3_tho40_d4d_cp.log`, `..._cp1400.log` |
| tho40 d5_d `26:10,8:0,2:32,6:39,33:24` | 5 | 35 | 6,125 | 130 (2.1%) | 810 | 98% on 8@0 | `logs/s3_tho40_d5d_jTX.log` |
| tho40 path, depth 12 (converged, 2,500 its) | 12 | 28 | 9,408 | 333 | 36 | - | `logs/s3_tho40_path_d12_conv_ex.log` |
| wil50 node 13 `21:28,8:41,28:19` | 3 | 47 | 6,627 | 144 (2.2%) | **15** | 60% on 28@19 | `logs/s3_wil50_n13_jTXm.log` |
@@N10ROW@@
@@TAIROW@@

A violated cut reads "given j at b (probability y), exchanging j with the fixed i gains on
average". The largest at tho40 d4_d: fixed 8@0 with 2@8, y = 0.17, conditional gain -1,367;
with 35@8, y = 0.09, -1,235; with 24@8, y = 0.14, -598. In the subtree `26:10,8:0` nearly all
of it sits on facility 8 at corner location 0: the relaxation says 8 would rather trade places
with whatever lands next to it. The triangle phase does not remove these violations (they
grow slightly), so the two families are complementary in principle. On wil50 the violations
are a hundred times smaller: 15 units in all against a gap of 219.

The elimination test (`-elim`: x_jb = 0 when even the best assignment of the rest leaves the
exchange gain negative, one LAP per cut) fires on 0 of 1,296 free assignments at the two depth-4
tho40 nodes and on 1 of 1,156 at depth 6: nothing at these depths (it fires on 16 of 64 at
nug12 depth 4).

### 3b. Exchange cuts alone against the plain solve, tho40 gate node [m]

`logs/s3_tho40_d4_plain.log`, `logs/s3_tho40_d4_ex.log` (target best + 50),
`logs/s3_tho40_gate_ex400.log` (best + 400): 300 iterations, then 8 warm rounds of 150.

| total iterations | 300 | 450 | 600 | 750 | 900 | 1,050 | 1,350 | 1,500 |
|---|---|---|---|---|---|---|---|---|
| plain | 239,727.7 | 240,072.7 | 240,076.9 | 240,080.7 | 240,083.5 | 240,084.4 | 240,085.8 | 240,085.8 |
| exchange cuts (best + 50) | | 240,114.9 | 240,155.5 | 240,186.0 | 240,214.5 | 240,225.0 | 240,235.4 | 240,233.9 |
| exchange cuts (best + 400) | | 240,109.7 | 240,127.7 | 240,209.8 | 239,309.2 | 239,850.5 | 240,226.9 | 240,221.5 |

**+150 over the converged plain value** (+140 with the larger target, which overshoots once).
About 10% of the violation mass at the plain iterate. The triangle CP phase gains +1,327 on
the same node and fathoms it 120 iterations after the first solve (`../tho40_nodes/`).

### 3c. With the triangles, tho40 far-miss nodes [m]

*Joint Polyak rounds* (`-joint`: the library's `dnn_cut_rounds` on the regional triangles,
1,000 new triangles a round, target best + 400; arm TX steps the exchange multipliers between
rounds through the hook, at most 0.5 per multiplier per step). Table by `s3_summary.py`
(`logs/s3_summary.txt`); certified bound after each round, shortfall to 240,514.05 in brackets.

| total its | d4_d: T | d4_d: TX | TX - T | d5_d: T | d5_d: TX | TX - T |
|---|---|---|---|---|---|---|
| 300 | 238,237.8 | 238,237.8 | | 237,994.8 | 237,994.8 | |
| 450 | 238,759.8 | 238,902.3 | +142 | 238,551.4 | 238,631.6 | +80 |
| 600 | 239,134.6 | 239,331.7 | +197 | 238,891.2 | 239,040.1 | +149 |
| 750 | 239,511.7 | 239,636.5 | +125 | 239,251.8 | 239,380.9 | +129 |
| 900 | 239,874.9 | 240,024.1 | +149 | 239,600.0 | 239,694.0 | +94 |
| 1,050 | 240,163.9 | 240,309.2 | +145 | 239,904.6 | 240,039.3 | +135 |
| 1,200 | 240,251.0 | 240,380.2 | +129 | 240,052.6 | 240,123.2 | +71 |
| 1,350 | 240,318.1 (196 short) | **240,535.1 fathomed** | +217 | 239,914.4 | 240,203.0 | +289 |
| 1,500 | 240,288.3 | - | | 240,094.2 (420 short) | 240,204.2 (310 short) | +110 |

Logs `logs/s3_tho40_d4d_jTd400.log`, `..._d4d_jTXm.log`, `..._d5d_jT.log`, `..._d5d_jTXm.log`.
**+70 to +290 at equal iterations, median +135**: the exchange cuts keep, on top of the
triangles, the +150 they give alone. At d4_d that turns a miss into a fathom (the triangles
alone end 196 short after 1,500 iterations; the CP form alone is 214 short after 1,700,
`logs/s3_tho40_d4d_cp1400.log`). The same pair without the step cap: d4_d fathoms at 1,350 too
(`..._d4d_jTXd400.log`, 240,522.7); d5_d leads by +105, +93, +194 and then a Polyak step ten
times the others drops the bound 3,600 (`..._d5d_jTX.log`): the multipliers need the cap.
@@D4EROW@@

*Exchange rounds first, then the CP form with the exchange term frozen*
(`logs/s3_tho40_d4d_xcp.log` against `logs/s3_tho40_d4d_cp1400.log`), d4_d:

| CP iteration | 100 | 200 | 400 | 600 | 800 | 950 | 1,400 |
|---|---|---|---|---|---|---|---|
| CP alone (300 + CP) | 238,831.0 | 239,172.5 | 239,583.7 | 239,839.1 | 240,008.6 | 240,105.2 | 240,299.8 |
| 3 exchange rounds, then CP (750 + CP) | 239,019.4 | 239,350.1 | 239,753.5 | 239,996.8 | 240,164.3 | 240,256.8 | |
| difference at equal CP iterations | +188 | +178 | +170 | +158 | +156 | +152 | |

At equal CP iterations +150 to +190; at equal total iterations (1,700) -43, because the 450
iterations spent finding the multipliers are worth about +195 to the CP phase alone at that
point. The gain is a level shift of about 150, not a faster climb: it pays only if the
multipliers come without iterations of their own (inside the phase, as in the joint rounds,
or inherited from the parent).

Side observation, not this agent's topic [m]: at these two far-miss nodes the rounds form with
`-delta 400` is ahead of the CP form at equal iterations (d4_d at 1,350 total: 240,318 against
about 240,160; d5_d at 1,200: 240,053 against 239,780 at 1,100 in `../tho40_nodes/`).

### 3d. wil50 [m]

Node 13 (`21:28,8:41,28:19`, `wil50.dat`, grid 5:10 on side 0, `-u 48816`), a shortened form
of the closed studies' protocol (400 iterations, 5 rounds of 200 at best + 50 uncapped, quiet
round of 300; 1,700 iterations per arm), `logs/s3_wil50_n13_jT.log`, `..._jTXm.log`:

| total its | 400 | 600 | 800 | 1,000 | 1,200 | 1,400 | quiet (+300) |
|---|---|---|---|---|---|---|---|
| T (triangles) | 48,596.8 | 48,656.0 | 48,681.5 | 48,697.6 | 48,709.2 | 48,716.3 | 48,729.6 |
| TX (+ exchange cuts) | 48,596.8 | 48,656.8 | 48,685.1 | 48,692.5 | 48,712.4 | 48,714.0 | 48,731.0 |
| TX - T | | +0.9 | +3.5 | -5.1 | +3.1 | -2.3 | **+1.4** |

**+1.4 at equal iterations, inside the round-to-round scatter of 5**, against a bar of +10.
With the uncapped Polyak step the bound falls 331 in the first round
(`logs/s3_wil50_n13_jTX_overshoot.log`: 15 units of violation cannot carry a step toward
best + 50). @@N10TEXT@@

@@TAITEXT@@

### 3e. Validity [m]

- Brute force, pairwise cuts and fixed-pair cuts: 84 nodes of seven n = 12 instances (depths
  2-4), 1,680 certified values against the enumerated optimum over exchange-optimal completions:
  max(certified - v_ex) = 0.0000 (`logs/validate_small.out`). 5 to 10 of every 12 nodes end
  **above the node's own optimum** v_node, which no inequality valid for all permutations can
  do; the relaxation captures a median of 21% (nug12) to 89% (scr12) of the exact rise.
- All families (pairwise, fixed-pair, rotations) on 112 on-path nodes of the same instances,
  two Polyak targets: 1,568 certified values, none above the optimum
  (`logs/validate_path_small.out`).
- tho40, the path to the best known solution (`../tho40_nodes/scripts/path.py`, cutoff 241,500
  so nothing stops early): exchange cuts alone at depths 4 and 8 end at 235,069 and 238,793
  (`logs/s3_tho40_path_d4_ex.log`, `..._d8_ex.log`); joint with the triangles at depth 12:
  240,255.8 (`..._path_d12_jTX.log`); from a converged plain solve at depth 12 (240,515.40 at
  2,500 iterations) six rounds with all families at a target of best + 2 stay at 240,437 to
  240,457, and at best + 400 fall to 128,446 and climb back to 202,793
  (`..._path_d12_conv_exrot_d2.log`, `..._path_d12_conv_ex.log`). **No certified value above
  240,516 in any path run.** @@PATH8@@
