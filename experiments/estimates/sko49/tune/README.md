# sko49: an estimate with sko42's best settings (cloud box, from Nov 12)

sko49 gets the settings that won on sko42 (`../../sko42/tune/README.md`), on its own
7 x 7 grid, and the same sampling plan as sko42. Below depth 3, sampled probes
(`../../probe.py`) stand in for the complete subtrees, which are too large here.

**The result (Nov 13).** Cold index, iso and fit models:
- under C: 60,000-66,000 cloud core-hours (95%: 23,000-150,000);
- under Cw (C with `cut.zone = 0.018`): 37,000-42,000 (95%: 17,000-66,000).

Cw / C is 0.62-0.63 in paired resampling, cheaper in 90-95% of draws. See "The
estimate under Cw" at the end.

## Settings: arm C (`C.par`, `C.sh`, `C.golden`)

`C.par` is `../sko49.par` (wil50's full run with `cut.grid = 7:7`, `cut.side = 1`)
with `cut.zone = 0.010`:
- **The grid.** `sko49t.dat`'s second matrix is the 7 x 7 rectilinear grid,
  row-major with 7 columns (checked again Nov 12, every entry). So the regions
  are cut on side 1, and `latest/SYM/sko49t.sym` (10 depth-1 orbits) is the
  symmetry.
- **The zone.** On sko42, zone 0.010 beat 0.0075 by 21-39% on five complete
  subtrees. On the 7 x 7 grid a cut phase gains +126 to +149 at depth 3
  (`../../../sko49_feasibility/NOTES.md` section 4), against 80-100 on sko42's
  6 x 7. So 0.010 is the same zone in units of the phase's gain on both
  instances, about 1.7 gains: 234 short here, 158 on sko42.
- **Everything else** is `sko49.par` as it stood on Nov 12. `C.golden`
  differed from `../settings.golden` only in the zone, the hashes, and the
  pool-check diagnostic, which `../settings.sh` turns on and sko42's arms leave
  off. From Nov 13, `sko49.par` is Cw (`cut.zone = 0.018`). Its golden comes
  from HEAD's build, which adds `dnn.giveup_warm`, `cut.giveup` and
  `QAP_DNN_TRACE_HIST` at their defaults, and renames `DNN_ADMM_CROSS_FIRST`
  to `dnn.cross_first`.

## Plan

| step | what | expected cloud-box time |
|---|---|---|
| stage A | every node at depths 0-2 (about 315 after symmetry), the depth-3 frontier | 3-4 h |
| stage B | 250 frontier nodes to the give-up point (`--converge 0`: the estimate indexes by give-up shortfall, as on sko42) | 2.5-3 h |
| anchors | the subtree cost by give-up shortfall: complete subtrees (`../../subtree.sh`) where they are affordable, probes (`../../probe.py`) across the rest, one node with both as the probe's check on sko49 | 6-10 h |
| estimate | `../../sko42/tune/estimate.py` with sko49's table, top and frontier | minutes |

## Stage A under C (`sampling_C/`, Nov 12)

```
SAMPLING_SETTINGS=experiments/estimates/sko49/tune/C.sh SAMPLING_PARAMS=experiments/estimates/sko49/tune/C.par \
  python3 experiments/wil50_feasibility/sampling/stage_a.py --run <dir> --slots 4 --cpus 0,1,2,3 --tau 0.025
```

The binary is the Nov 11 18:22 build (sources `dc5e5d7`), the one every sko42
arm ran on. The laptop's `dnn.giveup_warm` commits (`75fb712`, `714b199`, off by
default) are newer than it.

- **Cost.** 344 nodes and 68,096 CPU s (18.9 core-hours) on this box, 4.8 h
  on 4 slots: sko42's stage A under C took 58,946 s for 451 nodes.
- **Depth 1.** The 10 nodes are 451-642 short (median 608) after 300
  iterations, 115 s each, with 27 or 48 children. sko42's 12 were 231-393.
- **Depth 2.** 333 nodes, 1 fathomed in its first solve. The 54 that ended
  within the zone (234 short) got a cut phase:
  - the 19 from 2-156 short fathomed;
  - the 35 from 158-234 stalled 11-128 short, a gain of 130-160.
  - So 20 of 333 fathom (6%), against 165 of 438 (38%) on sko42.
- **The frontier.** 13,682 depth-3 nodes. Their inherited shortfall (the
  parent's best bound) has median 358 and 10th percentile 104. sko42's 10,704
  had median 181 and 10th percentile 25.
  - In zones that is 1.5 against 1.15.
  - In cut-phase gains it is 2.6 against 2.0.
- **The summary's wording.** `summarize_a.py` writes wil50's fixed text in
  three places, and none of them applies here: "the cut zone is 244", the
  200-244 / 244-300 buckets, and the expected 15 orbits and 635 depth-2 nodes.
  sko49's zone is 234 (0.010 x 23,384), and its symmetry gives 10 depth-1
  orbits.

## Stage B under C (`sampling_C/stageB_*`, Nov 12)

```
SAMPLING_SETTINGS=... SAMPLING_PARAMS=... python3 experiments/wil50_feasibility/sampling/stage_b.py \
  --stage-a <stage A dir> --run <dir> --n 250 --converge 0 --project-converge 40 --slots 4 --cpus 0,1,2,3
```

250 frontier nodes, the smallest seeded hashes, each solved warm to the give-up
point with no cut phase at depth 3. This took 21,542 CPU s, 1.5 h on 4 slots.

| | sko49 under C | sko42 under C |
|---|---|---|
| fathomed in the first solve | 52 of 250 (21%) | 140 of 250 (56%) |
| the misses' give-up shortfall: min / 10% / median / 90% / max | 10 / 45 / 164 / 303 / 396 | 5 / 14 / 68 / 129 / 189 |
| misses beyond the zone (234 / 158) | 57 (23% of the sample) | 2 |
| CPU per hopeless give-up | 104 s | 42 s |

Measured in zones, sko42's farthest miss is 1.2 zones out. Here 57 of the 250
are beyond 1.0, reaching 1.7. On sko42 a cut phase fathomed every depth-2 node
within 0.77 zones; here, within 0.67 (156). Those 57 nodes get no cut phase at
depth 3, and they are where the cost will be.

## Anchors (`anchors_C/`, Nov 12)

Nine probes ran on 3 pinned slots: `../../probe.py --k 16 --beta 0.025`, cold
from `-F` as sko42's anchors were. Each is at a quantile of the 198 misses'
give-up shortfalls:
- 45, 97, 164, 188, 232, 268, 303, 326 and 396;
- that is the 10th, 30th, 50th, 60th, 70th, 80th, 90th and 95th percentiles,
  and the maximum.

Why probes: the record's two sko49 depth-3 subtrees, 263-304 short converged,
were estimated at 47,000 laptop CPU s each under the Nov 5 settings
(`../../../sko49_feasibility/vs_wil50.md`), about 25 core-hours on this box.
sko42's largest was 1.3. A probe draws 16 of a level's open nodes, bounds
them, and weights them up, so it costs 1-80 nodes whatever the subtree's size.

| probe | stage B short (warm) | first solve short (cold) | phase at the node | levels | DNN nodes (est.) | DNN CPU s (est., SE) | the probe's own CPU s |
|---|---|---|---|---|---|---|---|
| Cp45 | 44.9 | 56.6 | fathomed | 3 | 1 | 172 | 172 |
| Cp97 | 96.8 | 105.9 | fathomed | 3 | 1 | 302 | 302 |
| Cp164 | 163.9 | 171.9 | to 10.4 | 3-4 | 47 | 1,536 (22) | 1,299 |
| Cp188 | 188.0 | 205.8 | to 55.7 | 3-4 | 47 | 1,671 (575) | 955 |
| Cp232 | 231.5 | 245.5 | none | 3-6 | 532 | 17,471 (10,012) | 3,872 |
| Cp268 | 267.6 | 280.1 | none | 3-5 | 176 | 10,082 (2,662) | 3,335 |
| Cp303 | 303.1 | 318.3 | none | 3-6 | 1,177 | 48,408 (33,036) | 7,527 |
| Cp326 | 325.5 | 349.1 | none | 3-7 | 5,514 | 170,185 (149,685) | 13,127 |
| Cp396 | 395.7 | 408.7 | none | 3-8 | 65,930 | 1,163,004 (790,289) | 11,957 |

- **The checks.** Two complete subtrees (`../../subtree.sh tune/C.sh`, one slot):
  - **s188** (`Cs188/`): 47 DNN nodes and 1,328 s against the probe's
    1,671 ± 575 s, within 0.6 SE. Every child of s188 fathoms, so this checks
    the probe's CPU, not its node count, which is exact by construction.
  - **s232** (`Cs232/`, 3.9 h): 271 DNN nodes at depths 3-6 (1, 46, 180, 44)
    and 13,785 s, against the probe's 532 ± 485 nodes and 17,471 ± 10,012 s,
    within 0.6 SE and 0.4 SE. 4 of the 46 depth-4 children stay open, where
    the probe drew 1 in 16; one depth-5 node stays open, where the probe's
    one open draw at depth 5 is weighted to about 8.
  - Read with s232 in place of its probe, the estimate below drops by 1-4%.
- **The cold start.** A probe's root is solved cold from `-F`: 8-24 short
  further than stage B's warm give-up. That matters at the zone's edge. s232
  is 231.5 short warm, inside the zone, but its cold 245.5 is outside, so it
  gets no phase. The estimate below indexes the anchors both ways:
  - by the warm shortfall (sko42's method), which overstates;
  - by the cold one, which matches the node the probe actually bounded.
- **Where the subtrees blow up.** Beyond the zone, a node's 46 children come
  out 85-300 short in their first solve. Below the five largest probes' roots:
  - 77 of the 86 cut phases on nodes within 135 short fathomed;
  - none of the 27 on nodes from 135 to 234 did: they stalled 9-121 short,
    after 300-800 s, with 44-45 children each;
  - 7 drawn nodes beyond 234 got no phase.
  
  Under s396, 12 of the 16 children drawn at depth 4 stay open, and the tree
  reaches depth 8.

## The estimate with nine probes (Nov 12)

```
python3 experiments/estimates/sko42/tune/estimate.py --table experiments/estimates/sko49/tune/sampling_C/stageB_summary.md \
  --top 344 --top-cpu 68096 --front 13682 anchors_C/Cp45=56.6 ... anchors_C/Cp396=408.7    # or =44.9 ... =395.7
```

| anchors indexed by | model | DNN nodes | CPU s (95%, stage B's sample only) | cloud core-hours |
|---|---|---|---|---|
| cold first solve | interp | 7.4e6 | 2.0e8 (1.0-3.3e8) | 56,000 |
| cold first solve | fit (x12.2 per 100 short) | 9.6e6 | 2.2e8 (1.4-3.2e8) | 61,000 |
| warm give-up | interp | 1.4e7 | 3.5e8 (1.9-5.5e8) | 98,000 |
| warm give-up | fit (x12.9 per 100 short) | 1.6e7 | 3.3e8 (2.1-4.8e8) | 91,000 |

- **Where the cost is.** The 9 sampled nodes beyond 330 short carry 69% of
  it, and the 19 beyond 305 carry 81-84%. The 193 below 234 carry 4%.
- **The anchors' own error is not in the interval.** Cp326 and Cp396, which
  set the tail's 69%, are ± 88% and ± 68%.
- **Against the record.** The Nov 6 estimate was 6.4e8 laptop CPU s, about
  3.6e5 cloud core-hours. This one is 3.6-6.3x less, and still out of reach of the
  laptop: at its 2x speed, 28,000-49,000 laptop core-hours.

**Replicates (running from Nov 12 17:51).** Six more probes, on other sampled
nodes in the tail, with the same k, beta and seed:
- `probes2`, on 3 slots: Cr311 `0:0,40:22,41:10`, Cr338 `30:25,41:8,46:35`,
  Cr368 `30:33,41:8,46:14` and Cr385 `12:0,41:2,46:1`;
- `probes3`, on one slot once s232 finishes: Cr283 `0:0,40:23,41:3` and Cr350
  `9:9,40:23,41:17`.

All six run the same binary as the nine probes (md5 `ee4c6ebd…`). The
"binary" a summary names is the repository's commit at launch, not the build.

**How the replicates enter (decided Nov 12 20:40, with only Cr338 in).**
- Cr338 is 350.2 short cold and cost 27,433 CPU s (SE 5,785 s). Cp326 is
  349.1 short and cost 170,185 s. Two nodes at the same shortfall differ 6x.
- interp runs through both points, so the curve from 350 to 409 short drops
  about 2.5x, and the total with it: 55,000 to 38,000 core-hours.
- `estimate_se.py` now also prints an iso model: interp through the anchors'
  isotonic fit (pool adjacent violators, so an anchor dearer than the next is
  averaged with it). Its prior is that the mean cost does not fall as the
  shortfall grows. With the nine probes it gives 56,600 core-hours, against
  interp's 55,500.
- With Cr338 added (cold index, Cs232 for Cp232):

  | model | CPU s | cloud core-hours (95%, sample and anchors) |
  |---|---|---|
  | interp | 1.37e8 | 38,000 (18,000-68,000) |
  | iso | 1.71e8 | 47,000 (17,000-83,000) |
  | fit | 1.80e8 | 50,000 (18,000-87,000) |

  The anchors are given in shortfall order. The intervals' ends move a few
  percent with the order, which changes the resampling draws.

- The final estimate reads iso and fit with all fifteen anchors. interp stays
  in the table but is not the headline once anchors have replicates.

**Thirteen anchors (Nov 12 21:30; Cr283 and Cr350 still running).**
`probes2` finished; its log is `anchors_C/probe2.log`.

| replicate | warm | cold | DNN nodes | DNN CPU s (SE) | depth reached |
|---|---|---|---|---|---|
| Cr311 | 311.5 | 326.1 | 3,211 | 93,371 (69,377) | 6 |
| Cr338 | 337.7 | 350.2 | 1,016 | 27,433 (5,785) | 5 |
| Cr368 | 367.5 | 380.1 | 107,424 | 1,795,628 (1,411,379) | 8 |
| Cr385 | 385.0 | 404.7 | 11,102 | 244,929 (133,950) | 7 |

- **The scatter.** Near 350 short, Cp326 and Cr338 cost 170,185 and 27,433 s.
  From 380 to 409 short, Cr368, Cr385 and Cp396 cost 1.8e6, 2.4e5 and
  1.2e6 s. So the cost at a given shortfall spreads over 5-7x from node to
  node, not just from the probes' own error.
- **iso's blocks** (cold): Cs232 and Cp268 pool at 262.8 short, Cp326 and
  Cr338 at 349.6, and Cr368 and Cr385 at 392.4 (1.02e6 s).

| index | model | DNN nodes | CPU s | cloud core-hours (95%, sample and anchors) |
|---|---|---|---|---|
| cold | interp | 1.0e7 | 2.38e8 | 66,000 (25,000-142,000) |
| cold | iso | 9.2e6 | 2.20e8 | 61,000 (23,000-150,000) |
| cold | fit | 1.1e7 | 2.58e8 | 72,000 (29,000-154,000) |
| warm | interp | 1.7e7 | 3.60e8 | 100,000 (34,000-236,000) |
| warm | iso | 1.6e7 | 3.51e8 | 97,000 (35,000-249,000) |
| warm | fit | 1.7e7 | 3.76e8 | 105,000 (42,000-226,000) |

With 13 anchors interp is back near iso, because the replicates are no longer
all on one side.

**Fourteen anchors (Nov 12 23:10; Cr350 still running).** Cr283 is 283.3 short
warm and 296.4 cold. It has 591 DNN nodes and cost 20,595 CPU s (SE 5,916), to
depth 5. Cr350 is 349.8 short warm and 360.8 cold.

| index | model | CPU s | cloud core-hours (95%, sample and anchors) |
|---|---|---|---|
| cold | interp | 2.38e8 | 66,000 (24,000-157,000) |
| cold | iso | 2.14e8 | 59,000 (22,000-166,000) |
| cold | fit | 2.42e8 | 67,000 (27,000-143,000) |
| warm | interp | 3.60e8 | 100,000 (34,000-242,000) |
| warm | iso | 3.45e8 | 96,000 (35,000-272,000) |
| warm | fit | 3.52e8 | 98,000 (38,000-211,000) |

## The estimate with fifteen anchors (Nov 12 23:27)

Cr350 finished last: 349.8 short warm and 360.8 cold. It has 6,190 DNN nodes
and cost 195,575 CPU s (SE 87,762), to depth 6. Its log is
`anchors_C/probe3.log` (with Cr283's).

The anchors, cold first solve: Cp45 56.6, Cp97 105.9, Cp164 171.9,
Cp188 205.8, Cs232 245.5, Cp268 280.1, Cr283 296.4, Cp303 318.3, Cr311 326.1,
Cp326 349.1, Cr338 350.2, Cr350 360.8, Cr368 380.1, Cr385 404.7 and
Cp396 408.7. The warm give-ups are in the tables above.

```
python3 experiments/estimates/sko49/tune/estimate_se.py --table experiments/estimates/sko49/tune/sampling_C/stageB_summary.md \
  --top 344 --top-cpu 68096 --front 13682 anchors_C/Cp45=56.6 ... anchors_C/Cp396=408.7    # in shortfall order
```

| index | model | DNN nodes | CPU s | cloud core-hours (95%, sample and anchors) | laptop core-hours |
|---|---|---|---|---|---|
| cold | interp | 1.1e7 | 2.51e8 | 70,000 (25,000-153,000) | 35,000 |
| cold | **iso** | 9.0e6 | 2.16e8 | **60,000** (23,000-150,000) | 30,000 |
| cold | **fit** (x12.1 per 100 short) | 9.4e6 | 2.39e8 | **66,000** (28,000-131,000) | 33,000 |
| warm | interp | 1.7e7 | 3.63e8 | 101,000 (35,000-245,000) | 50,000 |
| warm | iso | 1.5e7 | 3.35e8 | 93,000 (35,000-251,000) | 46,000 |
| warm | fit | 1.5e7 | 3.44e8 | 96,000 (40,000-192,000) | 48,000 |

- **The headline.** The cold index with iso and fit gives 60,000-66,000 cloud
  core-hours. The joint 95% range is 23,000-150,000.
- **Against the nine probes.** That estimate was 56,000-61,000. The
  replicates moved the point little, but they show the tail's scatter, which
  is now in the range.
- **Where the cost is** (cold; iso, then fit):
  - the 9 sampled nodes beyond 330 short carry 73% and 53%;
  - the 19 beyond 305, 84% and 71%;
  - the 38 from 234 to 305, 13% and 24%;
  - the 193 up to 234, 3% and 5%.
- **The laptop.** It is about 2x the cloud box per DNN iteration, so 30,000
  to 33,000 laptop core-hours: 5,000-5,500 hours, or 210-230 days, on its 6
  P-cores.

## The wider-zone test (arm Cw, from Nov 12 18:46)

`Cw.par` and `Cw.sh` are C with `cut.zone = 0.018`, which is 421 short. Under
it, every depth-3 node in stage B's sample runs a cut phase: the farthest is
409 short in its cold first solve. Under C the 57 sampled nodes beyond 234 run
none, and they hold 81-84% of the estimate. The question is whether a phase at
such a node pays for itself:
- it cannot fathom the node (a phase gains about 130-160);
- its children start from its pool;
- children from 234 to 421 short also get a phase.

`Cw.golden` comes from the run binary, the Nov 11 18:22 build (65 settings),
as `C.golden` does.

- **The prior.** On wil50 (Nov 5, `../../../wil50_feasibility/tail_probe/README.md`,
  "fullwide against full"), the same idea on a 455-short tail node was a wash:
  the root phase cost 800 iterations, level 4 was 8% cheaper and level 5 3%
  dearer, with the same fathoms. That was the rounds phase at zone 0.005
  against 0.010. sko49's CP phase on the 7 x 7 grid gains more.
- **The runs.** Probes with the same k, beta and seed on the two nodes that
  set the tail, against their C probes:
  - Zp396 (`9:30,40:22,41:17`, against Cp396) on one slot from 18:46;
  - Zp326 (`0:0,40:9,41:3`, against Cp326) from when the replicates on 3
    slots finish.
- **What to read.** The probe totals have SEs of 70-90%, so only a factor of
  about 2.5 between Cw and C would show. Level 4 is less noisy: compare its
  fathom fraction and its children's shortfalls, 16 drawn of 46 in each.
- **Slots.** probes3 (Cr283, Cr350) was drained to free a slot for Zp396. It
  resumes on 2 slots when the replicates finish.

**Zp396 through level 4 (Nov 12 20:17).**
- **The root.** Its phase took it from 408.7 to 279.0 short, in 245 s on top
  of the 89 s solve. As expected, that is far short of fathoming, and all 46
  children stay open.
- **The draws differ.** Each level's random start is seeded by the probe's
  name, and the root's phase also reorders the children by predicted bound.
  Still, 10 of the 16 children drawn at depth 4 are the same nodes in both
  probes (`anchors_Cw/cmp_level.py`). On those 10:
  - each child's first solve starts 60-95 closer than under C;
  - the same 4 fathom in both;
  - C's 3 children without a phase (246-270 short) get one under Cw. It
    stalls at 125-152 short.
  - the 7 phases that run in both cost 22-43% less under Cw. The 3 of those
    that stall end 1-8 closer.
  - all 10 cost 2,855 s under Cw and 3,162 s under C, 10% less, the 3 extra
    phases included.
- **The other 6 draws.** Cw's 6 include two near-fathom children (7 and 63
  short) that fathom; C's 6 fathom none. So level 4's 6/16 against 4/16 is the
  draw, not the setting.
- **What it means for the rest.** Each probe's open depth-4 children end where
  they hand off to their own children:
  - under Cw, all 10 ran a phase and end 18-153 short;
  - under C, 7 ran a phase and end 18-103, but 5 had none and stay at 246-303.

  If Cw pays, it shows at levels 5-6, which are still running. Through level 4
  (the root and the 16 draws), Cw costs 4,991 s and C 5,265 s.

**Zp396 complete (Nov 12 21:32; `anchors_Cw/Zp396`).**

| level | C nodes (est.) | C CPU s per node | C CPU s | Cw nodes (est.) | Cw CPU s per node | Cw CPU s |
|---|---|---|---|---|---|---|
| 3 | 1 | 91 | 91 | 1 | 334 | 334 |
| 4 | 46 | 323 | 14,874 | 46 | 291 | 13,388 |
| 5 | 680 | 306 | 208,330 | 1,273 | 89 | 112,889 |
| 6 | 3,149 | 71 | 223,622 | 6,700 | 48 | 323,781 |
| 7 | 18,819 | 27 | 504,884 | 19,110 | 12 | 223,733 |
| 8 | 43,235 | 4.9 | 211,202 | - | - | - |
| all | 65,930 | | 1,163,004 (SE 790,289) | 27,130 | | 674,126 (SE 457,027) |

- **The total.** Cw is 0.58x, but that is only about 0.5 SE of the
  difference.
- **The mechanism is clear at level 5.**
  - Under Cw, 10 of the 16 drawn nodes fathom at their first solve (median
    first solve 0 short, range -9 to 119), and 6 run a phase.
  - Under C, 3 fathom at their first solve (median 72 short, range -15 to
    229), and 13 run a phase.

  Every Cw node at depth 4 ran a phase, so its children start from a better
  pool. Per node, levels 5-7 cost 1.5-3.5x less, and the tree stops at depth 7
  instead of 8.
- **Node counts are noisier.** Cw's estimate has more nodes at levels 5-6.
- **What it would mean.** If the tail nodes all came out like this, the 81-84%
  of the estimate beyond 305 short would drop by about 40%, and the total by
  about a third. That is one pair. Zp326 and Zr368 test it.

**A third pair (from Nov 12 21:34).** Zr368 is Cw on Cr368's node
(`30:33,41:8,46:14`, 380.1 short cold), the dearest anchor at 1.8e6 s. It runs
on cpu 3, which Zp396 freed. Its log is `zoneC`.

**Zp326 complete (Nov 12 23:10; `anchors_Cw/Zp326`).**

| level | C nodes (est.) | C CPU s per node | C CPU s | Cw nodes (est.) | Cw CPU s per node | Cw CPU s |
|---|---|---|---|---|---|---|
| 3 | 1 | 89 | 89 | 1 | 335 | 335 |
| 4 | 46 | 422 | 19,417 | 46 | 285 | 13,116 |
| 5 | 135 | 295 | 39,905 | 1,056 | 20 | 20,982 |
| 6 | 1,588 | 49 | 77,305 | 746 | 7.0 | 5,231 |
| 7 | 3,744 | 8.9 | 33,469 | - | - | - |
| all | 5,514 | | 170,185 (SE 149,685) | 1,850 | | 39,663 (SE 19,782) |

- **The total.** Cw is 0.23x, about 0.9 SE of the difference.
- **The root.** Its phase took it from 349.1 to 219.6 short.
- **Level 4.** Only 4 of the 16 draws are the same nodes in both. On those 4,
  the first solves start 0-68 closer under Cw, and the 2 phases that run in
  both cost 16-24% less. The depth-4 first solves have a median of 50 short
  under Cw and 133 under C.
- **Level 5.**
  - Under Cw, 15 of the 16 drawn nodes fathom at their first solve (median -1
    short, range -108 to 45).
  - Under C, none do: all 16 run a phase (median 74 short, range 15 to 211),
    and 5 stay open.

  Per node, levels 5-6 cost 7-15x less, and the tree stops at depth 6 instead
  of 7.

**Two pairs so far.** Cw costs 0.58x and 0.23x, by the same mechanism: every
depth-4 parent runs a phase, so its children fathom at their first solve. The
extra phases are cheap. A phase at a depth-3 node costs about 245 s. Under
Cw, the 23% of the 13,682 frontier nodes that are beyond 234 short would run
one, about 7.7e5 CPU s in all, or 0.4% of the estimate.

**A Cw estimate (`zoneD`, from Nov 12 23:12).** The first solve, which indexes
the anchors, does not depend on the zone. Under Cw, a node behaves as under C
only if it and all its descendants are within 234 short. So the anchors up to
206 short (cold) carry over: no node drawn in their subtrees is farther than
its root. Only those beyond 234 need Cw probes:
- Zp326, Zp396 and Zr368 (`zoneA`-`zoneC`);
- in `zoneD`, the other eight: Zr385, Zr350, Zr338, Zr311, Zp303, Zr283, Zp268
  and Zp232 (on s232's node), the tail first.

`zoneD` starts on cpu 0. It drains and resumes on more cpus as Cr350 (cpus 1
and 2) and Zr368 (cpu 3) finish, so the box runs at most 4 solvers.

**Zr368 complete (Nov 12 23:46; `anchors_Cw/Zr368`).**

| level | C nodes (est.) | C CPU s per node | C CPU s | Cw nodes (est.) | Cw CPU s per node | Cw CPU s |
|---|---|---|---|---|---|---|
| 3 | 1 | 85 | 85 | 1 | 311 | 311 |
| 4 | 46 | 299 | 13,748 | 46 | 287 | 13,192 |
| 5 | 391 | 255 | 99,814 | 1,346 | 70 | 94,292 |
| 6 | 6,721 | 52 | 348,282 | 11,559 | 17 | 199,524 |
| 7 | 19,343 | 30 | 576,909 | 14,717 | 5.9 | 86,994 |
| 8 | 80,923 | 9.4 | 756,791 | - | - | - |
| all | 107,424 | | 1,795,628 (SE 1,411,379) | 27,669 | | 394,314 (SE 228,215) |

- **The total.** Cw is 0.22x, about 1.0 SE of the difference.
- **The root.** Its phase took it from 380.1 to 268.6 short.
- **Level 4.** The first solves have a median of 73 short under Cw and 128
  under C. 4 of the 16 draws are the same nodes in both:
  - one has no phase under C (252 short) and stalls at 161 under Cw (256 s);
  - the two phases that run in both cost 18-22% less under Cw and end at the
    same shortfall;
  - one fathoms at its first solve in both.
- **Level 5.** Under Cw, 10 of 16 fathom at their first solve (median -0.5
  short), against 3 of 16 under C (median 71). Per node, levels 5-7 cost
  3.0-5.1x less, and the tree stops at depth 7 instead of 8.

**Three pairs.**

| node (cold) | C CPU s | Cw CPU s | Cw / C | difference in SEs | level 5 fathomed at the first solve, Cw against C |
|---|---|---|---|---|---|
| 0:0,40:9,41:3 (349.1) | 170,185 | 39,663 | 0.23 | 0.9 | 15 / 16 against 0 / 16 |
| 30:33,41:8,46:14 (380.1) | 1,795,628 | 394,314 | 0.22 | 1.0 | 10 / 16 against 3 / 16 |
| 9:30,40:22,41:17 (408.7) | 1,163,004 | 674,126 | 0.58 | 0.5 | 10 / 16 against 3 / 16 |
| all three | 3,128,817 | 1,108,103 | 0.35 | 1.2 | 35 / 48 against 6 / 48 |

- **No total is significant alone.** Probe totals are dominated by the noisy
  node counts at the deep levels.
- **The mechanism is consistent and strong.** At level 5, 35 of 48 draws
  fathom at their first solve under Cw, against 6 of 48 under C. Each probe's
  per-node cost at levels 5-7 is lower under Cw.
- **The Cw estimate** (zoneD) prices this across the whole tail; see below.

**The other eight (`zoneD`, done Nov 13 02:31; `anchors_Cw/Z*`, log `anchors_Cw/zoneD.log`).**
The container restarted at about 00:40 UTC on Nov 13 and killed `zoneD` at
level 4. It resumed at 00:47 on cpus 0-3 with the same binary. `probe.py` kept
every finished task and reran the four cut-off ones from their parents' states.
Each task's CPU comes from its own complete run, so the results are as if
there had been no restart.

**Eleven pairs.** C's probe against Cw's on each anchor node beyond 234 short
(for 245.5 short, C's probe Cp232; the estimate uses the complete subtree
Cs232). "First-solve fathoms" counts the 16 draws at a level that fathom at
their first solve, Cw against C (`anchors_Cw/pairs.py`, which reads the
uncommitted run directories).

| node (cold short) | C CPU s | Cw CPU s | Cw / C | Cw - C in SEs | level 4 first-solve fathoms | level 5 first-solve fathoms |
|---|---|---|---|---|---|---|
| 23:0,40:8,41:10 (245.5) | 17,471 | 5,834 | 0.33 | -1.1 | 13 against 4 | 16 against 14 |
| 9:34,40:26,41:16 (280.1) | 10,082 | 3,492 | 0.35 | -2.2 | 12 against 4 | Cw stops at 4; C 16 |
| 0:0,40:23,41:3 (296.4) | 20,595 | 17,974 | 0.87 | -0.3 | 7 against 3 | 15 against 15 |
| 30:17,41:8,46:1 (318.3) | 48,408 | 36,415 | 0.75 | -0.3 | 7 against 4 | 15 against 5 |
| 0:0,40:22,41:10 (326.1) | 93,371 | 56,366 | 0.60 | -0.5 | 7 against 2 | 13 against 1 |
| 0:0,40:9,41:3 (349.1) | 170,185 | 39,663 | 0.23 | -0.9 | 5 against 2 | 15 against 0 |
| 30:25,41:8,46:35 (350.2) | 27,433 | 98,222 | 3.58 | +1.2 | 5 against 2 | 12 against 15 |
| 9:9,40:23,41:17 (360.8) | 195,575 | 267,900 | 1.37 | +0.5 | 1 against 0 | 7 against 1 |
| 30:33,41:8,46:14 (380.1) | 1,795,628 | 394,314 | 0.22 | -1.0 | 3 against 2 | 10 against 3 |
| 12:0,41:2,46:1 (404.7) | 244,929 | 382,100 | 1.56 | +0.6 | 2 against 1 | 7 against 5 |
| 9:30,40:22,41:17 (408.7) | 1,163,004 | 674,126 | 0.58 | -0.5 | 0 against 0 | 10 against 3 |
| all eleven | 3,786,681 | 1,976,406 | 0.52 | | 62 against 24 of 176 | 120 against 62 of 160 |

- **The totals.** Cw is cheaper in 8 of the 11 and dearer in 3. Only one
  difference reaches 2 SEs (280.1 short, in Cw's favour). The geometric mean
  of the ratios is 0.66.
- **The mechanism holds.**
  - At level 4, more draws fathom at their first solve under Cw in 10 of the
    11 pairs, and the other is a tie (0 and 0).
  - At level 5, more do in 8 of the 10 pairs that have one, with a tie in 1.
    Fewer do in 1: 350.2 short, 12 against 15.
- **The three dearer pairs.** Levels 3-4 cost about the same in both arms
  (15,400-16,200 CPU s under C, 14,200-15,500 under Cw), so the extra cost is
  at levels 5-6. There a few open draws set the node counts:
  - 350.2: C's level 5 fathoms all 16 and stops at depth 5. Under Cw, 2 stay
    open, and level 6 comes out at 2,130 nodes and 51,446 s.
  - 360.8: level 5 costs 121 s a node under Cw against 260 under C, but its
    node estimate is 1,281 against 542.
  - 404.7: Cw stops at depth 6, C at depth 7. Cw's level 6 comes out at
    15,431 nodes, against 3,379 under C.

## The estimate under Cw (Nov 13 02:40)

The anchors are C's four up to 206 short (cold) and the eleven Cw probes:
Cp45, Cp97, Cp164 and Cp188, then Zp232, Zp268, Zr283, Zp303, Zr311, Zp326,
Zr338, Zr350, Zr368, Zr385 and Zp396. The index values are C's: the cold
first solve does not depend on the zone, and the warm give-ups come from
stage B, which ran under C.

```
cd experiments/estimates/sko49/tune
python3 estimate_se.py --table sampling_C/stageB_summary.md --top 344 --top-cpu 68096 --front 13682 \
  anchors_C/Cp45=56.6 anchors_C/Cp97=105.9 anchors_C/Cp164=171.9 anchors_C/Cp188=205.8 \
  anchors_Cw/Zp232=245.5 anchors_Cw/Zp268=280.1 ... anchors_Cw/Zp396=408.7    # in shortfall order
```

| index | model | DNN nodes | CPU s | cloud core-hours (95%, sample and anchors) | laptop core-hours | C, same model |
|---|---|---|---|---|---|---|
| cold | interp | 6.0e6 | 1.34e8 | 37,000 (16,000-69,000) | 19,000 | 70,000 |
| cold | **iso** | 5.0e6 | 1.34e8 | **37,000** (17,000-64,000) | 19,000 | 60,000 |
| cold | **fit** (x11.0 per 100 short) | 6.2e6 | 1.51e8 | **42,000** (23,000-66,000) | 21,000 | 66,000 |
| warm | interp | 8.6e6 | 1.97e8 | 55,000 (24,000-104,000) | 27,000 | 101,000 |
| warm | iso | 7.8e6 | 2.01e8 | 56,000 (27,000-92,000) | 28,000 | 93,000 |
| warm | fit | 9.7e6 | 2.12e8 | 59,000 (32,000-93,000) | 29,000 | 96,000 |

**Cw against C, paired.** `paired.py` takes the two arms' anchors node by node.
Its interval resamples stage B's sample and the anchor nodes with the same
draws for both arms, so the scatter both share (which nodes are dear) mostly
cancels.

```
python3 paired.py --table sampling_C/stageB_summary.md --top 344 --top-cpu 68096 --front 13682 \
  anchors_C/Cp45,anchors_C/Cp45=56.6 ... anchors_C/Cs232,anchors_Cw/Zp232=245.5 ... anchors_C/Cp396,anchors_Cw/Zp396=408.7
```

| index | model | C core-hours | Cw core-hours | Cw / C (95%, paired) | draws with Cw cheaper |
|---|---|---|---|---|---|
| cold | iso | 60,080 | 37,231 | 0.62 (0.32-1.16) | 90% |
| cold | fit | 66,400 | 41,966 | 0.63 (0.40-1.07) | 95% |
| warm | iso | 92,982 | 55,891 | 0.60 (0.29-1.14) | 93% |
| warm | fit | 95,536 | 58,752 | 0.61 (0.37-1.08) | 95% |

- **One pair carries much of the gap.** Dropping any one pair (cold) leaves
  Cw / C at 0.56-0.68, except the 380.1-short pair (Cr368 against Zr368).
  Without that pair it is 0.95 (iso) and 0.81 (fit), so it carries most of
  the iso gap and half the fit gap. Cr368's own SE is 79% of its value
  (1,795,628 +- 1,411,379).
- **The tail is tighter under Cw.** Its dearest anchor is 674,126 s, against
  C's 1,795,628, and the 95% range (iso and fit) narrows from 23,000-150,000
  to 17,000-66,000.
- **Where the cost is** (cold; iso, then fit). The shares are close to C's:
  - the 9 sampled nodes beyond 330 short carry 69% and 51%;
  - the 19 beyond 305, 82% and 69%;
  - the 38 from 234 to 305, 14% and 25%;
  - the 193 up to 234, 4% and 6%.
- **What the estimate leaves out.**
  - Stage A and stage B ran under C. Under Cw, stage A's depth-2 nodes from
    234 to 421 short would also run a phase, so their depth-3 children would
    start closer than their cold index says. Leaving that out errs high.
  - The extra root phases, about 7.7e5 CPU s across the frontier (0.4%), are
    in the anchors' level 3.
  - The draws differ between arms (each probe's seed has its name), so a pair
    compares two estimates of the same subtree, not the same nodes.
- **The laptop.** At about 2x the cloud box per DNN iteration, the cold
  estimate is 19,000-21,000 laptop core-hours: 3,100-3,500 hours, or 129-146
  days, on its 6 P-cores. Under C it was 210-230 days.
- **What to use.** Cw. Its mechanism shows in nearly every pair, its extra
  phases cost 0.4%, and without its best pair it still comes out at
  0.81-0.95x of C.

## Replicates on the 380.1-short node (from Nov 13 02:49)

Without the 380.1-short pair, Cw / C is 0.81-0.95, not 0.62-0.63. Both of
its probes are noisy: C's SE is 79% of its value and Cw's 58%. So the node
gets two more probes in each arm, with `--seed 2` (each probe's draws are
seeded by the seed, its name and the level):
- C: Cr368b and Cr368c on cpus 0 and 1 (`C.sh`, `C.par`);
- Cw: Zr368b and Zr368c on cpus 2 and 3 (`Cw.sh`, `Cw.par`).

All four use the run binary (md5 `ee4c6ebd…`) and k = 16, beta = 0.025, as
before.

**How they enter (decided before any finished).** Each arm's value at the
node becomes the mean of its three probes, each an unbiased estimate of the
same subtree's cost. `estimate_se.py` and `paired.py` read that mean as
`anchors_C/Cr368+anchors_C/Cr368b+anchors_C/Cr368c=380.1`. Everything else is
unchanged.
