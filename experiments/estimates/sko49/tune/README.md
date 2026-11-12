# sko49: an estimate with sko42's best settings (cloud box, from Nov 12)

sko49 gets the settings that won on sko42 (`../../sko42/tune/README.md`), on its own
7 x 7 grid, and the same sampling plan as sko42. Below depth 3, sampled probes
(`../../probe.py`) stand in for the complete subtrees, which are too large here.

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
- **Everything else** is `sko49.par` as it stands. `C.golden` differs from
  `../settings.golden` only in the zone, the hashes, and the pool-check
  diagnostic, which `../settings.sh` turns on and sko42's arms leave off.

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
