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

Why probes: the record's two sko49 depth-3 subtrees, 263-304 short converged,
were estimated at 47,000 laptop CPU s each under the Nov 5 settings
(`../../../sko49_feasibility/vs_wil50.md`), about 25 core-hours on this box.
sko42's largest was 1.3. A probe bounds about 24 nodes a level and weights them,
so it costs about 100-150 nodes whatever the subtree's size.
