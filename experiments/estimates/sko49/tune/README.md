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

Why probes: the record's two sko49 depth-3 subtrees, 263-304 short converged,
were estimated at 47,000 laptop CPU s each under the Nov 5 settings
(`../../../sko49_feasibility/vs_wil50.md`), about 25 core-hours on this box.
sko42's largest was 1.3. A probe bounds about 24 nodes a level and weights them,
so it costs about 100-150 nodes whatever the subtree's size.
