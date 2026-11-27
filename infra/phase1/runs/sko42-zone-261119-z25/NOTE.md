# sko42-zone-261119-z25: its worker went from on demand to Spot at about 04:30 UTC, Nov 19

One of the two arms of the `cut.zone` A/B (`infra/phase1/specs/sko42-zone.spec`; the other arm is
`sko42-zone-261119-z10`). Both were launched at 02:22 UTC with one on-demand `c7a.16xlarge` worker each
(`MARKET=on-demand`), next to their on-demand `c7a.2xlarge` coordinators (7 slots).

- **What happened:** the on-demand vCPU quota fit the four coordinators then running (these two and the
  two estimates'), not the workers. Every launch of either arm's worker group failed with
  `VcpuLimitExceeded`, from 02:31 to 04:03 (the group's activities, `collect.sh`). Until the switch both
  arms ran on their coordinators' 7 slots only. z10 at 04:27 had 117 of 265 units merged, 12.9 CPU h.
- **The switch:** at about 04:30 UTC both groups' mixed-instances policy was set to Spot by hand
  (`aws autoscaling update-auto-scaling-group ... OnDemandPercentageAboveBaseCapacity 0`,
  price-capacity-optimized over c7a, m7a and r7a.16xlarge, as `launch.sh --workers` makes a Spot group).
  Nothing else changed: same recipe, coordinator, slots and units.
- **When the workers come:** not before the Spot quota has room. The tho40-est2 and tai50b-est groups
  held all 256 Spot vCPUs until those estimates ended. Each group retries on its own, so the two arms'
  workers may start some minutes apart.
- **What it changes in the comparison:** nothing in the tree or in the CPU totals, which the A/B
  compares; the arms' wall times are not comparable (they were not before either). A Spot interruption
  would add requeued units and lost attempts to that arm's record; `collect.sh` shows them.
- **Why it could happen:** `start.sh` checked the worker count against the quota per instance and, with
  a spec, launched regardless. Since ff210e95 it adds up the whole launch's needs per quota and does not
  launch past one with `--yes`.
