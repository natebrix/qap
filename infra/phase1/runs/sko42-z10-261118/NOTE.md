# sko42-z10-261118: the cut.zone 0.010 arm, launched without --trace

One of two arms on sko42 that differ in `cut.zone`: this one 0.010 (`experiments/estimates/sko42/sko42.par`),
`sko42-z25-261118` 0.025 (`sko42_zone25.par`). Same commit (bbf15b89), instance type (c7a.16xlarge, Spot,
eu-north-1c, 64 slots), incumbent (15814) and unit sizing, launched three minutes apart on Nov 18.

The unit-sizing line typed at start.sh's prompt was cut short, so this arm's COORD lacks
`--ckpt-upload 300 --trace` and the other arm's has them. Kept as is rather than relaunched. So:

- no node or DNN traces for this arm: compare total CPU and nodes (and the per-level table), not per-node
  mechanisms such as the children's first-solve fathom rate;
- the z25 arm pays the tracing cost (a few % of CPU on nug30, README.md), which favours this arm slightly;
- the checkpoint upload interval is the default here, 300 s there.
