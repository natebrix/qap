- **The root's plateau** (section 3): 6,470 iterations, still climbing, no curve.
- **A separate random sample of 40 or more depth-2 nodes, and depth 3 below them**: the dives' start nodes are
  that sample instead (the table of section 6 counts them by depth; `logs/dives_SET.txt` has each one).
- **Thirty complete dives through the web.** The dives number well over 200, but most start and end at
  one fathoming node; the ones that matter (a branch at depths 4 to 6, section 6) are a handful, and levels
  4, 5, 6 and 8 have few dives. E1 and E2's intervals understate that.
- **Any tuning A/B**: one settings arm (`par/A.par`) throughout. The phase limit, the stall rule, the region
  set, `pool.cp`, `dnn.amx` and the root's length are all open (section 9).
- **The 32 unsampled depth-1 nodes**, among them `31:0` (a twin by the landscape) and `31:42`.
- **The spine below depth 13 and any twin below depth 4**; the path check stops there too.
- **End-to-end validity on solved instances other than tai30b** (tai35b, tai40b with their own region
  files), which the regions session listed as the next check and which should come before a tai50b run.
- **The run itself.** At this size the next step is not a better estimate but the search: `bnbdrv` from
  the saved long root with `par/A.par` (phase limit raised), under the coordinator or `bench/parallel.py`.
