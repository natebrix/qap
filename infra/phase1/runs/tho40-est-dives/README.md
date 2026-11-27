# tho40-est-dives: the dives run of the tho40-est estimate (Nov 18)

Gathered from the store on Nov 18 after both instances were terminated (docs/review_nov18_runs.md, item 1):
524 of the 552 dives units merged; the other 28 (all small lane, depth 5-6) were queued but never leased.

- `dives.jsonl.gz`: every merged unit's `units/UID/aATTEMPT/dives.jsonl`, concatenated in merge order
  (11,168 lines: `sample`, stratum and dive records, as coord/dives.py reads them)
- `units.tsv`: one row per unit, merged or not; `line_from`/`line_to` are its 0-based line range in
  dives.jsonl.gz, `has_dives` whether the store held its file
- `estimate.txt`, `estimate.json`: `qapctl estimate tho40-est-dives`, part-way (28 units pending)
- `run.json` (the run's ledger record, without its solver entry), `status.json`, `source_totals.json`,
  `params.par`, `incumbent.json`: as the store keeps them
