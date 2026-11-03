# Fourth attempt: incomplete checkpoint

Checked November 3, 2026 at 13:40 Pacific. The persistent terminal session returned `Unknown process id 7892`. The last saved heartbeat was 13:14:20 Pacific, approximately 26 minutes earlier. The manifest still says running, but this stale label is not a liveness check. No solver failure or terminal result appears. Executor loss is suspected; current OS process liveness cannot be verified here.

Four C runs remain incomplete: d1/d2/d4 reached round 8 (2,800 cumulative iterations); s5 reached round 6 (2,200). Every observed numerical prefix exactly matches the earlier saved run, excluding timings. Binary, parent-pool, instance and original-plan hashes match. The 15 completed runs remain byte-identical to the prior snapshot and reparse consistently with the manifest. Their frozen-pool gates remain unchanged; the full study is not qualified.

No further attempts were launched. No terminal CPU measurements are available for these four partial attempts. Logged DNN call CPU sums are only lower bounds on work spent (see interruption_audit.json); elapsed wall time is not substituted for CPU. Historical partial attempts and completed results are retained separately.

The entire output directory is preserved here. This checkpoint does not imply completion or silently turn partial rounds into final certificates. A reliable execution environment or resumable solver state is needed before another restart is sensible.
