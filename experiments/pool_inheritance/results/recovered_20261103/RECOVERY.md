# Recovery checkpoint — 2026-11-03

Remote head remained 1f811332e51d7891f9e70ba7bfac4023a66d4fe6. The original workspace was accessible, revealing uncommitted final-retry progress. This isolated branch preserves a separate snapshot of its runs, analysis edits as an unapplied patch, and an audit. The original workspace was not modified.

The last recorded heartbeat is 2026-11-03 19:44:17 UTC. Snapshot capture began at 20:06 UTC. Process visibility is isolated; lack of visible processes does not establish that the original executor stopped. These files establish past activity, not current liveness.

Fifteen completed runs were reparsed and agree with their manifest. All six A/B pairs survive; their numerical gates pass (median gain 54.364110; d2 frozen bound 48780.010403). This is NOT a qualifying completed study: four C runs lack final results. The complete-study analyzer correctly rejects this manifest. Plan, instance and parent-pool hashes match.

Recovered C progress: d1, d4 and s5 reached round 11 (3700 iterations); d2 reached round 10 (3400 iterations). All numerical prefixes match the original interrupted attempt. No completed C results, terminal CPU, or target crossings have been invented. Pending final quiet solves could materially change the outcome.

No solvers were launched. PROTOCOL.md explicitly ends retries after this final attempt and retains the original six-hour deadline (21:04:53 UTC). The saved files contain logs and a parent pool, not resumable child Y/Z solver states. Restarting these C runs would be a new attempt, not continuation from round 11, and cannot be silently counted under the original protocol. Completion would require either retrieving the final retry's terminal output or explicitly documenting a separate follow-up execution with a new deadline and attempt accounting. The user's takeover authorization enabled this recovery; the original study remains honestly incomplete.

The unapplied analysis patch is preserved as evidence, not endorsed as validated code. audit.json reports only checks actually performed.
