# Continuation (arm C) runs, fourth attempt, on the Fable cloud box (Nov 3)

The four C runs of `../../PROTOCOL.md` (d1, d2, d4, s5) lost their executor
three times (`../../recovery_plan3.json`). At Nathan's request they were run
here with the same commands, the same source (this branch at 1f81133, built in
a separate worktree) and the committed parent pool `../runs/parent.pool`:
`../../run_recovery3_cloud.sh`, four single-threaded processes on four cores,
22:05-22:48 UTC, no interruption.

Host: 4 x Intel Xeon @ 2.1 GHz (AMX-capable, unused by this source),
scipy-openblas32 0.3.34 (the original host: Xeon Platinum 8573C, OpenBLAS
0.3.30, eight processes). Reproduction check first (`s11_C_repro.log`):
s11's C run again loads the same 12,358 kept / 642 dropped cuts and fathoms
at the same iteration 160, certificate 48,819.83 against 48,819.77 in
`../runs/s11_C.log`; the difference is the BLAS version. Certificates and
rounds are comparable with the original runs; CPU is not paired with them.

## Results

Iterations: the 400-iteration first solve with the inherited pool, twelve
300-iteration Polyak rounds, then the quiet round (at most 1,500). CPU is
whole-process user + system from the shell's timer (`PROCESS_CPU` line).

| child | start (round 0) | best Polyak round | quiet round | retained certificate | fathomed | iterations | process CPU | first reaching the target |
|---|---|---|---|---|---|---|---|---|
| d1 | 48,753.31 | 48,769.52 (round 12) | 48,814.26 at 790 | **48,814.26** | yes | 4,790 | 2,419 s | the cutoff 48,814.05: quiet round |
| d2 | 48,725.84 | 48,737.76 (round 12) | 48,797.50 at 1,500 | **48,797.50** | no | 5,500 | 2,805 s | the own-phase value 48,789.9: quiet round (no Polyak round reached it) |
| d4 | 48,748.44 | 48,758.97 (round 12) | 48,814.59 at 1,020 | **48,814.59** | yes | 5,020 | 2,542 s | the cutoff: quiet round |
| s5 | 48,779.92 | 48,800.48 (round 12) | 48,814.88 at 400 | **48,814.88** | yes | 4,400 | 2,281 s | the cutoff: quiet round (no own-phase target) |

With arm B of `../runs/frozen_checkpoint.json` (one frozen quiet solve of at
most 1,500 iterations on the inherited pool):

| child | B certificate | B fathomed | B iterations | C certificate | C fathomed | C iterations |
|---|---|---|---|---|---|---|
| d1 | 48,814.09 | yes | 1,490 | 48,814.26 | yes | 4,790 |
| d2 | 48,780.01 | no | 1,500 | 48,797.50 | no | 5,500 |
| d4 | 48,808.51 | no | 1,500 | 48,814.59 | yes | 5,020 |
| s5 | 48,814.35 | yes | 760 | 48,814.88 | yes | 4,400 |

## Observations (for the study's report; no gate depends on them)

- **The Polyak rounds add little on top of the inherited multipliers; the
  quiet round does the work.** The best Polyak round is only +10 to +21
  above each child's start, and every certificate that crosses a target does
  so in the quiet round, which restarts from the best round's multipliers.
- **The first multiplier steps overshoot.** Every child's bound drops below
  its start for several rounds after a few good ones: d1 by 1,399 at round 4
  (recovering to +16 by round 12), s5 by 181 at round 3, d2 by 43 at round 5,
  d4 by 23 at round 4. The step factor halves on each drop, so the rounds
  recover, but a third to a half of them are spent recovering. The Polyak
  target (best + 50, uncapped) is calibrated for a phase that starts from an
  uncut multiplier, not from a pool already near its optimum.
- **Against B,** C fathoms d4 (B stopped 5.5 short) and lifts d2 by 17.5 to
  48,797.50, above its own-phase value 48,789.9, at about 3.3 times B's
  iterations. d1 and s5 were already fathomed by B at 1,490 and 760
  iterations.
