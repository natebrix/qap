# Certificate study: running notes (Nov 4, from 20:45 UTC)

Crash-recovery log; the findings are in REPORT.md.

## 20:45-21:10 reading, code, smoke tests (box busy until ~22:00)

- Worktree branch was at master 16a72e8; reset to cognac 99097cc (clean).
- Certificate anatomy (latest/dnn.C certify): lb = s + M00 + GL(M) - (m+1) scale
  lambda - margin, M = L + scale Zc, Zc = Z + beta (Y - VRV') (PRSM half-step
  multiplier; V'ZcV = beta * negative part of W exactly). GL inner LAPs: row-col
  reduction, exact Hungarian only within 1e-3 |cutoff| (48.8 units) of the cutoff.
- Margins in the wil50 logs: 3.6e-6 to 4.95e-6 units, refine=0, cert_fail=0
  (amx/logs/val_d4n3_amx1.log, warm_y/runs/child_9.log). The PSD term and the
  rounding margins have nothing to give.
- bnb.C never sets DNNParams::child_bounds: child_lb is unused (confirmed).
- warm_y/runs/child_9.log (uncut, warm child, pass z): primal ~48,813 from
  iteration 150 on, certificate 48,729.7 (150), 48,780.1 (300), 48,804.3 (600),
  48,810.4 (900). The lag is the dual converging slowly, not the primal.
- Traces: every cut-phase node of pool_search / coord_search had its uncut
  primal BELOW the cutoff at give-up (-1.9 to -112): no certificate change on
  the uncut solve could have spared any of the 13 phases.
- coord_search CPU (18,496 s): fathoming first solves ~9,200 s (cold d4/d5
  4,100; pooled 4,770), first-solve misses ~2,700, cut phases ~7,000 (from
  pool_search's split).
- Pooled children first-solve iterations (pool_search): 6 at 10, 47 at 20, 40
  at 30, ..., tail to 500. A free child bound can save at most the first check.

## Instrumentation (latest/dnn.C, DNN_CERTX, opt-in, default path unchanged)

Variants at the same iterate: 1 exact GL, 2 Hahn-Grant skew polishing, 4 line
search Zc + t (Z - Zc), 8 average of last 4 Zc, 16 dump M. dnntest -K takes a
';' list. Smoke: nug12 to convergence all variants agree with the DNN value
(567.990); instrumented and plain runs give bit-identical RESULT lines;
brute-force n = 8 (scripts/validity_bruteforce.py): no variant above the
enumerated optimum on 12 (sub)problems. child_gl.py's exact parent GL equals
dnn.C's to 4 decimals (nug12 dump).

nug20 node 0:1, 100 iterations: HG +4 to +40 over std, exact +0.6 to +7,
line search never above t = 0 (concave, decreasing), average below std while
the certificate rises.

## 21:00-21:20 more variants, validity

- Bit 32 (HGR): HG with row-column reductions in the inner blocks (O(m^4) a
  pass) and one exact pass at the end. nug20 node: 12+1 passes reach the 8-exact
  -pass value (2,524.91 vs 2,524.82 at it 100) at about half the time at m = 19.
- Bit 64 (warm skew K carried across checks): WORSE than plain exact GL
  (2,496.9 vs 2,521.2 at it 100, nug20): a K tuned to one check's M is stale at
  the next. The invariance claim itself holds: at the same check
  GL_exact(M + K) >= polish value (+0.17 to +7.1, DNN_CERTX_WTEST). Dead.
- validity_qaplib.py (all variants, every check, 1,500 its): nug12, had12,
  chr12a, tai12a, tai12b: no value above the optimum. tai12b is a sharp case:
  HG lifts the converged certificate 39,464,306 -> 39,464,921 (optimum
  39,464,925). logs/valid_*.log.
