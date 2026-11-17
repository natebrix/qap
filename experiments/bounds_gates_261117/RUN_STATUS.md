# Final status

All experiments and verification jobs have completed. No production run or setting was changed. The completed evidence does not support a new whole-instance solve-time estimate.

- Coarse-label cuts: `26260591`, `research/coarse-label-cuts-261117`.
- CP scaling: `79c5dc37`, `research/cp-congruence-261117`.
- Regional equalities: `57e6a285`, `research/cp-region-equalities-261117`.
- Shared evidence and independent checks: `research/bounds-gates-261117`; see README and RESULTS_LEDGER.

The full baseline suite and the three implementations' correctness checks passed, with the optional AWS mock skipped for missing dependencies. The complete equality subtree evaluated identical ordered nodes and used more total DNN iterations. All tested coarse pools lost at equal CPU, including the longer 600-step tai60b check. Scaling showed no earlier endpoint-target crossing in its three screens.

Native cold tai60b stock screening was interrupted after 905 seconds elapsed without a completed node trace. No candidate native run followed it. Bounded tai60b comparisons instead used the same explicit 100-ordinary/200-CP seed and clearly labeled continuation budgets.

The seed patch is only for reconstructing comparison inputs. It is not in candidate solver binaries. Large generated warm blobs remain reproducible via REPLAY.md. Publication of the four local research branches is pending explicit approval after an automatic approval rejection; it is the only unfinished external action.
