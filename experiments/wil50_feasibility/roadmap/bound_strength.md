# wil50 in a week: the bound-strength angle (depths 3-5)

Memo, Nov 2 2026. Nothing was run.

**Update Nov 4:** a three-agent study (`new_cuts/SYNTHESIS.md`) closed item 3
(the regional M2 block: every single threshold is already saturated after the
triangle phase) and item 5 (mixed-region and rectangle triangles: 0-5), and
found no family that adds +40 on top of the triangles. Item 1 (inheritance)
passed its gate (Astra, `research/pool-inheritance`). Every number is either **[m]**, measured and
read from the named file (paths relative to `experiments/wil50_feasibility/stronger_bounds/`
unless they start with `latest/`), or **[e]**, my estimate. One level is +130
units at depths 3-5, worth 47x on the tree (`../NOTES.md` section 5). A gain of
d units on every node is worth 47^(d/130): +10 is 1.34x, +20 is 1.8x, +40 is 3.3x.
A gain that only acts where the cut phase runs is worth much less.

## Summary

The certified bound at depths 3-4 is close to the ceiling of the linear
regional families. The triangle pool is within about 10 of its dual optimum
[m, `dualtest/README.md`]. The 5-point inequalities [m, step 7], the triples
[m, `triples/README.md`] and the local-7 cuts [m, `local7/README.md`] add
nothing on top of it. After a 1,500-iteration quiet round the certificate
trails the Lagrangian primal by only **3-8** [m: `dual/tri_{10,13,d2}.log`
3.8/4.7/7.2, `dualtest/slow30_*.log` 5.3-7.9, `local7/*_{13,d2}.log` 3.0-5.9].
So neither a better dual method nor a better certificate is hiding another
level.

I do not expect another full level (+130) from new inequalities within a
month; I put it below 10% [e]. The one strength candidate with a plausible
+30 is the regional moment block, and it needs a gate test first.

The largest bound-side lever I can find is not a new inequality. It is
**inheriting the parent's cut pool and multipliers in the children.** Today
the children keep the uncut multiplier (`latest/bnb.C` around line 660). The
argument rests on one fact. Under a fixed multiplier mu, a child's exact
Lagrangian value is at least its parent's: the child's DNN feasible set maps
linearly (Y_p = T Y_c T') into the parent's DNN feasible set, with the same
objective. So each child of a phase that ended X short starts at most X short,
before its own one-level gain. In the node-10 subtree that would reach the
~317 children of the 8 failed phases, about 47% of the cut arm's CPU
(derived below). It would also reach the 47 depth-4 children of the
depth-3 clear miss itself.

## Ranked table

| # | idea | expected gain | cost per node | effort | risk | evidence | cheapest first test / kill result |
|---|---|---|---|---|---|---|---|
| 1 | **Inherit cut pool + multipliers into children**; also run a phase (no quiet round needed for fathoming) at depth-3 clear misses outside the zone, to equip their 47 children | **1.5-4x on the dominant class** on top of the measured 1.81x [e] | +2% per iteration for `L_add` [m, README step 1]; removes most child phases (~900 s each [m, derived]) | 2-4 days, including state-file v5 for the pool | medium: the transfer is bracketed but not measured | exact inequality above; node 10's cut value 48,743.9 [m, `dual/tri_10.log`] against uncut 48,655.9; children's uncut gains +96/+49/+92 (d1/d2/d4) [m, `depth4/uncut_*.log` via `research_bound_quality.md` 1.2] | dump node 10's pool; reload it frozen on d1, d2, d4 and 3 random siblings (`-r 0`, quiet 1,500); about 1 laptop-hour. **Kill** if median (inherited - child uncut converged) < +40, or d2 < its own phase (48,789.9) - 30 |
| 2 | **Skew-symmetric polishing of the certificate** (Hahn-Grant / RLT1 dual ascent on the residual L+Z') | flips nodes ≤8 short after the quiet round; mainly shorter quiet rounds: **1.05-1.15x** on the phase side [e]; more if applied to plain nodes (speed angle) | 10-30 GL evaluations per certificate, about 0.2-1.5 s [e] against about 4 s of iterations between certificates [e, memo 3.3] | 2-4 days | medium: a near-optimal residual may already be nearly symmetric-tight | lag at 300-iteration rounds 25-90 [m, memo 1.2]; quiet rounds are 20-54% of a fathoming phase's iterations [m, `laptop/README.md` section 2] | dump L+Z' at the end of round 6 of d6 and round 3 of d1; ascent offline in Python (`scipy linear_sum_assignment`). **Kill** if it recovers < 20% of primal minus certified |
| 3 | **Regional M2 block as Lagrangian decomposition** (= SOS-certified moment cuts, Astra 10b): DNN unchanged, one certified regional SDP per mid region per round | +10 to +40 at the phase [e]; P(≥ +65) about 10-15% [e]; net 0.8-1.5x [e] | 4 regions x 50-300 ADMM iterations x ~0.3 s eig (n = 1,129) per round: +40-200% on a phase [e] | 1 week certificate and validation, 2-3 weeks in all | high | local-7 gate violations 0.13-0.20 [m]; but 5-point wash, local-7 pool negative, triples tiny [m] | SCS projection of node 13's **triangle-capped** quiet iterate (col<5, row<3) onto the M2 set, alongside the projection onto MET. **Kill** if the excess normalized violation < 0.03, or the fixed-mu gain on top of the triangle pool's best multipliers is < +10 |
| 4 | Per-cut normalized Polyak step (unit norm in moment space), then the local-7 closing check (b) under it | +5 to +20 over triangles alone [e]; **1.05-1.2x** | none | 0.5-1 day | low effort, low value | local-7 alone reaches 84-87% of the triangles' gain; the joint pool lands 15-25 below triangles alone [m, `local7/README.md`] | slow30 joint pool, normalized, on 13 and d2. **Kill** at ≤ triangle slow30 + 5 (48,762 / 48,795) |
| 5 | Mixed-region and rectangle triangles (the only regional linear family not yet tried in a bound) | +5 to +15 [e]; 1.05-1.15x | separation +5-10 s per round [e] | 1-2 days | low value | same max violation 0.091 as the threshold family, near-duplicate top cuts [m, memo 1.3] | 12-round protocol on 13 and d2. **Kill** at < +5 |
| 6 | Cut-phase policy from round-2 data: abort, or extend the quiet round to 2,500 | abort: **≤ 3-5%** [e]; extend: flips the 0-10-short endings | none | 0.5 day plus per-phase logging | a false abort costs about 4-5k s | see Q5 | log per-phase lines in the next search run (no compute) |
| 7 | Cuts as an in-ADMM block (MADAM-style) | no value above the pool's optimum (≤ +10 [m, dualtest]); speed only | +5-10% per iteration [e] | 1-2 weeks | high: PRSM's entrywise Y-step becomes a QP; `fast60` lost [m] | `dualtest/README.md` | defer |
| 8 | Literature relaxations (Lagrangian DNN, RLT2 on subsets, vertex-fixing lifts, RLT3) | 0 against DNN plus regional cuts [e] | - | - | - | Q3 | none |

## Top three

**1. Inheritance.** Bracket for the test, using the node-10 children. Each
inherited value is at least node 10's exact L(mu), about 48,744-48,748
[m: certificate 48,743.9 and quiet residual 3.8, `dual/tri_10.log`]. It is at
most about each child's own pool optimum: d2 reaches 48,789.9 [m,
`dual/tri_d2.log`], and d1 and d4 fathom [m, `laptop/README.md`]. If the
children's uncut level gains (+96, +49, +92) carried over additively, d1
and d4 would start above the cutoff and d2 would start about 21 short, close
to its own phase result. That would mean the parent's phase buys each child
its level for the price of one plain solve.

Mechanics in the library. A cut on facilities other than the newly fixed
facility f keeps its meaning: for i ≠ f, x_ai = 0, so z_i^S is unchanged.
A cut touching f collapses to a constraint on the complement moment, which
`Y >= 0` already implies, so it can be dropped. Children should also
warm-start from the parent's final cut-phase Z, not the uncut one.

Cost structure today [m, derived from `search/wil50n10ext_cut.log`]: 76,755
CPU-s in all, 41,796 in `dnn_bound` (365 nodes, 115 s each), and about
34,959 outside it (cut phases plus QPB), about 900 s per phase. About 317
of the 365 nodes sit below the 8 failed phases (365 - 1 - 47 depth-4), about
36k s.

My end-to-end estimate for the node-10 subtree with inheritance, including a
phase at node 10 itself: one phase (~2,300 Xeon-s [e], twice the 1,159
laptop-s of `dual/tri_10.log`), plus 47 children at ~100-120 s, plus 15-30%
residual misses. That comes to 15-35k s against 76,755: 2-5x. I quote
1.5-4x to allow for multiplier mismatch [e]. BiqMac/BiqBin also pass the
parent's cuts to the children (from memory; check before citing).

**2. Certificate polishing.** For any skew-symmetric K, <K, Y> = 0 on
symmetric Y. So GL(L + Z' + K) - (m+1)·lambda bound is a valid certificate
for every K. The spectral term is unchanged, because V'KV is skew and
lambda_max is taken of the symmetric part. The polyhedral set becomes
P_GL ∩ {Y = Y'}, the Adams-Johnson RLT1 LP. Hahn-Grant ascent climbs it with
GL-cost passes: about 20-50 ms each at m = 47 [e, memo 3.3], 2,209 inner
LAPs.

At the dual optimum it adds nothing, since every P that contains the DNN set
has the same Lagrangian dual. Its value is entirely at non-converged Z',
which is exactly the 300-iteration lag of 25-90 and the 600-iteration
quiet-round residual of 25-48 [m: `dual/both_*.log`, `dual/quiet4_*.log`].
The same lag exists at the plain-DNN give-up: a cold 400-iteration bound
sits 48-55 below its converged value at m = 46 [m, README step 4]. If the
offline test recovers half of it, the speed agent should own the plain-node
use. This angle claims only the shorter quiet round.

**3. The regional M2 block, gated.** Dualize the links u_i = p_i and
u_ij = W_ij (Astra REPORT section 7). The DNN side then stays exactly as
today, `L_add = sum_S A_S^* lambda_S`. Each region contributes a certified
min over its moment set, b'y + tau_S·lambda_min(slack), with tau_S =
1 + k + C(k,2), which is 277 at k = 23. lambda_min is bounded by the same
Rump Cholesky the DNN certificate uses, so a 1e-6 error costs about 3e-4
units.

At fixed lambda this is one SOS-certified quadratic cut per region. It is
therefore the same object as Astra's 10b separation, and it fits the existing
pool and Polyak loop.

What it can add beyond triangles: degree-4 PSD consistency over all m
members of a region. It does **not** imply the 5-point inequalities (Laurent
2003) or the local-7 polytope. All three measured local families have now
come out flat or overlapping, so I hold the median at about +20 [e].

## The five questions, short answers

**Q1, heterogeneous families.** The joint pool contains the triangle pool,
so its optimum is at least the triangles'. The 15-25 loss is dynamics [m].
In order of effort:

- unit-norm cuts or a per-cut diagonal metric, which divides the step by
  ||G_t||²: 0.5 day;
- per-family step blocks: 1 day;
- a proximal bundle with block-diagonal prox: 3-5 days;
- an in-ADMM block: 1-2 weeks, with no convergence proof for our splitting.

Value if well handled: the local-7 excess over the triangles only.
Sevenonly reaches 84-87% of the triangles' gain, so the overlap is large.
Scaling the fixed-mu +12 [m] by the triangles' fixed-to-optimized 4x, then
by a 15-50% non-overlap, gives **+5 to +20, about +10** [e]. That is a
third of the hoped +40-50, and 1.05-1.2x [e].

**Q2, moment block.** See item 3. In-ADMM consensus would cost +0.15x eig
per region, +0.6x for four regions, about +36% per iteration [e: memo 3.5,
with eig at about 60% of an iteration from the explainer]. There is no
cheaper exact relaxation that captures more:

- a 128-weight local-polytope block in the Y-step equals the local-7 family
  at its dual optimum, which is item 4's ceiling;
- a 29x29 level-2 block per 7-subset is weaker than that polytope;
- regional SOC or [1 p'; p W] ⪰ 0 is a compression of Y and already implied.

Only global-per-region (M2) or cross-region (item 5) consistency is new.

**Q3, other relaxations.** At n = 50 nothing in 2020-2026 beats the DNN:

- QAPLIB 48,121 and Graham et al. 48,126 [m, as quoted in Astra's
  `REPORT.md` section 2], against our root 48,135;
- Fujii et al.'s Lagrangian DNN is the same relaxation, with a weaker
  bracketed bound;
- Mittelmann-Peng 46,526;
- RLT2 is below the DNN in every computed value at n ≥ 25: nug30 5,846
  against 5,944 [m, `../NOTES.md` section 5]. RLT2 on subsets is the
  triples [m, retired].
- de Klerk-Sotirov vertex-fixing lifts amount to our depth-1 branching; the
  grid's symmetry group has order 4;
- RLT3 needs 10.1 TB.

None is new against DNN plus regional cuts at depth 3. Those regional cuts
are already beyond any published n = 50 bound.

**Q4, certificate.** After 1,500 quiet iterations at most 3-8 remains
[m]. A certificate fix flips only nodes ending ≤8 short. Of the
1,500-quiet endings in the logs, none sits at 0-8 short: d2 is 24.2 and
18.3; `laptop/dual/short6_d1.log` is 16.4 from last-round multipliers [m].
The prize is time. Quiet rounds are 960-1,570 iterations of a fathoming
phase, and item 2 is the candidate. An end-of-phase LP over P_GL ∩ pool
(memo candidate 4) is bounded by the residual plus the multiplier slack,
≤ 8 + 10. That is minutes of LP for at most the 3-25-short tail: do it
only if item 2 fails.

**Q5, predicting failures.** Under the search protocol (300-iteration rounds,
delta 50, from the give-up point), the shortfall after round 2 separates
the 8 nodes with known outcomes [m, parsed from
`laptop/dual/{short*,clear_*}.log`, `dual/tri_{13,14,d2}.log`]:

| outcome | nodes | shortfall after round 2 | primal - cutoff at round 2 |
|---|---|---|---|
| fathomed | d1, d4, d6, d11, 14 | 67-78 | -27 to +1 |
| failed | d2, 10, 13 | 102-148 | ≤ -43 |

Two nodes are undecided: d7 at 74.5 and d14 at 81.6 (`depth4/cuts_{7,14}.log`,
no quiet round). So a rule of "abort if > 95 after round 2" fits the sample.
But the stall rule already stops failures at round 4-5 [m,
`laptop/dual/clear_*.log`], so aborting saves 2-3 rounds, 320-480 Xeon-s
[e]. Over 8 failures that is **3-5%** of 76,755 [e]. One false abort
forfeits a fathom worth 45 children, about 4-5k s [e]. The failures' real
cost is their ~36k s of children, which item 1 targets.

Better uses of the round-2 signal: give predicted ≤15-short nodes a
2,500-iteration quiet round (d1 was 0.009 short at 1,500 [m]), and add
per-phase lines to the search log. `search/wil50n10ext_cut.log` has no
per-phase data, so the 39-phase sample is lost.

## Interactions

- **Speed.** Phase economics scale with the per-iteration cost: a 2x faster
  eig halves both the phases and the plain nodes, so the ranking stays. Item
  2 is a plain-node lever, since the give-up lag is 48-55 at m = 46. Item 1
  makes every child a cut node: `L_add` is +2% per iteration [m] and a 36 MB
  dense matrix at m = 46 [e]. Carry the pool as tuples (about 100 KB), not
  the matrix. As `latest/state.C` v4 showed, anything carried by the warm
  start must also go into the state file.
- **Overlap with `node_speed.md` item 5.** The speed memo proposes the same
  inheritance, for the children of failed phases only, and prices it as a
  speed gain of 1.1-1.25x. The difference here is the bound bracket (each
  child starts at least at the parent's cut value) and the policy of running
  a phase at out-of-zone depth-3 clear misses to equip their children. That
  policy is where the 1.5-4x comes from. One test, the node-10 children with
  the parent's pool, settles both estimates.
- **Tree.** The 130 units per level and 47x per level apply to uniform
  gains. Items 3-5 act only at phase nodes, so their factors are lower.
  Item 1 changes which nodes are near misses at depth 4-5, and with it the
  zone. It also suggests a branching tie-break: fixing a facility drops only
  the cuts that touch it, so branching on a facility with little cut
  multiplier mass preserves more of the inherited dual [e]. A one-level
  gain reshapes the tree agent's per-level counts below depth 3, so the
  inheritance test should feed their model before either of us extrapolates.
- **Budget.** Items 1, 2 and 4 together are 1.6-5x [e], mostly from item 1.
  That covers part of the 4-12x needed. A +130 level from new inequalities
  is not in sight, and I would not plan the week on one.
