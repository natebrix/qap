# Params files, version 2: named blocks, bounds by reference, per-block DNN settings (design, Nov 25 2026)

Status: design, not implemented. To be built after wil50's run and, if its estimate is reasonable, sko49's
(their recipes passed their path checks in the current format). Nathan's proposal (Nov 24; his sketch is the
"Origin" section) with the changes agreed in that discussion. Decided: when it lands, every params file moves to
version 2 and the version-1 reader is removed; the last version-1 commit is tagged, so old runs and studies
reproduce from git.

## Why

1. **Positional blocks.** `latest/util.C` `qap_read_params` reads each block as 2 + 7 numbers, one per line, the
   rest of the line a comment. A missing line shifts every value after it (a validator now catches the worst
   cases, `util.C` after the read loop). The comments drift: `experiments/estimates/*/` files still say "all
   nodes thru depth 6 use rule 8" over a rule-9 block.
2. **Block selection has surprises.** `pick_param` (`latest/bnb.C:670`) takes the first block whose gap threshold
   the node's gap exceeds, then walks forward to the first block whose depth limit covers the node, without
   checking the gap again; the root is forced to block 0. Seen this week:
   - `infra/phase1/recipes/sko42_d2r10_sb.par` meant rule 10 at depth 2 only, but its depth-1 nodes with a gap
     at most 0.20 also fall into the rule-10 blocks; harmless only because strong branching overrides the line
     at depth 1.
   - `sko49_d2r10_sb.par` and `wil50_d2r10.par` needed 8 blocks to change depth 2 alone and keep depths 3-4 as
     they were (the check: an emulation of `pick_param` over depths 0-11 and gaps 0.001-1.5).
3. **"Where in the tree" lives in many places, each with its own syntax.** The blocks vary the branching rule
   and QPB's budgets by gap and depth; everything else varies by depth through keys: `dnn.depths` (D1:D2),
   `dnn.ext` (levels below D2), `dnn.maxgap`, `cut.maxdepth`, `cut.giveup_depth`, `sb.depths` (D1:D2), `sb.all`
   (a depth), `sb.power` (D:P,...), `sb.gap` (D:G,...).
4. **Three different gaps.** Blocks use (inc - max(inherited bound, predicted bound)) / (inc - root QPB), known
   before the node is bounded; `dnn.maxgap` uses the node's own QPB bound; `sb.gap` its DNN bound. All have the
   same denominator, the root's QPB bound (`stats->root_bound`, set at depth 0 from QPB, `bnb.C:1060`, also
   when the root gets a DNN bound).

## Origin

Nathan's sketch (Nov 24): a block names its branching rule, its primary bound and its fathoming bound; the QPB
parameters live in named sections (`[QPB-A]`, `[QPB-B]`, `[QPB-C]`) that blocks reference, as `[settings]`
holds the rest; a `[DNN]` section can override `[settings]`. Two questions in it are answered below: "-" or
NONE for no fathoming bound (the `bounds` list), and where `cut.giveup_depth` goes (it disappears).

Three changes from the sketch, agreed in the discussion:
- **Explicit block conditions** in place of the gap-then-depth walk (item 2 above).
- **An ordered list of bounds** in place of primary and fathoming. The DNN bound does more than fathom: a node's
  bound is max(QPB, DNN), its children inherit it and it sets their gap; rules 9 and 11 take candidate lines from
  the DNN's marginals; children warm-start from the parent's DNN state; strong branching needs it.
- **Run-wide settings kept out of blocks.** Some keys cannot differ between a parent and its children (below).

## The format

A version-2 file has no positional lines. It is a list of sections, `[type NAME]`, with `key = value` lines and
`#` comments. Names are unique across the file. The reader tells the versions apart by the first non-comment
line (`[` for version 2) and, after the migration, refuses version 1 with a message that names the tag.

| section | holds | referenced by |
|---|---|---|
| `[settings]` | run-wide keys, and the default of every per-node key | |
| `[block NAME]` | `when`, `branch`, `strong`, `bounds` | the search, in file order |
| `[branch NAME]` | the line rule and its parameters | a block's `branch` |
| `[sb NAME]` | DNN strong branching's parameters | a block's `strong` (optional) |
| `[qpb NAME]` | QPB's parameters | a block's `bounds` |
| `[dnn NAME]` | per-node DNN and cut-phase keys that differ from `[settings]` | a block's `bounds` |

`dnn` with no section is a built-in name: the DNN bound with `[settings]`' values.

### Blocks and how a node finds its block

```
[block d2-far]
when   = depth 2, gap > 0.12
branch = r10
bounds = qpb150, dnn
```

- `when = depth A` or `depth A:B` (inclusive; `A:` for A and deeper), and optionally `gap > G` and/or
  `gap <= H`. A block without `when` matches every node; exactly one such block, and it is last.
- A node takes the **first** block, in file order, whose condition it meets. The gap is the inherited one (as
  today; the block is chosen before the node is bounded).
- No block is forced: the root is depth 0 with gap 1. A translated file starts with `when = depth 0`.
- The reader reports a block that no node can reach (wholly shadowed by earlier ones) as an error.
- `--print-settings` prints every block resolved (its rule, QPB, DNN and cut values) after the run-wide keys.

### Bounds

`bounds = QPBNAME[, DNNNAME]`: tried in order until one fathoms. QPB comes first and is required, since the
rules use its U (rules 2, 4, 8) or fall back to it (rules 9, 11 where a node has no DNN state). A block without a
DNN name runs QPB alone: no `-` or NONE token. Other bounds (GLB first, say) are out of scope; the list leaves
room for them.

DNN on a node then depends on its block and on `dnn.maxgap` (a `[dnn]` key, the node's own QPB gap, as today).
`dnn.ext` stays run-wide with its meaning (a child of a DNN miss keeps the bound for `dnn.ext` more levels, the
count carried in the warm state, `bnb.C:1084-1090`, `1500`), anchored at the deepest depth at which any block lists
DNN (today: `dnn.depths`' D2).

### Keys

Version-1 fields and keys to version 2:

| version 1 | version 2 | section | note |
|---|---|---|---|
| block field 1, 2 (gap, depth) | `when` | block | conditions, not a walk |
| rule | `rule` | branch | |
| Frank-Wolfe iterations, for bound | `qpb.fw` | qpb | Nathan's name |
| Frank-Wolfe iterations, when we can't fathom | `qpb.fw_giveup` | qpb | Nathan's name |
| which bound (0-7) | `qpb.kind` | qpb | qp, qp-param, qp-glb, qp-imp, glb, pb, evb3, buff (`util.h:60`); every recipe uses 1, qp-param |
| step for parametric bound | `qpb.param_step` | qpb | |
| Frank-Wolfe iterations, for row selection | `lookahead_fw` | branch | Nathan's `qpb.select`: the QPB solves of candidate lines' children (`branch.C:1352-1760`), a branching parameter; they use the block's `qpb.kind` |
| compute B for best _ rows of U | `candidates` | branch | Nathan's `qpb.strong_best`; rule 9 reads it only where it falls back to rule 4; `all` for rule 10's 99 |
| `branch.la_exp`, `branch.xbar_lines`, `branch.weighted`, `branch.weight` | the same keys | branch | per block |
| `sb.depths`, `sb.all` | a block's `strong`; `sb.all = 1` | block; sb | |
| `sb.power` (D:P,...) | `sb.power = P` | sb | one value; the depth split moves to blocks |
| `sb.gap` (D:G,...) | `sb.gap = G` | sb | one value |
| `sb.k`, `sb.iters` | the same keys | sb | |
| `dnn.depths` | which blocks list DNN | block | removed as a key |
| `dnn.maxgap`, `dnn.giveup`, `dnn.giveup_warm`, `dnn.gu_zone`, `dnn.beta_sched`, `dnn.retry`, ... | the same keys | settings or dnn | per-node (below) |
| `cut.maxdepth` | a `[dnn]` section with `cut.zone = 0` | dnn | 0 is off (`bnb.C:1251`) |
| `cut.giveup_depth` | a `[dnn]` section with `cut.giveup` set, others without | dnn | removed as a key |
| `cut.zone`, `cut.giveup`, `cut.handover`, `cut.rounds`, `cut.iters`, `cut.cp`, ... | the same keys | settings or dnn | per-node |
| header: symmetry flag, depth limit | none; `-s` and `-d` | | every params file has the flag 0; 2 set a depth |

### Run-wide and per-node keys

The settings registry (`latest/settings.C`) gains a scope per key. A run-wide key in a `[dnn]`, `[qpb]`,
`[branch]` or `[sb]` section is an error. The rule: a key is run-wide if a node's value reaches another node's
computation beyond the bound it hands down: state formats, the pool's structure, the cut geometry, the warm chain.

- Run-wide: `ship.bits`, `dnn.warm`, `dnn.ext`, `pool.on`, `pool.cp`, `dnn.amx`, `cut.grid`, `cut.side`,
  `cut.regions`, `dnn.warm_y`, the `est.*` keys, `par.file`.
- Per-node: the give-up, zone, phase, penalty and retry keys above, `dnn.cheb`, `dnn.refresh`, `dnn.cross`,
  `dnn.lazycert`, `dnn.beta`, `dnn.gamma*`, `dnn.tight_lam`.
- To settle key by key in the implementation: `cut.zone > 0` also decides whether a node keeps its DNN state for
  its children (`bnb.C:1135`). With `cut.zone` per block, the state must still be kept wherever a child's block
  lists DNN.

### Hash and identity

Today the policy hash includes `par.hash`, a hash of the params file's bytes (`settings.C:108`,
`driver.C:592,628`), so a comment edit changes it. Version 2 hashes the canonical resolved form: the run-wide
policy keys, then each block in order with its condition and resolved branch, sb, qpb and dnn values.
`par.hash` stays in the dump for information and leaves the hash. A file and its faithful translation then
hash alike, and comments are free.

- State files (QAP-STATE 10) record the settings hash (`state.C:321`); `-r` refuses a different one without
  `--allow-settings-change`. A version-2 build does not resume a version-1 state: it says so and names the tag.
- The coordinator's root cache keys on the solver binary and `--print-settings` (`coord/coordinator.py`,
  around line 530): any new build misses it, so the migration costs one root per problem, as any rebuild does.
- A run in flight is not touched: its workers build at the run's launch commit (`launches.tsv`).

## Migration

1. **Tag** the last version-1 commit `params-v1`.
2. **Reader** in `util.C`/`settings.C`: sections, scopes, the checks above; the resolved blocks in a structure
   the search uses. `pick_param` becomes the first-match rule. The code that read settings once into
   `static const` locals (about 40: `bnb.C` 20, `estimate.C` 8, `branch.C` 7, `dnn.C` 4, `state.C` 1) reads a
   per-node key from the node's block instead. This is most of the work.
3. **Translator** `tools/params_v1_to_v2.py` (as `bench/envpar.py` was for the environment): version 1 in,
   version 2 out, behavior-preserving. It enumerates the old rule over depth and gap and emits one block per
   distinct (depth range, gap range, rule, QPB values, DNN on or off, sb values) in first-match order, merging
   neighbors that agree. It is the definition of equivalence.
4. **Equivalence checks**, all to pass before converting:
   - Every tracked params file (about 274; the other `.par` files are 35 settings-only `dnntest` files and about
     128 one-line estimator files): the old rule (a Python port of `qap_read_params` and `pick_param`, checked
     against the `params-v1` build on a sample) against the new reader on its translation, at depths 0 to n and a
     grid of gaps, for the resolved rule, QPB values, DNN on or off and every per-node key.
   - The same node trace (fixes and outcome per node) from the old build and the new one, on nug20 and nug22
     under `infra/phase1/recipes/nug30.par`; `tho40c_sb.par` and `sko42_d2r10_sb.par` to depth 2; `tai50b.par`
     (regions) to depth 2; `wil50_d2r10.par` to depth 1; one GLB and one PB file from `latest/PARAM/`.
   - `make test`: the rule-10 brute force, `sb_test`, ship, parallel, coordinator, node trace.
5. **Convert** every params file in one mechanical commit, run by the translator, with no hand edits; hand
   cleanups (merging blocks that only differ in a field the rule ignores) come after, each its own commit.
6. **Callers**:
   - `-x` and `-w` (overrides of `dnn.depths` and friends): removed. Live users: `tests/node_trace_test.py:361`
     and the `--extra` examples in `coord/qapctl.py` (convert to params files). About 60 scripts under
     `experiments/` pass `-x`; they are records, reproducible from the tag.
   - `estimate.C` and `console.C` call `pick_param` (`estimate.C:207,442`, `console.C:897`): same rule.
   - brixius-web `lib/qap/params.ts` (the strategies table on run pages) and `lib/qap/log.ts` (the solver's
     printed blocks) read version 2 as well. They keep the version-1 parser: published bundles of past runs carry
     version-1 `params.par` and do not change.
   - `RUNNING.md` (the format, and the generated settings table with its scope column), `--settings-map`,
     the `experiments/estimates/*/settings.golden` dumps (regenerated in the conversion commit).

Estimate [e]: 3-5 days, most of it the per-node settings plumbing and the equivalence runs.

## Tests the implementation adds

- Parser errors, one test each: a key in the wrong scope; an unknown key; a missing or misplaced catch-all; an
  unreachable block; an unknown or wrongly typed reference; a duplicate name; rule 10 or 11 without
  `branch.la_exp`; `strong` on a block without DNN; a version-1 file after the migration (the message names the
  tag).
- The hash: a comment edit leaves it unchanged; a translated file hashes as the resolved old one.
- Per-node keys take effect: a two-block nug20 file whose blocks differ only in `cut.giveup` gives the phase
  lengths of each block in the DNN trace.

## Worked example: `sko49_d2r10_sb.par` in version 2

Its eight version-1 blocks, as the old rule resolves them (the Nov 23 emulation): depth 0, rule 9 (lookahead 25,
candidates 30, QPB 150/100); depth 1 by gap above 0.20, 0.12-0.20, at most 0.12: rule 9 with (25, 30, 150/100),
(12, 3, 150/100), (25, 5, 100/50); depth 2 above or at most 0.12: rule 10 with QPB 150/100 or 100/50; depths 3-4
above 0.12: (12, 3, 150/100); everything else (25, 5, 100/50). Strong branching at depths 0 and 1, DNN at every
depth.

```
[settings]
dnn.warm = 1
ship.bits = 32
pool.on = 1
cut.grid = 7:7
cut.side = 1
cut.zone = 0.028
dnn.giveup = 300
dnn.giveup_warm = 100
cut.handover = 27:0:5
dnn.retry = 1
dnn.beta_sched = two,f=1.6667,t=60,w=0.25,wc=1.4
# ... the rest of sko49_d2r10_sb.par's [settings], less dnn.depths, sb.depths, sb.all and branch.la_exp

[qpb qpb150]
qpb.fw = 150
qpb.fw_giveup = 100
qpb.kind = qp-param
qpb.param_step = 30

[qpb qpb100]
qpb.fw = 100
qpb.fw_giveup = 50
qpb.kind = qp-param
qpb.param_step = 30

[branch r9-wide]
rule = 9
lookahead_fw = 25
candidates = 30

[branch r9-narrow]
rule = 9
lookahead_fw = 12
candidates = 3

[branch r9]
rule = 9
lookahead_fw = 25
candidates = 5

[branch r10]
rule = 10
lookahead_fw = 25
candidates = all
branch.la_exp = 32

[sb sb-root]
sb.all = 1
sb.power = 8

[sb sb-d1]
sb.power = 1

[block root]
when = depth 0
branch = r9-wide
strong = sb-root
bounds = qpb150, dnn

[block d1-far]
when = depth 1, gap > 0.20
branch = r9-wide
strong = sb-d1
bounds = qpb150, dnn

[block d1-mid]
when = depth 1, gap > 0.12
branch = r9-narrow
strong = sb-d1
bounds = qpb150, dnn

[block d1]
when = depth 1
branch = r9
strong = sb-d1
bounds = qpb100, dnn

[block d2-far]
when = depth 2, gap > 0.12
branch = r10
bounds = qpb150, dnn

[block d2]
when = depth 2
branch = r10
bounds = qpb100, dnn

[block d3-4-far]
when = depth 3:4, gap > 0.12
branch = r9-narrow
bounds = qpb150, dnn

[block rest]
branch = r9
bounds = qpb100, dnn
```

Still eight blocks, because the translation keeps every case the old file had; but each says what it is, and a
new recipe writes only the cases it means. The depth-1 variants differ only in fields strong branching overrides
(the line) or barely uses, a cleanup for after the conversion. `cut.giveup_depth`'s job, in a recipe that has it
(`experiments/nug30_record/par/g5e_r11x2_gu83_rt.par`, the laptop record: `cut.giveup = 100:8,300:3` from depth 2;
`infra/phase1/recipes/tai60b.par`: `1000:0` from depth 1), becomes a `[dnn deep]` section with `cut.giveup` that
the blocks from that depth list, and plain `dnn`, without it, above.

## Open questions for the implementation

- Whether `strong` on a block whose nodes may lack a DNN state is an error or a no-op (today the run stops at
  startup when `dnn.depths` does not cover `sb.depths`).
- Whether blocks may also gate on the node's own QPB gap after bounding (today only `dnn.maxgap` does), or that
  stays a `[dnn]` key.
- `qpb.kind = pb` sets the Frank-Wolfe budgets to 1 today (`util.C`, after the reader's loop): keep, or make
  `qpb.fw` required and checked.
