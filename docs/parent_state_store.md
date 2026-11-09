# A parent-state store for the coordinator (review item C2)

Design note, Nov 8 2026. **Implemented Nov 9** as designed, with these
choices: the setting is the run-class environment variable `QAP_DNN_WARM_DIR`
(the registry, `RUNNING.md` "Settings"); the file name is the 64-bit FNV-1a of
the token (dedup only; a collision costs a warm start, never a bound); a
reference is relative to the state file's directory when the store is under
it, absolute otherwise; `parallel.py` resolves references on read, writes
task files with them, links or copies the referenced files into a
checkpoint's own `warm/` directory, deletes from the store at every
checkpoint the files no open node references and no running worker can still
be writing (older than the oldest running task's start), and has
`--checkpoint-inline`. `make shipcheck` (`ship64store`) checks a stop and
resume through the store is the same search node for node. The `coord/`
coordinator does not set the variable and so keeps self-contained unit files;
giving it a store is a separate step (its backends are S3 and local files).

## The problem

Since Nov 8 the solver ships each branched node's final DNN state to its
children (`QAP_DNN_SHIP=32` is the default; `state.C`, token
`dnn32:<m>:<scale>:<hasY>:<base64>`). At m = 40 the token is the upper
triangles of Z and Y, 2 x 1601 x 1602 / 2 fp32 values, 10.3 MB, 13.7 MB as
base64 (at 64 bits, the exact test setting, the whole matrices: Z and Y are
symmetric to rounding only since `DNN_CORE_NOSYMV`). The warm line carries it every time the parent is written:

- `bench/parallel.py` writes a task file for every launch. A task's nodes
  usually share a few parents, but with 349 tasks on nug30 and thousands on
  tho40 the same parent is re-emitted many times (every task is `-r`/`-o`
  on its own file, which the worker rewrites on exit or on a steal).
- The coordinator keeps every open node's warm text in memory (the `Warm`
  objects), and `State.write` builds the whole file as a list of strings,
  then joins it, then writes: about three copies of the checkpoint in
  memory. The tho40 depth-2 checkpoint was 3.0 GB on disk and the
  coordinator peaked at 9 GB.
- Dedup is by Python object identity. When a parent's nodes are split over
  two tasks, each worker rewrites that parent into its output file, the
  coordinator parses two `Warm` objects with the same text, and the next
  checkpoint holds the token twice.
- Every parse of a task file (each exit, each steal) reads and splits
  megabytes of base64 that the coordinator never looks at.

On 1,000 cores with wil50 (m up to 48 at depth 2: n1 = 2305, a token of
21 MB fp32, 28 MB base64) the task files alone would be tens of gigabytes
per hour of churn.

## The design: one file per parent, referenced by content hash

The DNN token moves out of the state file into a store directory; the warm
line keeps a reference. State version 9 (8 is the settings record,
`docs/settings_plan.md` Phase 1):

```
warm <m> <p[0..n-1]> <X...> <dnn_ext> <pool> @warm/3f2a...9c.dnn
```

- `@<path>` in the `<dnn>` field means: the token is the contents of that
  file, path relative to the directory of the state file being read. The
  file holds exactly today's token (`dnn32:...`), so the reader is the
  present `dnn_from_token` after one `ifstream` read; the `-` and inline
  forms stay valid, so every version 6 and 7 file still loads.
- The file name is the SHA-1 of the token. Two workers that produce the
  same parent (a stolen task redone after a resume) write the same bytes
  to the same name; `rename` of a temp file makes the write atomic and the
  second writer a no-op.
- X stays inline. At m = 40 it is 1,600 doubles, 32 KB of text, and the
  coordinator's parser needs nothing from the DNN field, so the warm line
  stays opaque to `parallel.py` (it already treats the token as text).

### Solver side (`state.C`, `assign.h`)

- `QAP_DNN_WARM_DIR=<dir>`: when set, `write_state` writes each parent's
  token to `<dir>/<sha1>.dnn` (if that file does not exist) and emits
  `@<dir relative to the state file>/<sha1>.dnn`. When unset, the token is
  inline as today, so a plain `bnbdrv -o F` run still makes a self-contained
  file.
- `WarmStart` gets a `std::string dnn_ref`: a parent read by reference keeps
  its reference and is re-emitted by reference without re-encoding (today
  a worker re-encodes every parent it writes back, 10 MB of base64 each).
  A parent created in this process has no ref until its first write, when
  the hash is computed once and stored.
- Hashing 10 MB is about 10 ms; the DNN bound that produced the state took
  minutes. SHA-1 from CommonCrypto on macOS and OpenSSL or a 150-line
  implementation elsewhere; a 64-bit FNV of the bytes would do as well for
  dedup, since a collision only costs a re-encoded token, never a wrong
  bound (every bound is certified independently of where the ADMM starts).

### Coordinator side (`bench/parallel.py`)

- `--warm-dir`, default `<workdir>/warm`, exported to the workers as
  `QAP_DNN_WARM_DIR`. Task files and the frontier live in `workdir`, so
  references are `warm/<hash>.dnn`.
- The checkpoint is written next to its own `warm/` directory: the
  coordinator writes `<checkpoint>` with references into
  `<dirname(checkpoint)>/warm/`, and hard-links (or copies, across file
  systems) the referenced files there from the work store. A checkpoint
  plus its `warm/` directory is the unit to keep or move; `RUNNING.md`
  says so. For a self-contained file there is `--checkpoint-inline`, which
  does today's work (reads every referenced file and inlines the tokens)
  for the one-off case of handing a state to another machine.
- Garbage collection at every checkpoint: the files referenced by the
  pool and by the in-flight tasks' nodes are the live set; everything else
  in the store is deleted. A file no open node references cannot be needed
  again from that checkpoint on, and a running worker only references the
  parents of its own nodes (in the live set) or files it is itself
  creating.
- The coordinator never reads a token. Memory per open node is the node
  line plus a shared 32 KB warm line, so a 300,000-node tho40 frontier is
  about 150 MB instead of several GB, and checkpoints take seconds.

### What stays the same

- Node lines, statistics, the pool token, `-r`/`-o`/`-t`/`-k`, the steal
  protocol, and the resume path: `parallel.py --resume` on a version 9
  checkpoint reads references like any other warm text.
- Certification: the shipped state only moves where a child's ADMM starts.
  A missing or corrupt store file should be a hard error on read, as a
  truncated inline token is today (`read_state` fails), not a silent cold
  start, so that a run never quietly loses its warm starts.

## Cloud coordination

The cloud session's stage-B tooling (the depth-2 frontier for wil50 and
the planned live-status path in `docs/qap-run-output.md`) parses state
files. Two things to agree before landing this:

1. Any reader of warm lines must treat a `<dnn>` field starting with `@` as
   a reference and copy `warm/` along with the file. The Python `State`
   class needs no change for parsing; only tools that decode tokens do.
2. The hash and file-name scheme above, so that a store written on one
   machine is valid on another (the names carry no host or pid).

## Steps and checks

1. `state.C`: read `@path`; write under `QAP_DNN_WARM_DIR`; `dnn_ref` on
   `WarmStart`; version 9 in `state.h`'s comment and `STATE_VERSION`.
   Check: the resume test with `QAP_DNN_SHIP=64` through the store is the
   same run bit for bit (as the inline test is today); a version 7 file
   still loads.
2. `parallel.py`: `--warm-dir`, the checkpoint link/copy, GC,
   `--checkpoint-inline`. Check: the coord test; nug30 with and without the
   store gives the same 3,716 nodes; checkpoint size and coordinator RSS on
   the tho40 depth-2 frontier (`experiments/tho40_estimate/depth2`), before
   and after.
3. `RUNNING.md`: the checkpoint-plus-`warm/` rule and the two new flags.

## The fallback without a solver change

If the state format must stay as it is, the coordinator alone can still
stream: keep `Warm` as a (hash, offset) pair into a spill file, write task
files and checkpoints by streaming warm texts from that file, and dedup by
hash. That fixes the coordinator's memory and the double-written parents
but not the worker I/O: every task file still carries its parents' tokens
in full, which is the larger cost at m = 40 and above.
