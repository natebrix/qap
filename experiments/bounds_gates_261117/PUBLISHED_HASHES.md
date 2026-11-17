# Published commits for the hashes these notes cite

The four research branches were published with rewritten history: every commit was recreated
(author times 11:53-12:02, Nov 17 2026). So the hashes the notes cite are not in the published
repository. They appear in `README.md`, `RUN_STATUS.md`, `REPLAY.md`,
`verification/certificates.json`, `seed_capture/*/run.json` and the experiment branches' own READMEs.
This table maps each one to its published commit.

How each was matched:

- **checksum:** `verification/certificates.json` records the SHA-256 of the six source files of each
  build it checked (`build.source_sha256`), and the published commit's files have those SHA-256s.
- **subject:** only the commit message and the files the commit changes match.
- **branch head:** the published branch's last commit. No checksum covers a whole tree.

| Cited | What | Published | Branch | Matched by |
|---|---|---|---|---|
| `6e53f4ff` (`6e53f4ff38710ce1...`) | coarse-label cuts: the measured solver | `f4487612` | coarse-label-cuts | checksum, all 6 files. No later published commit changes them. |
| `26260591` | coarse-label cuts: final head | `505fc306` | coarse-label-cuts | branch head |
| `d38cacbb` | the coarse `dnncuts.C` used as the tai60b checkpoint reader (`experiments/cp_scaling_261117/README.md` on cp-congruence) | `experiments/dnn_cuts/dnncuts.C` at `f4487612` | coarse-label-cuts | not verifiable: no checksum recorded. The published branch has one version of this file, the one the `6e53f4ff` build used (SHA-256 `105b1a22...`). |
| `f6b2bbe0` (`f6b2bbe0aeafd8dd...`) | CP congruence scaling: primary implementation | `e7e11bb8` | cp-congruence | checksum, all 6 files (from `1e119c63` on, `latest/dnn.C` differs) |
| `d175ac1f` | CP scaling: homogenizing-weight fallback | `1e119c63` | cp-congruence | subject. It is the only later published commit that changes solver code (`latest/dnn.C`, `tests/dnn_continuation.C`). |
| `79c5dc37` | CP scaling: final head | `3395b8bb` | cp-congruence | branch head |
| `0fc7508b` | regional equalities: primary implementation | `c8e79e32` | cp-region-equalities | subject. Its `latest/settings.C` has the `f734af31` build's checksum and its `latest/dnncut.C` does not, as expected when the fix changes `dnncut.C` only. |
| `f734af31` (`f734af31b50893a1...`) | regional equalities: imported-region capacity fix | `215e07cb` | cp-region-equalities | checksum, all 6 files |
| `57e6a285` | regional equalities: final head | `5c35a805` | cp-region-equalities | branch head |
| `7fafc560` (`7fafc560e436f12c...`) | the seed build for the shared warm states (`REPLAY.md`, `seed_capture/*/run.json`) | never published, by design | none | cognac `7fdac382` plus cognac's `experiments/tai60b_lines/force_line.patch` (added in `2a5a2e77`), whose SHA-256 `b3c0e372...` is the one `REPLAY.md` records. `REPLAY.md`'s commands recreate it as a commit with another hash and the same `latest/bnb.C`. |

`9d030c2ecf7f` in the logs and manifests (`/workspace/scratch/9d030c2ecf7f/...`) names a scratch
directory, not a commit. `ba47c211` and `7fdac382` are cognac commits and resolve as cited.
