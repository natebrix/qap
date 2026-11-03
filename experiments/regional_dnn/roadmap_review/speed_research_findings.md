# Numerical speed and fleet audit: wil50 / sko49

Research against `qap_latest_research` commit ca50467, 2026-11-02 (America/Los_Angeles). No repository edits, installs, rented hardware, or remote writes. New experiments are isolated scratch copies. Distinguish **M** (measured), **C** (code observation), **A** (arithmetic conditional on measured shares), and **H** (hypothesis requiring a gate). The 76,755 versus 139,149 CPU-second subtree cut saving is already in the 2.2e9–7e9 baseline and is not multiplied in again.

## Main finding: preserve the actual same-node state between cut rounds

**This is a new, testable algorithmic lever absent from the roadmap.** `latest/dnncut.C:394` invokes `dnn_bound` afresh for each cut round and the quiet round. `latest/dnn.C:1260–1266` resets Y to zero except Y00=1 and reads only `start->Z`. `DNNWarm::Y` is saved at lines 1560–1567, but used only for separation. Every round discards the primal state and spectral tracker, even when the quiet round uses the identical objective as the best preceding round.

The negative result documented for **conditioning a parent's Y on a newly fixed assignment** does not cover same-node continuation. Here the dimension and feasible set are unchanged. First implement explicit same-node warm-Y semantics in the cut driver only, keeping the child warm-start path as-is. A second stage can retain the workspace, eigenbasis, refresh counter, objective scale and best-bound history between solves; an exact uninterrupted replay additionally needs the certificate shift state and check/give-up schedule. A best-round snapshot must match its multipliers and objective scale, not the final round's state. Changing L may require rescaling Z as current code does and validating whether the old basis remains useful.

**M, actual solver experiment:** A scratch copy of dnn.C differs only by an environment-switched copy of start->Y into Y. It links the installed SciPy OpenBLAS and the original LowerSum implementation. Nug20 root, PRSM gamma1=gamma=1, d3/L2, no early stopping, 200 iterations then 100:

| R-step | first 200 certified LB | restart with Z only, +100 | restart with Y,Z, +100 | uninterrupted 300 |
|---|---:|---:|---:|---:|
| tracked subspace | 2504.691671824542 | 2505.521291810191 | 2505.779138483394 | 2505.728789003385 |
| full single-precision spectrum | 2504.867417074598 | 2505.776906191506 | 2505.853335152050 | 2505.853335152050 |

A second **nug30 root/subspace** check gave first200=5938.534252345584; Z-only +100=5945.774320821572; YZ +100=5945.983810369371; uninterrupted300=5945.786097337961. This again supports a restart penalty without estimating wil50 gains. Full outputs are in `speed_continue_results.log`.

The full-spectrum YZ continuation matches uninterrupted bound and primal to printed precision exactly (primal 2499.603615416); Z-only restarts lose progress. With the subspace solver, the reset spectral tracker changes the trajectory, so YZ-only is not exact continuation. This validates the mechanism, **not a wil50 speedup**. Nug20 CPU times for the two full-spectrum continuation arms were 1.0354 and 1.0334 seconds: the benefit is iteration progress, not a faster iteration. The subspace arms differed in timing but were too noisy to infer a factor.

**Decisive next test:** replay saved same-node state at every round boundary on d1/d4/d6/d11 and clear misses, compare Z-only versus YZ versus full-state continuation, with identical cut policies and a fixed incumbent. First hold the cut pool/multipliers frozen to isolate restart loss; then let separation evolve. Record certified LB at equal work, CPU to fathom, failures, and total subtree CPU. Gate on at least 10–15% phase CPU reduction with no costly new failures; the only honest current range includes no benefit. If phase time drops by 10%, 25%, or 50%, run_ab3's 46% phase share gives whole-search factors **1.048, 1.130, 1.299**, respectively (A), before tree changes. Do not multiply this independently with shorter-round schedules: they affect the same iterations.

## Ranked follow-up candidates

| Priority | Candidate | What is known | Cheap decisive test / honest bound |
|---|---|---|---|
| 1 | Same-node YZ/full-state continuation above | Code proof plus small actual solver demonstration | Frozen-pool replay then paired cut subtrees. No wil50 factor established. |
| 2 | Adopt already measured d3/L2, then eliminate its dead Weyl pass | All32 CPU 5070 -> 4219 seconds, factor 1.20, same 18 fathoms; run_ab3 still used d4. Current code computes Weyl norm and Wprev copy immediately before overwriting lam_hi with the power estimate | Existing gain, not new research. Confirm on cut phases and sko49 before applying global factor. Guard Weyl work with tight_lam==0; check bitwise outputs. Approximately 16 ms/530 ms historical budget, at most ~3% if reproduced. |
| 3 | Compute residual bookkeeping only when consumed | Y/Z-step writes rs,ds and computes two dot products every iteration; pres/dY only read at checks unless Halpern/momentum is on | Skip on non-check iterations with accelerators off. Scratch extracted-loop experiment confirms bitwise Y/Z/Zc and check residual equality over 180 updates, but timing is noisy. At most the entire historical 30ms/530ms YZ share, hence <1.06 total, realistically a fraction. |
| 4 | AMX with explicit precision boundary and packing reuse | Intel officially supports BF16 packed GEMM. Fixed W is reused by d+1 large products per iteration, so split/pack W once; pack changing block operands per multiplication | Measure end-to-end relevant shapes including conversion, pack, Gram/QR/RR and full occupancy. Strong kernel candidate, no new measured solver factor. |
| 5 | Adaptive spectral accuracy / refresh schedule, preserving FP32 QR | d3/L2/r60 is untested; r120/Weyl was weaker. New angle is event-driven refresh based on block residual or certificate progress, rather than blindly long cadence | Compare checkpoint LB and CPU on both long quiet rounds and 150–300 iteration fathoms. Reduce refresh only after spectrum stabilizes, escalate when residual/negative count changes. No forecast; 7% pure refresh share bounds a pure refresh-saving gain to 1.075 absent iteration effects. |
| 6 | Workspace and memory reductions for full-host throughput | Current Work allocates double EVD workspace even on the subspace path; all four tracker matrices allocate full K² although active block is ~Kx600 | Lazy double workspace, capacity-sized tracker blocks with explicit growth. Record peak RSS and loaded socket throughput. This is primarily memory capacity/bandwidth opportunity, not predicted per-node acceleration. |
| 7 | GPU port only after a complete iteration throughput benchmark | GEMM peak is not iteration throughput; eigensolvers, host certificates and transfers constrain pipeline | Batched/multistream complete iteration at K=2117, kb=600; time 1/2/4/8/16 nodes with bounded workspaces, CUDA events plus wall time and host CPU. Do not select “8–16 streams” as a fact. |
| 8 | Asynchronous checkpoints plus measured slot placement | Existing runner really has a global stop-and-drain barrier; old roadmap already proposed leases | Remove barrier in a small runner first or follow plan. Compare completed-node throughput per physical host, not busy thread counts. This changes walltime/utilization, not work budget. |

## Important corrections to prior hardware claims

### Precision and AMX

1. The old 10-bit mantissa experiment is **200 iterations of rounded operands**, not convergence clearance. Its patch switches filter products, optionally W-times-Y in RR; it does not establish that Gram products, Cholesky QR, small eigensolve and all RR products tolerate this arithmetic. Gram matrices are especially sensitive to low precision and conditioning. Keep CholQR2/Gram and RR in FP32 initially, accelerate the big filter products first.
2. FP16 and TF32 share a 10-bit fraction, **not exponent range**. The rounding experiment preserves the FP32 exponent and does not emulate FP16 over/underflow. Chebyshev vectors can grow. A one-pass AMX-FP16 design needs per-block rescaling/range checks and a fallback, then convergence testing.
3. Two BF16 components per operand generally imply **three leading products**: AhBh + AhBl + AlBh. Merely saying “two-pass BF16 carries 16 bits” conflates components and matrix products. Omitting AlBl gives a normwise second-order approximation, not a relative elementwise guarantee under cancellation. Full FP32 emulation takes more work; Henry/Tang/Heinecke establish the component-product approach. Benchmark the exact chosen expansion.
4. The current environment is not the old one: `lscpu` reports Xeon Platinum 8573C, AVX512_BF16/FP16 but **no amx_tile/amx_bf16/amx_int8 flags**. Nine online CPUs, cgroup quota 800000/100000 = **eight CPU-seconds per wall second**. No GPU is presented. Therefore no honest AMX/GPU speed claim can be measured here.
5. oneMKL's packed BF16 compute API allows reuse of W's packed representation; Intel says packing and computing must use matching thread settings. Include packing cost amortized across the three/four filter products and RR, not a fresh conversion per product.

With accelerated share f=0.56, local kernel factors 2,4,8 imply whole-iteration factors **1.389,1.724,1.980**. Infinite speed on that share caps at **2.273** (A). If only the filter subset is accelerated, f is smaller. CPU trims and AMX overlap: never multiply their advertised whole-iteration factors without rebuilding a component cost model. The old 1.8–2.2 combined number is an unmeasured target, not evidence.

### Certificate validity boundary

The high-level certificate remains valid for arbitrary finite approximate iterates because certify recomputes V'ZcV and charges a Cholesky-based spectral penalty. The code actually checks the transformed multiplier in double; it does not blindly trust W-FF'. FP32 R can make it indefinite, causing a larger shift/refinement. Low-precision iteration work is therefore a sound design direction **if the existing certificate path and its rounding assumptions are preserved**.

However “any precision is safe” should not be used as permission to change the certificate transform, dot reductions, dpotrf or floating-point model unreviewed. `certify`'s margin assumes its double Householder transform with `64*m*u*||Z||` and a Rump-style Cholesky error allowance. A GPU implementation using a different transform, reduced-precision library math, FMA/reduction algorithm, or flush-to-zero policy needs the numerical bound re-established, or simply copy Zc back and run the present host certificate. Preserve symmetry and detect NaN/Inf. A small brute-force check does not itself prove the new rounding allowance.

### GPU throughput and host bottlenecks

A100's official dense TF32 peak is 156 TFLOP/s (312 is structured-sparse; this workload is dense). Peak gives a lower bound on GEMM time, not an achieved timestep. The old 3–8ms *aggregate node-iteration* estimate needs measured complete kernels and concurrency; latency per node and aggregate time per node must not be conflated.

The roadmap proposes moving both the 2117-square refresh and certificate to 12–24 host vCPUs. At old measured costs that consumes per node-iteration:

- refresh: 0.734/20 = **0.0367 host-core seconds**;
- certificate: 0.143/10 = **0.0143 host-core seconds**;
- sum **0.051 host-core seconds**, before the ~600-square RR eigensolve, transfers, separation, dispatch and bookkeeping.

Thus 125–333 node-iterations/s (8–3ms aggregate) needs **6.4–17.0 fully available reference host cores** for those two tasks alone. If 12–24 vCPUs are SMT threads corresponding to 6–12 physical cores, the fast end is already impossible at those host costs. Extra streams cannot remove this throughput limit. On-GPU EVD and certificate reduce it but need measurement and numerical validation. Also query actual device/host workspace requirements; generic cuSOLVER dispatch does not establish performance at our shapes.

The earlier 0.57GB node-state estimate is stale for current code. Direct allocation count at m=47 is **611MB** for dominant Work arrays plus **78MB** per warm YZ state, excluding vector overhead, cuts, best snapshots, BLAS workspaces and any simultaneous solver state. Several warm snapshots in cut rounds matter. GPU memory fitting many nodes does not show they should all execute concurrently.

100 GPUs in one week at 95% useful utilization would require **38.3–121.8 reference-core equivalents per GPU** for the current 2.2e9–7e9 work. With 0.53s per reference iteration, these correspond to **13.8–4.35ms aggregate node-iteration** thresholds. This is useful benchmark acceptance arithmetic, not a prediction and not a reduction in algorithmic reference-core work.

## Fleet and accounting audit

`dnn.C:221` measures CLOCK_THREAD_CPUTIME_ID. With one BLAS thread this is useful, although SMT and host contention affect what a “reference core” means. It omits helper-thread CPU when BLAS uses multiple threads, and GPU work/asynchronous waits. Preserve **wall elapsed**, **whole-process plus child CPU**, device elapsed/occupancy, hardware identity, physical-core/vCPU mapping, BLAS thread count and completed decisions. Otherwise an accelerator can appear to eliminate CPU work merely because it is not counted.

run_ab3 cut arm recorded 76755.31 CPU s /21933.11 wall s = 3.50 occupied CPU threads, **87.5% of four slots**, not evidence for a universal 70–85% utilization. Base arm ~3.90/4. CPU/wall does not diagnose whether losses are checkpoints, tail, throttling or waits, and 1000-worker tails must be measured. At 95% efficiency a 1000-core week supplies **574,560,000 core-seconds**, not 600 million. Reducing 2.2e9–7e9 to fit therefore needs **3.83–12.18x** work/throughput improvement; at ideal efficiency it is 3.64–11.57x.

Do not charge a coordinator's utilization improvement against a baseline already assuming 100% busy. Checkpoint elimination can reduce actual CPU work only when it prevents lost warm starts/recomputation, which must be separately measured. The per-node checkpoint/lease granularity must accommodate long cut phases and work stealing without repeatedly stripping the new cut/YZ state.

Socket contention: the old transform “6GB/s” is per isolated streaming phase. Multiplying it by 64 processes as if each process streamed constantly is not a valid average bandwidth estimate: transforms occupy only a fraction of each iteration; GEMMs themselves also access DRAM/cache. Measure steady-state node throughput across 1,2,4,... physical cores and then SMT with CPU pinning/NUMA placement, and pick workers per host by throughput/$ and RSS. Nodes have independent matrices, so GEMM reuse across nodes is absent, but staggered phases may improve utilization. Concurrent vs synchronous cut rounds change burstiness. Current virtual topology even reports inconsistent nine CPUs/five L2 instances; do not infer physical core counts from that VM alone.

## Sources checked directly (primary)

- Intel BF16 GEMM reference: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/cblas-gemm-bf16bf16f32.html
- Intel packed BF16 API: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/cblas-gemm-bf16bf16f32-compute.html
- Intel moderator confirmation AMX dispatch and packing: https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/AMX-bf16-use-in-MKL-CBLAS/m-p/1641052?profile.language=en
- Henry, Tang, Heinecke, ARITH26 paper: https://arxiv.org/abs/1904.06376
- NVIDIA cuBLAS precision, Tensor Core alignment and reproducibility documentation: https://docs.nvidia.com/cuda/cublas/
- NVIDIA cuSOLVER API and device-memory/host-workspace model: https://docs.nvidia.com/cuda/cusolver/index.html
- NVIDIA A100 specifications: https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/a100/pdf/a100-80gb-datasheet-update-nvidia-us-1521051-r2-web.pdf

Relevant retrieval IDs for root to open before final citations: turn39search0, turn39search2, turn39search3, turn39search4, turn40view0, turn40view1, turn43search24.

## Reproduction artifacts

Scratch only: `speed_continue.cpp`, `speed_dnn_experiment.C`, `speed_continue_results.log`, executable `speed_continue`, `speed_yz_microbench.cpp` and executable `speed_yz_microbench` alongside this report. Compile flags `-O3 -march=native -ffp-contract=off -DBNB_STANDALONE -DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS`; link existing scipy.libs/libscipy_openblas-6cdc3b4a.so. Harness directly constructs MAT views of QAPLIB matrices and runs the original certified bound. The microbenchmark timings varied significantly; no speedup extrapolated from them. No solver files were edited and no full-instance search was attempted.
