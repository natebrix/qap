# nug30: does today's cut phase beat the record recipe? (bound_lab agent `nug30_phase`, Nov 16 2026)

IN PROGRESS. Baseline measured [m] (`runs/base/log.txt`): the record recipe (rule 9, DNN at depths 1-7, give-up 200,
no cut phase), whole search, 4 workers on this box (`lscpu.txt`): 3,688 nodes, 110,926 DNN iterations, 3,030 CPU s,
771 s wall. Cut arms z003/z010/z030/z060 running (`queue.log`); `python3 summarize.py` prints the table.
