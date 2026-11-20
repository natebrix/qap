"""Q4: the same node under two or more lines, from every trace that has one.  -> out/q4_obs.tsv, out/q4.txt

An observation is (instance, node as a SET of fixes, line, b, source): a branched node whose children are all in
that source's trace. Sources (each linked on its own): production traces, the c47 descent, every linetest arm,
tai60b_lines' depth-2 arms, the tai50b runs (rule 2 bottom, rule 9, the signal run, gz0nm), tai35b/tai40b's
line_signals arms. Nodes seen under different lines are the counterfactuals; nodes seen twice under the SAME line
(different sources: cold/strong parents, other builds) measure how much b moves without a change of line.

Line features from the instance and the node alone (fixconv.py: fix f:l means facility f, the index of the first
matrix A, at location l):
  column c_j: ocl = A-cluster of the incumbent's facility at j; cand = free facilities of that cluster
  row r_i:    ocl = i's cluster; cand = free locations whose incumbent facility is in i's cluster
  tight       mean relative L1 distance between two members' rows of A (clusters.py), small = interchangeable
"""
import os, glob, json, itertools
from collections import defaultdict
import numpy as np
from common import *

CL = json.load(open('out/clusters.json'))
LS = f'{WT}/experiments/line_signals/runs'


def sources():
    S = {}
    S[('tai60b', 'prod')] = RUNS['tai60b-261117']['files']
    S[('tai60b', 'c47down')] = RUNS['tai60b-c47']['files']
    for d in sorted(glob.glob(f'{LT}/runs/scan/*') + glob.glob(f'{LT}/runs/arm17/*') + glob.glob(f'{LT}/runs/d9x/*')):
        S[('tai60b', 'lt:' + d.split('/runs/')[1])] = [d + '/nodes.txt']
    for a in ['sc54', 'c54', 'c47', 'r31', 'kids']:
        S[('tai60b', 'lt:' + a)] = [f'{LT}/runs/{a}/nodes.txt']
    S[('tai60b', 'lt:sc54')] += sorted(glob.glob(f'{LT}/runs/lvl/c54r9/nodes/*'))   # the strong c54 arm's depth-8 lines (rule 9) and their children (lvl.py)
    for a in ['c26', 'r1', 'r45', 'r8', 'rule9']:
        S[('tai60b', 'tai60b_lines:' + a)] = [f'{WT}/experiments/tai60b_lines/runs/{a}/nodes.txt']
    S[('tai50b', 'prod-r2')] = RUNS['tai50b-261116']['files']
    S[('tai50b', 'r9')] = RUNS['tai50b-r9']['files']
    S[('tai50b', 'ls')] = RUNS['tai50b-ls']['files']
    g = glob.glob(f'{WT}/infra/laptop/runs/tai50b-gz0nm/nodes.txt*')
    if g:
        S[('tai50b', 'gz0nm')] = g
    for a in ['tai35b_base', 'tai35b_base2', 'tai35b_la1p', 'tai35b_la4d2', 'tai35b_la4p']:
        S[('tai35b', a)] = [f'{LS}/{a}/nodes.txt.gz']
    for a in ['tai40b_base', 'tai40b_la4p']:
        S[('tai40b', a)] = [f'{LS}/{a}/nodes.txt.gz']
    return {k: [f for f in v if os.path.exists(f)] for k, v in S.items()}


def read_source(files):
    nodes = {}
    for f in files:
        with opener(f) as fh:
            for line in fh:
                t = line.split()
                if len(t) >= 12 and t[0] == 'node':
                    n = Node(t)
                    nodes[n.key] = n
    for n in nodes.values():
        if n.fixes:
            pk = ','.join(f'{a}:{b}' for a, b in n.fixes[:-1]) or '-'
            if pk in nodes:
                nodes[pk].kidsl.append(n)
    return nodes


class Inst:
    def __init__(self, inst):
        self.n, self.A, self.B = read_dat(inst)
        self.obj, self.p = read_sln(inst)
        self.q = [0] * self.n
        for f, l in enumerate(self.p):
            self.q[l] = f
        self.fcl = {}
        self.cls = CL[inst]['fac_clusters']
        for k, c in enumerate(self.cls):
            for f in c:
                self.fcl[f] = k
        m = np.array(self.A, dtype=float); mean = m.mean()
        self.tight = {}
        for k, c in enumerate(self.cls):
            if len(c) > 1:
                self.tight[k] = float(np.mean([np.abs(np.delete(m[i], [i, j]) - np.delete(m[j], [i, j])).mean() / mean
                                               for i, j in itertools.combinations(c, 2)]))
            else:
                self.tight[k] = float('nan')

    def feat(self, fixes, line):
        side, idx = line[0], int(line[1:])
        ff = set(f for f, l in fixes); fl = set(l for f, l in fixes)
        if side == 'c':
            ocl = self.fcl[self.q[idx]]
            cand = sum(1 for f in self.cls[ocl] if f not in ff)
        else:
            ocl = self.fcl[idx]
            cand = sum(1 for l in range(self.n) if l not in fl and self.fcl[self.q[l]] == ocl)
        return side, ocl, cand, self.tight[ocl]


def main():
    S = sources()
    obs = []
    inst_cache = {}
    for (inst, src), files in S.items():
        if not files:
            continue
        nodes = read_source(files)
        I = inst_cache.setdefault(inst, Inst(inst))
        for n in nodes.values():
            if n.out != 'B' or n.kids is None or len(n.kidsl) != n.kids:
                continue
            b = sum(1 for c in n.kidsl if c.out == 'B')
            side, ocl, cand, tight = I.feat(n.fixes, n.line)
            fs = ','.join(f'{a}:{c}' for a, c in sorted(n.fixes)) or '-'
            path = int(all(I.p[f] == l for f, l in n.fixes))
            obs.append(dict(inst=inst, src=src, node=fs, depth=n.depth, line=n.line, b=b, nkids=len(n.kidsl),
                            nC=sum(1 for c in n.kidsl if c.out == 'C'), short=n.short, side=side, ocl=ocl,
                            cand=cand, tight=tight, path=path, level=n.level))
    with open('out/q4_obs.tsv', 'w') as f:
        keys = list(obs[0].keys())
        f.write('\t'.join(keys) + '\n')
        for o in obs:
            f.write('\t'.join('' if o[k] is None else (f'{o[k]:.4f}' if isinstance(o[k], float) else str(o[k]))
                              for k in keys) + '\n')
    print(len(obs), 'observations')


if __name__ == '__main__':
    main()
