"""Why the test must be strict: count nodes whose best pair min Delta is exactly 0 (a non-strict test would drop
them) and that hold an optimum; and for symmetric instances, an explicit pair of twin nodes that would drop each other.
Output: nonstrict.out"""
import itertools, numpy as np
from qapdata import read_dat
from dom import Inst
from validate import perms_of, all_costs
for nm in ['esc08a', 'esc08e', 'nug08', 'had10']:
    A, B = read_dat(nm); inst = Inst(A, B); n = inst.n
    P = perms_of(None, None, {}, n); costs = all_costs(A, B, P); opt = costs.min()
    zero_opt = 0; tested = 0; twins = None
    # all depth-2 nodes {0: a, 1: b} and depth-3 nodes {0:a,1:b,2:c}
    for d in (2, 3):
        for cols in itertools.permutations(range(n), d):
            F = {r: c for r, c in zip(range(d), cols)}
            mins, pairs = inst.pair_mins(F); tested += 1
            sel = np.all(P[:, :d] == np.array(cols, dtype=np.int8), 1)
            if mins.max() == 0 and costs[sel].min() == opt:
                zero_opt += 1
                if twins is None and d == 2:
                    G = {0: cols[1], 1: cols[0]}
                    m2, _ = inst.pair_mins(G)
                    sel2 = np.all(P[:, :2] == np.array([cols[1], cols[0]], dtype=np.int8), 1)
                    twins = (F, G, int(m2.max()), bool(costs[sel2].min() == opt))
    print(f'{nm}: nodes tested {tested}; best pair min Delta == 0 and holding an optimum: {zero_opt}; '
          f'example twin pair (F, its exchange image, image best min, image holds an optimum): {twins}')
