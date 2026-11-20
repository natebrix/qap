"""One row per branched node, for every run in common.RUNS: out/bn_<run>.tsv.gz.

Columns (all from the trace, plus the instance's best known solution and clusters):
 run depth key rule level line side idx  nkids b nD nC nQ nP nE  sub_n sub_cpu cpu
 short rel_short qpb_gap path mirror ham  parent_key parent_line  bshorts (the branching children's shortfalls,
 ';'-joined, descending)  bfree (for a column: the branching children's facilities; for a row: their locations)
 bclus (cluster of each branching child's facility, column lines only on tai-b) nclus_b maxclus_b
 near_fath (children fathomed by D or C within 5% of the node's own shortfall above the cutoff... see below)
 lsig (tai50b-ls: 1 where the landscape signal chose the line)

short       cutoff - the node's best bound (max of qpb, dnn, cut)
rel_short   short / the root's short
qpb_gap     (cutoff - qpb) / cutoff
path        1 if every fix agrees with the best known permutation (fix f:l means p[f] = l; fixconv.py)
mirror      tai60b: 1 if every fix agrees with the permutation with facilities 43 and 44's locations exchanged
ham         number of the node's fixes that disagree with the best known permutation (min over p and mirror)
near_fath   children fathomed (D, C) whose bound exceeded the cutoff by less than 1% of the parent's short
"""
import gzip, json, sys
from common import *

CL = json.load(open('out/clusters.json'))
ROUNDN = lambda x: '' if x is None else f'{x:.0f}'


def ls_chosen():
    s = set()
    for line in open(f'{LAP}/runs/tai50b-ls/ls.log'):
        if 'chose' not in line:
            continue
        fx = line.split('fixes')[1].split('|')[0].split()
        s.add(','.join(fx) if fx else '-')
    return s


def build(run):
    R = RUNS[run]
    inst = R['inst']
    nodes, tops = load(run)
    try:
        obj, p = read_sln(inst)
    except Exception:
        p = None
    pm = None
    q = None
    if p:
        q = [0] * len(p)
        for f, l in enumerate(p):
            q[l] = f
    if inst == 'tai60b' and p:
        pm = list(p); pm[43], pm[44] = p[44], p[43]
    fcl = {}
    if inst in CL:
        for k, c in enumerate(CL[inst]['fac_clusters']):
            for f in c:
                fcl[f] = k
    root = [t for t in tops if t.depth == 0]
    root_short = root[0].short if root else max(t.short for t in tops if t.short is not None)
    if run == 'tai60b-c47':
        root_short = 2899357.0   # tai60b-261117's root (same cutoff), so rel_short compares with production
    lsset = ls_chosen() if run == 'tai50b-ls' else set()
    rows = []
    for n in nodes.values():
        if n.out != 'B':
            continue
        cnt = {k: 0 for k in 'BDCQPE'}
        for c in n.kidsl:
            cnt[c.out] = cnt.get(c.out, 0) + 1
        side, idx = n.line[0], int(n.line[1:])
        bk = sorted([c for c in n.kidsl if c.out == 'B'], key=lambda c: -(c.short or 0))
        bfree = []
        for c in bk:
            f, l = c.fixes[-1]
            bfree.append(f if side == 'c' else l)
        bclus = [fcl.get(f, -1) for f in bfree] if (side == 'c' and fcl) else []
        nclus = len(set(bclus)) if bclus else ''
        maxclus = max([bclus.count(k) for k in set(bclus)]) if bclus else ''
        # cluster features (tai-b): free facilities per cluster at this node; for a column c_j the cluster of the
        # incumbent's (= best known's) facility at j and how many of that cluster are free; for a row r_i the
        # cluster of i, how many free locations host a facility of i's cluster in the incumbent, and how many
        # of the branching children are such locations
        ocl = fopt = hostfree = bhost = dclus = fdom = ''
        if fcl and p:
            fixedf = set(f for f, l in n.fixes); fixedl = set(l for f, l in n.fixes)
            freec = {}
            for f, k in fcl.items():
                if f not in fixedf:
                    freec[k] = freec.get(k, 0) + 1
            if side == 'c':
                ocl = fcl[q[idx]]; fopt = freec.get(ocl, 0)
                if bclus:
                    dclus = max(set(bclus), key=bclus.count); fdom = freec.get(dclus, 0)
            else:
                ocl = fcl[idx]
                hostfree = sum(1 for l in range(len(q)) if l not in fixedl and fcl[q[l]] == ocl)
                bhost = sum(1 for l in bfree if fcl[q[l]] == ocl)
        sh = n.short
        near = 0
        if sh and sh > 0:
            for c in n.kidsl:
                if c.out in 'DC' and c.short is not None and -c.short < 0.01 * sh:
                    near += 1
        if p:
            path = int(all(p[f] == l for f, l in n.fixes))
            mir = int(bool(pm) and all(pm[f] == l for f, l in n.fixes))
            ham = sum(1 for f, l in n.fixes if p[f] != l)
            if pm:
                ham = min(ham, sum(1 for f, l in n.fixes if pm[f] != l))
        else:
            path = mir = ham = ''
        rows.append([run, n.depth, n.key, rule_of(run, n), n.level, n.line, side, idx,
                     len(n.kidsl), n.kids, cnt['B'], cnt['D'], cnt['C'], cnt['Q'], cnt['P'], cnt['E'],
                     n.sub_n, f'{n.sub_cpu:.1f}', f'{n.cpu:.1f}',
                     ROUNDN(sh), f'{sh / root_short:.4f}' if sh is not None else '',
                     f'{(n.cutoff - n.qpb) / n.cutoff:.5f}' if (n.qpb is not None and n.cutoff) else '',
                     path, mir, ham,
                     n.parent.key if n.parent else '', n.parent.line if n.parent else '',
                     ';'.join(ROUNDN(c.short) for c in bk), ';'.join(map(str, bfree)), ';'.join(map(str, bclus)),
                     nclus, maxclus, near, int(n.key in lsset), ocl, fopt, dclus, fdom, hostfree, bhost])
    hdr = ('run depth key rule level line side idx nkids kidsf b nD nC nQ nP nE sub_n sub_cpu cpu short rel_short '
           'qpb_gap path mirror ham parent_key parent_line bshorts bfree bclus nclus_b maxclus_b near_fath lsig ocl fopt dclus fdom hostfree bhost').split()
    with gzip.open(f'out/bn_{run}.tsv.gz', 'wt') as f:
        f.write('\t'.join(hdr) + '\n')
        for r in sorted(rows, key=lambda r: (r[1], r[2])):
            f.write('\t'.join(map(str, r)) + '\n')
    tot_n = len(nodes)
    tot_cpu = sum(x.cpu for x in nodes.values())
    print(f'{run}: {tot_n} nodes, {tot_cpu / 3600:.1f} CPU h, {len(rows)} branched, tops {len(tops)} '
          f'(depths {sorted(set(t.depth for t in tops))[:5]}), root short {root_short:.0f}')


if __name__ == '__main__':
    for run in (sys.argv[1:] or list(RUNS)):
        build(run)
