"""t0.py NAME [FIXES ...]: at the root (FIXES '-') and at given nodes of NAME: the AJ LP (HiGHS IPM, no crossover,
and the Lagrangian bound of its dual), the converged DNN (dnntest -H 0, 20,000 iterations, no cutoff: the solver's
certified value without the cut phase), the solver's QPB (bnbdrv node trace, PARAM/default.par) and GLB (ajlp.glb).
Appends one line per node to out/t0_NAME.tsv."""
import sys, os, re, subprocess, time, resource, numpy as np, ajlp
WT = '/home/user/wt-lp/latest'
name = None


def dnn(fixes, iters=20000):
    cmd = [f'{WT}/dnntest', '-H', '0', '-i', str(iters), '-c', '1e15']
    if fixes != '-': cmd += ['-x', fixes]
    o = subprocess.run(cmd + [f'{WT}/../qaplib/{name}.dat'], capture_output=True, text=True).stdout
    m = re.search(r'RESULT .*lb=(\S+) primal=(\S+) status=(\S+) iters=(\d+) cpu=(\S+)', o)
    return float(m.group(1)), float(m.group(2)), m.group(3), int(m.group(4)), float(m.group(5))


def qpb(fixes):
    tr = f'/tmp/claude-0/qpbtr_{os.getpid()}.txt'
    if os.path.exists(tr): os.remove(tr)
    d = 0 if fixes == '-' else len(fixes.split(','))
    cmd = [f'{WT}/bnbdrv', '-p', f'{WT}/PARAM/default.par', '-u', '1e15', '-d', str(d)]
    if fixes != '-': cmd += ['-F', fixes]
    subprocess.run(cmd + [f'{WT}/../qaplib/{name}.dat'], capture_output=True, text=True, env=dict(os.environ, QAP_NODE_TRACE=tr), cwd=WT)
    ln = open(tr).readline().split()
    return float(ln[3])


def main():
    global name
    name = sys.argv[1]; nodes = sys.argv[2:] or ['-']
    A, B, opt, perm = ajlp.inst(name)
    os.makedirs('out', exist_ok=True)
    out = open(f'out/t0_{name}.tsv', 'a')
    if os.path.getsize(f'out/t0_{name}.tsv') == 0:
        out.write('inst\tfixes\tm\topt_node\tGLB\tQPB\tLP\tLP_lagr\tDNN\tDNN_primal\tDNN_iters\tDNN_cpu\tLP_cpu\tLP_status\tipm_iters\n')
    for fx in nodes:
        Ar, Br, C, s, I, J = ajlp.reduce(A, B, fx)
        m = len(I)
        # the node's optimum: global opt if the optimal permutation satisfies the fixes, else '-'
        F = ajlp.parse_fixes(fx)
        on = opt if all(perm[r] == v for r, v in F) else '-'
        lp = ajlp.LP(Ar, Br, C); r = lp.solve('ipm')
        lag, _ = lp.lagr(r['u'])
        g = ajlp.glb(Ar, Br, C) + s
        q = qpb(fx)
        dv = dnn(fx)
        row = [name, fx, m, on, f'{g:.1f}', f'{q:.2f}', f'{r["obj"] + s:.2f}', f'{lag + s:.2f}', f'{dv[0]:.2f}', f'{dv[1]:.2f}',
               dv[3], f'{dv[4]:.1f}', f'{r["cpu"]:.1f}', r['status'], r['iters'][1]]
        out.write('\t'.join(map(str, row)) + '\n'); out.flush()
        print('\t'.join(map(str, row)), flush=True)


if __name__ == '__main__':
    main()
