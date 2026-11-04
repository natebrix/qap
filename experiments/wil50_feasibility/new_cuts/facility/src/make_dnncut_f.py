#!/usr/bin/env python3
"""Make dnncut_f.C from a copy of latest/dnncut.C (dnncut_f.C.orig): the same
code with the public function renamed dnn_cut_rounds_f and
 (1) an export hook (CutExport, fdrv.h): the round-0 iterate, and at the quiet
     round the triangle-only Lagrangian term at the multipliers the quiet round
     uses, its constant, and the ADMM state it starts from; optionally stop there;
 (2) extra regions (g_fdrv_regions, fdrv -regions FILE [-nogrid]) on either
     side, each region carrying its own side, so location regions other than the
     grid prefixes (parity classes) and facility regions can join the pool. The
     cut algebra (aggregate, add_rcut, the four triangle types) is unchanged; only
     the lifted index map is taken per region."""
s = open('dnncut_f.C.orig').read()


def rep(old, new):
    global s
    assert s.count(old) >= 1, old[:60]
    s = s.replace(old, new, 1)


# drop the DNNCutParams constructor (it stays in latest/dnncut.o)
i = s.index('DNNCutParams::DNNCutParams()'); j = s.index('}\n', i) + 2
s = s[:i] + s[j:]
rep('#include "dnncut.h"', '#include "dnncut.h"\n#include "fdrv.h"')
rep('DNNCutResult dnn_cut_rounds(MAT *A', 'DNNCutResult dnn_cut_rounds_f(MAT *A')
rep('int n_orig, const DNNCutParams &cp, DNNWarm *warm, double first_lb)',
    'int n_orig, const DNNCutParams &cp, DNNWarm *warm, double first_lb, CutExport *ex)')
rep('    const std::vector<double> &Y = warm->Y;\n',
    '    const std::vector<double> &Y = warm->Y;\n    if (ex && r == 0) ex->Y0 = warm->Y;\n')
hook = '''    if (ex) {
      ex->Ltri.assign((size_t)n1 * n1, 0.0); ex->hsum_tri = 0; ex->abs_tri = 0; ex->active = 0;
      for (const Cut &t : cuts)
        if (t.mu > 0) {
          ex->abs_tri += t.s < 0 ? add_cut(t, ex->Ltri, n1, t.mu) : add_rcut(t, regions[t.s], ex->Ltri, lift, n1, t.mu);
          ex->hsum_tri += t.mu * cut_rhs(t); ex->active++;
        }
      ex->state = *warm; ex->best_mu_lb = best_mu_lb; ex->have = 1;
      if (ex->stop_before_quiet) { if (out.best_lb > cutoff) out.status = DNN_FATHOMED; return out; }
    }
    DNNResult res = solve(1, cp.freeze_iters);
    if (ex) { ex->final_state = *warm; ex->quiet_lb = res.lb; ex->quiet_primal = res.primal; }
'''
rep('    DNNResult res = solve(1, cp.freeze_iters);\n', hook)

# (2) extra regions
rep('struct Region { std::string name; std::vector<int> loc; };',
    'struct Region { std::string name; std::vector<int> loc; int side = -1; };')
rep('void aggregate(const Region &S, const std::vector<double> &Y, const Lift &lift, int n1,',
    'void aggregate(const Region &S, const std::vector<double> &Y, const Lift &lift0, int n1,')
rep('  const int m = lift.m;\n  p.assign(m, 0.0);',
    '  const Lift lift = {lift0.m, S.side < 0 ? lift0.side : S.side};\n  const int m = lift.m;\n  p.assign(m, 0.0);')
rep('double add_rcut(const Cut &t, const Region &S, std::vector<double> &L, const Lift &lift, int n1,\n                double mu)\n{\n',
    'double add_rcut(const Cut &t, const Region &S, std::vector<double> &L, const Lift &lift0, int n1,\n                double mu)\n{\n'
    '  const Lift lift = {lift0.m, S.side < 0 ? lift0.side : S.side};\n')
rep('}  // namespace\n\nDNNCutResult dnn_cut_rounds_f(',
    '}  // namespace\n\nFdrvRegions g_fdrv_regions;\n\nDNNCutResult dnn_cut_rounds_f(')
extra = '''    else if (cp.prefix) regions = prefix_regions(m);
    if (g_fdrv_regions.nogrid) regions.clear();
    for (const FdrvRegion &e : g_fdrv_regions.extra) {
      // original indices of side e.side -> local indices; skip empty and full sets
      const std::vector<int> &fo = g_fdrv_regions.free_orig[e.side];
      Region R; R.name = e.name; R.side = e.side;
      for (size_t t = 0; t < fo.size(); t++)
        if (std::find(e.orig.begin(), e.orig.end(), fo[t]) != e.orig.end()) R.loc.push_back((int)t);
      if (R.loc.empty() || (int)R.loc.size() >= m) continue;
      regions.push_back(R);
    }
'''
rep('    else if (cp.prefix) regions = prefix_regions(m);\n', extra)
open('dnncut_f.C', 'w').write(s)
print('ok')
