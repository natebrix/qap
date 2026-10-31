"""Reconstruct the two fixed-multiplier triangle experiments from a checkpoint."""
import argparse,json
from pathlib import Path
import numpy as np
from region_cuts import node_indices,triangle_separate,cut_matrix
ap=argparse.ArgumentParser();ap.add_argument('prefix');ap.add_argument('--fixes',required=True)
ap.add_argument('--keep',type=int,default=100);ap.add_argument('--multipliers',default='2,10')
args=ap.parse_args();rows,cols=node_indices(args.fixes);m=len(rows);N=m*m+1
Y=np.fromfile(args.prefix+'.Y.bin').reshape(N,N)
cuts,stats=triangle_separate(Y,rows,keep=args.keep);L=np.zeros((N,N));h=0
for c in cuts:
    G,hh=cut_matrix(c,m);L+=G;h+=hh
for raw in args.multipliers.split(','):
    mu=float(raw);path=args.prefix+f'_tri{args.keep}_mu{mu:g}.L.bin';(mu*L).tofile(path)
    print(path,'shift_add',-mu*h)
Path(args.prefix+f'_tri{args.keep}_spec.json').write_text(json.dumps(dict(cuts=cuts,sum_h=h,
   sum_initial_violation=sum(c['violation'] for c in cuts)),indent=2))
