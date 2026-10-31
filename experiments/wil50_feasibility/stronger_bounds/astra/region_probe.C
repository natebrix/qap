// Research driver: native row-major X (A indices first, B indices second).
// Usage: region_probe instance fixes prefix iterations [Ladd.bin shift_add warm.bin]
#include "dnn.h"
#include "qap.h"
#include "assign.h"
#include "util.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
int main(int argc, char **argv) {
  if (argc < 5) return 2;
  setvbuf(stdout,NULL,_IOLBF,0);
  dnn_set_blas_threads(getenv("PROBE_THREADS")?atoi(getenv("PROBE_THREADS")):1);
  QAP *q = readQAPLIB(argv[1]); if (!q) return 2;
  set_fathom_tolerance(q); symmetrize_QAP(q);
  int n=q->A->m; QAPAssignment node(n); std::string s=argv[2];
  if (s!="-") for(size_t p=0;p<s.size();) {
    size_t e=s.find(',',p); if(e==std::string::npos)e=s.size();
    int i,j; if(sscanf(s.substr(p,e-p).c_str(),"%d:%d",&i,&j)!=2 ||
       i<0 || j<0 || i>=n || j>=n || node.p[i]!=NONE || node.pinv[j]!=NONE) return 2;
    node.fix(i,j);p=e+1;
  }
  QAP qr; reduceQAP(q,&qr,node); int m=n-node.nfix,N=m*m+1;
  if(m<2){fprintf(stderr,"This generalized-objective research probe requires residual order >= 2.\n");delete q;return 2;}
  size_t sz=(size_t)N*N; std::string pre=argv[3];
  DNNParams p; p.max_iter=atoi(argv[4]);p.hopeless_after=0;
  p.gamma=1;p.gamma1=1;p.cheb_degree=4;p.verbose=1;
  if(getenv("PROBE_ADMM")){p.gamma=1.618;p.gamma1=0;}
  p.inner_lap=DNN_LAP_EXACT; p.check_every=25;
  // Double eigensolver is available for tiny-instance validation.
  if(getenv("PROBE_DOUBLE"))p.eig_mode=DNN_EIG_DOUBLE;
  double shiftadd=0;std::vector<double> L;DNNWarm start,save;
  if(argc>5 && std::string(argv[5])!="-") {
    if(argc<7){fprintf(stderr,"Ladd requires shift_add.\n");delete q;return 2;}
    L.resize(sz);FILE*f=fopen(argv[5],"rb");if(!f)return 2;
    if(fread(L.data(),sizeof(double),sz,f)!=sz)return 2;fclose(f);
    p.L_add=L.data();for(double v:L)p.L_add_abs+=fabs(v);
    shiftadd=atof(argv[6]);
  }
  if(argc>7 && std::string(argv[7])!="-") {
    FILE*f=fopen(argv[7],"rb");if(!f)return 2;start.m=m;start.Z.resize(sz);
    if(fread(&start.scale,sizeof(double),1,f)!=1)return 2;
    if(fread(start.Z.data(),sizeof(double),sz,f)!=sz)return 2;fclose(f);
  }
  DNNResult r=dnn_bound(qr.A,qr.B,qr.C,qr.shift+shiftadd,HUGE_VAL,p,
                       start.Z.empty()?NULL:&start,&save);
  printf("RESULT m=%d lb=%.17g primal=%.17g status=%s iters=%d cpu=%.6f margin=%.9g\n",
    m,r.lb,r.primal,dnn_status_name(r.status),r.iterations,r.cpu,r.margin);
  FILE*f=fopen((pre+".Y.bin").c_str(),"wb");if(!f)return 2;
  fwrite(save.Y.data(),sizeof(double),sz,f);fclose(f);
  f=fopen((pre+".warm.bin").c_str(),"wb");if(!f)return 2;
  fwrite(&save.scale,sizeof(double),1,f);fwrite(save.Z.data(),sizeof(double),sz,f);fclose(f);
  delete q;return r.status==DNN_ERROR?1:0;
}
