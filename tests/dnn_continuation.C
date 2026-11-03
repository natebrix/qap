// Same-node continuation: identity tests, exact-enumeration safety, and paired probes.
#include "dnn.h"
#include "dnncut.h"
#include "qap.h"
#include "assign.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct Matrix {
  MAT a{}; std::vector<double> values; std::vector<double*> rows;
  explicit Matrix(int n): values(n*n), rows(n) {
    a.m=a.n=n; for(int i=0;i<n;i++) rows[i]=values.data()+i*n; a.me=rows.data();
  }
};
static void require(bool ok, const char* what) { if(!ok) throw std::runtime_error(what); }
static double distance(const std::vector<double>& a,const std::vector<double>& b) {
  require(a.size()==b.size(),"state dimensions"); double d=0;
  for(size_t i=0;i<a.size();i++) {require(std::isfinite(a[i])&&std::isfinite(b[i]),"finite state");d=std::max(d,fabs(a[i]-b[i]));}return d;
}
static DNNParams params(int mode,int scheme) {
  DNNParams p;p.eig_mode=mode;p.gamma=scheme?1:1.618;p.gamma1=scheme?1:0;
  p.cheb_degree=3;p.tight_lam=2;p.check_every=5;p.tol=0;p.stall_window=0;
  p.hopeless_after=0;p.inner_lap=DNN_LAP_EXACT;return p;
}
static double optimum(const MAT* A,const MAT* B,const MAT* C,double shift) {
  int n=A->m;std::vector<int> p(n);std::iota(p.begin(),p.end(),0);double best=HUGE_VAL;
  do {double v=shift;for(int i=0;i<n;i++) {v+=C->me[i][p[i]];for(int j=0;j<n;j++)v+=A->me[i][j]*B->me[p[i]][p[j]];}best=std::min(best,v);}while(std::next_permutation(p.begin(),p.end()));return best;
}
static void safe(const DNNResult&r,double opt,int&count) {
  require(r.status!=DNN_ERROR,"solver error");require(std::isfinite(r.lb),"finite certificate");
  require(r.lb<=opt+1e-7,"certificate above exact optimum");require(r.status!=DNN_FATHOMED,"false fathom at optimum");count++;
}
static std::vector<double> valid_cut(int n,double mu) {
  // -z_0 <= 0, z_0 = membership of facility 0 in locations {0,1}.
  int N=n*n+1;std::vector<double>L(N*N,0);for(int a=0;a<2;a++){int i=1+a;L[i*N+i]=-mu;}return L;
}
static void checks() {
  std::mt19937 rng(20261102);int cases=0,certs=0,negative=0,rescales=0,cutcases=0;
  double maxstate=0,maxlb=0,maxlegacy=0;
  for(int n=3;n<=8;n++) for(int seed=0;seed<3;seed++) {
    Matrix A(n),B(n),C(n);double shift=seed-3;
    for(int i=0;i<n;i++)for(int j=i;j<n;j++) {
      A.a.me[i][j]=A.a.me[j][i]=(int)(rng()%16)-5;
      B.a.me[i][j]=B.a.me[j][i]=(int)(rng()%13)-4;
    }
    for(auto&v:C.values)v=(int)(rng()%11)-5;
    double opt=optimum(&A.a,&B.a,&C.a,shift);
    for(int mode=0;mode<=1;mode++)for(int scheme=0;scheme<=1;scheme++) {
      DNNParams p=params(mode,scheme);DNNWarm whole,part,zonly;
      p.max_iter=45;auto all=dnn_bound(&A.a,&B.a,&C.a,shift,HUGE_VAL,p,nullptr,&whole);safe(all,opt,certs);
      double best=-HUGE_VAL;DNNResult r;
      for(int steps:{10,15,20}) {
        p.max_iter=steps;p.warm_primal=1;
        r=dnn_bound(&A.a,&B.a,&C.a,shift,HUGE_VAL,p,part.Z.empty()?nullptr:&part,&part);
        safe(r,opt,certs);best=std::max(best,r.lb);
      }
      double dy=distance(whole.Y,part.Y),dz=distance(whole.Z,part.Z);
      maxstate=std::max(maxstate,std::max(dy,dz));maxlb=std::max(maxlb,fabs(all.lb-best));
      require(dy<=1e-12&&dz<=1e-12,"split Y/Z differs from uninterrupted full-spectrum state");
      require(fabs(all.lb-best)<=1e-7,"split running certificate differs");
      // Negative control: same chunks, legacy Z-only restarts.
      p.warm_primal=0;
      for(int steps:{10,15,20}){p.max_iter=steps;r=dnn_bound(&A.a,&B.a,&C.a,shift,HUGE_VAL,p,zonly.Z.empty()?nullptr:&zonly,&zonly);safe(r,opt,certs);}
      double d=std::max(distance(whole.Y,zonly.Y),distance(whole.Z,zonly.Z));maxlegacy=std::max(maxlegacy,d);negative+=d>1e-6;cases++;
      // A changed cut objective must preserve Y and rescale Z by old/new scale.
      auto L=valid_cut(n,1024);p.L_add=L.data();p.L_add_abs=2048;p.warm_primal=1;p.max_iter=5;
      DNNWarm scaled=whole,out1,out2;auto q=dnn_bound(&A.a,&B.a,&C.a,shift,opt,p,&whole,&out1);safe(q,opt,certs);
      require(out1.scale!=whole.scale,"rescaling fixture did not change scale");
      for(auto&v:scaled.Z)v*=whole.scale/out1.scale;scaled.scale=out1.scale;
      auto q2=dnn_bound(&A.a,&B.a,&C.a,shift,opt,p,&scaled,&out2);safe(q2,opt,certs);
      require(distance(out1.Y,out2.Y)<=1e-12&&distance(out1.Z,out2.Z)<=1e-12,"objective rescaling mismatch");rescales++;
    }
    // Real changing-cut rounds, including best-round quiet restoration, both axes.
    for(int mode=0;mode<=2;mode++)for(int side=0;side<=1;side++) {
      auto p=params(mode,1);DNNCutParams cp;cp.warm_primal=1;cp.prefix=true;cp.side=side;
      cp.first_iters=20;cp.iters=20;cp.rounds=3;cp.kcuts=10;cp.delta=1;
      cp.cap_at_cutoff=0;cp.freeze_iters=30;cp.freeze_best=1;
      std::vector<int> free(n);std::iota(free.begin(),free.end(),0);DNNWarm w;
      auto cr=dnn_cut_rounds(&A.a,&B.a,&C.a,shift,opt,p,free,n,cp,&w);
      require(cr.status!=DNN_ERROR&&cr.status!=DNN_FATHOMED,"cut continuation false fathom/error");
      require(std::isfinite(cr.best_lb)&&cr.best_lb<=opt+1e-7,"cut certificate above optimum");cutcases++;
    }
  }
  require(negative>0,"negative control not sensitive");
  // Malformed explicitly requested same-node starts must fail closed.
  Matrix A(3),B(3),C(3);auto p=params(0,1);p.max_iter=5;p.warm_primal=1;DNNWarm w;
  auto ok=dnn_bound(&A.a,&B.a,&C.a,0,HUGE_VAL,p,nullptr,&w);require(ok.status!=DNN_ERROR,"zero fixture");
  for(int k=0;k<6;k++) {
    DNNWarm bad=w;if(k==0)bad.Y.pop_back();if(k==1)bad.Z.pop_back();if(k==2)bad.scale=0;
    if(k==3)bad.Y[0]=std::numeric_limits<double>::quiet_NaN();if(k==4)bad.m=4;
    if(k==5){bad.scale=1e308;bad.Z[0]=1e308;}
    auto r=dnn_bound(&A.a,&B.a,&C.a,0,HUGE_VAL,p,&bad);require(r.status==DNN_ERROR,"bad start accepted");
  }
  // Child transfer remains legacy Z-only: opting into same-node Y must reject it.
  DNNWarm child;require(dnn_warm_child(w,0,0,child),"child mapping");Matrix D(2),E(2),F(2);
  auto bad=dnn_bound(&D.a,&E.a,&F.a,0,HUGE_VAL,p,&child);require(bad.status==DNN_ERROR,"child warm accepted as same-node Y");
  printf("PASS identity_cases=%d safe_certificates=%d rescaling_cases=%d cut_pool_cases=%d malformed_cases=7 max_state_error=%.17g max_certificate_error=%.17g legacy_different=%d legacy_max_state_error=%.17g\n",cases,certs,rescales,cutcases,maxstate,maxlb,negative,maxlegacy);
}
static bool fixes(const char* spec,QAPAssignment&node,int n) {
  std::string s=spec;if(s=="-")return true;
  for(size_t p=0;p<s.size();){size_t e=s.find(',',p);if(e==std::string::npos)e=s.size();int i,j;
    if(sscanf(s.substr(p,e-p).c_str(),"%d:%d",&i,&j)!=2||i<0||j<0||i>=n||j>=n||node.p[i]!=NONE||node.pinv[j]!=NONE)return false;
    node.fix(i,j);p=e+1;}return true;
}
static void print(const char*arm,const DNNResult&r,const DNNWarm&w) {
  printf("ARM %s lb=%.17g primal=%.17g iterations=%d cpu=%.6f scale=%.17g status=%s\n",arm,r.lb,r.primal,r.iterations,r.cpu,w.scale,dnn_status_name(r.status));fflush(stdout);
}
int main(int argc,char**argv) {
 try {
  setvbuf(stdout,nullptr,_IOLBF,0);dnn_set_blas_threads(1);
  if(argc==1){checks();return 0;}
  require(argc==8||argc==10,"usage: dnncontinuation --probe instance fixes first rest eig scheme [L.bin shift_add]");
  require(std::string(argv[1])=="--probe","bad mode");
  QAP*q=readQAPLIB(argv[2]);require(q!=nullptr,"read instance");symmetrize_QAP(q);
  int n=q->A->m;QAPAssignment node(n);require(fixes(argv[3],node,n),"bad fixes");QAP qr;reduceQAP(q,&qr,node);
  int first=atoi(argv[4]),rest=atoi(argv[5]),mode=atoi(argv[6]),scheme=atoi(argv[7]);
  require(first>0&&rest>0&&first%5==0&&rest%5==0,"positive check-aligned iterations required");
  DNNParams p=params(mode,scheme);std::vector<double>L;
  if(argc==10){size_t N=(n-node.nfix)*(n-node.nfix)+1;L.resize(N*N);std::ifstream f(argv[8],std::ios::binary);f.read((char*)L.data(),L.size()*sizeof(double));require((bool)f,"L file");p.L_add=L.data();for(auto v:L)p.L_add_abs+=fabs(v);qr.shift+=atof(argv[9]);}
  printf("PROBE instance=%s fixes=%s first=%d rest=%d eig=%d scheme=%d\n",argv[2],argv[3],first,rest,mode,scheme);
  DNNWarm initial,z,y,all;p.max_iter=first;auto r=dnn_bound(qr.A,qr.B,qr.C,qr.shift,HUGE_VAL,p,nullptr,&initial);print("first",r,initial);double lb0=r.lb;
  p.max_iter=rest;r=dnn_bound(qr.A,qr.B,qr.C,qr.shift,HUGE_VAL,p,&initial,&z);print("Z_only",r,z);double lbz=std::max(lb0,r.lb);
  p.warm_primal=1;r=dnn_bound(qr.A,qr.B,qr.C,qr.shift,HUGE_VAL,p,&initial,&y);print("YZ",r,y);double lby=std::max(lb0,r.lb);
  p.max_iter=first+rest;r=dnn_bound(qr.A,qr.B,qr.C,qr.shift,HUGE_VAL,p,nullptr,&all);print("uninterrupted",r,all);
  printf("COMPARE best_Z_only=%.17g best_YZ=%.17g YZ_minus_Z=%.17g Y_error=%.17g Z_error=%.17g\n",lbz,lby,lby-lbz,distance(y.Y,all.Y),distance(y.Z,all.Z));
  delete q;return 0;
 } catch(const std::exception&e){fprintf(stderr,"FAIL %s\n",e.what());return 1;}
}
