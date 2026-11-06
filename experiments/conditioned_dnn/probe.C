// Complete assignment disjunctions and identity-checked child Y/Z reuse.
// Research driver only: the production solver and certificate are unchanged.
#include "dnn.h"
#include "qap.h"
#include "assign.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static void need(bool b, const char *s) { if (!b) throw std::runtime_error(s); }
static double seconds() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static int integer(const char *s) {
  size_t k; int v=std::stoi(s,&k); need(k==strlen(s),"invalid integer"); return v;
}
static DNNParams params(int it, int eig) {
  need(it>0 && it%10==0,"iterations must be positive multiples of 10");
  need(eig>=0 && eig<=2,"invalid eigensolver");
  DNNParams p; p.max_iter=it; p.eig_mode=eig; p.gamma=p.gamma1=1;
  p.beta=.3; p.refresh=40; p.cheb_degree=3; p.tight_lam=2;
  p.check_every=10; p.inner_lap=DNN_LAP_EXACT;
  p.tol=0; p.hopeless_after=0; p.stall_window=0; return p;
}
static void fixes(const char *arg,QAPAssignment &a) {
  std::string s=arg; if(s=="-")return;
  for(size_t k=0;k<s.size();) {
    size_t e=s.find(',',k); if(e==std::string::npos)e=s.size();
    int i,j,used=0; std::string t=s.substr(k,e-k);
    need(sscanf(t.c_str(),"%d:%d%n",&i,&j,&used)==2 && used==(int)t.size(),"invalid fixing");
    need(i>=0&&i<a.n&&j>=0&&j<a.n&&a.p[i]==NONE&&a.pinv[j]==NONE,"conflicting fixing");
    a.fix(i,j); k=e+1;
  }
}
// Explicit original-index labels plus all reduced data bind a state to one node.
static std::vector<double> identity(const QAP &q,const QAPAssignment &a) {
  std::vector<double> v{double(a.n),q.shift};
  for(int i=0;i<a.n;i++)v.push_back(a.p[i]);
  for(const MAT *m:{q.A,q.B,q.C})
    for(int i=0;i<(int)m->m;i++)for(int j=0;j<(int)m->n;j++)v.push_back(m->me[i][j]);
  return v;
}
struct Snapshot {
  DNNWarm warm; double lb=-HUGE_VAL; int eig=0;
  std::vector<double> id,child_lb;
};
template<class T> static void put(std::ostream &f,const T &v) {f.write((const char*)&v,sizeof(v));}
template<class T> static void get(std::istream &f,T &v) {f.read((char*)&v,sizeof(v));need(bool(f),"truncated snapshot");}
static void putv(std::ostream &f,const std::vector<double>&v) {
  uint64_t n=v.size();put(f,n);f.write((const char*)v.data(),n*sizeof(double));
}
static void getv(std::istream &f,std::vector<double>&v,uint64_t limit) {
  uint64_t n;get(f,n);need(n<=limit,"invalid snapshot size");v.resize(n);
  f.read((char*)v.data(),n*sizeof(double));need(bool(f),"truncated snapshot vector");
  for(double d:v)need(std::isfinite(d),"nonfinite snapshot vector");
}
static void save(const std::string &path,const Snapshot &s) {
  if(path=="-")return;
  std::string tmp=path+".tmp"; std::ofstream f(tmp,std::ios::binary);
  need(bool(f),"cannot open snapshot"); uint64_t magic=0x43444e4e00000001ULL;
  put(f,magic);put(f,s.warm.m);put(f,s.warm.scale);put(f,s.lb);put(f,s.eig);
  putv(f,s.id);putv(f,s.child_lb);putv(f,s.warm.Y);putv(f,s.warm.Z);
  f.close();need(bool(f),"snapshot write failed");need(rename(tmp.c_str(),path.c_str())==0,"snapshot rename failed");
}
static Snapshot load(const std::string &path,const std::vector<double>&id,int eig) {
  std::ifstream f(path,std::ios::binary);need(bool(f),"cannot read snapshot");
  uint64_t magic;get(f,magic);need(magic==0x43444e4e00000001ULL,"snapshot version");
  Snapshot s;get(f,s.warm.m);get(f,s.warm.scale);get(f,s.lb);get(f,s.eig);
  need(s.warm.m>=2&&s.warm.m<=100,"snapshot dimension");
  need(std::isfinite(s.lb)&&std::isfinite(s.warm.scale)&&s.warm.scale>0,"snapshot scalar");
  getv(f,s.id,id.size());need(s.id==id,"snapshot node/objective mismatch");
  need(s.eig==eig,"snapshot solver configuration mismatch");
  size_t N=s.warm.m*s.warm.m+1;
  getv(f,s.child_lb,s.warm.m*s.warm.m);getv(f,s.warm.Y,N*N);getv(f,s.warm.Z,N*N);
  need(s.warm.Y.size()==N*N&&s.warm.Z.size()==N*N,"snapshot state shape");
  need(f.peek()==EOF,"snapshot trailing bytes");return s;
}
static void valid(const DNNResult &r) {
  need(r.status!=DNN_ERROR&&std::isfinite(r.lb),"DNN failed to return finite certified bound");
}
static void emit(const std::string &mode,int i,int j,const DNNResult&r,const Snapshot&s,
                 double inherited,double elapsed,double io) {
  std::cout<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"facility\":"<<i
    <<",\"location\":"<<j<<",\"m\":"<<s.warm.m<<",\"lb\":"<<s.lb
    <<",\"raw_lb\":"<<r.lb<<",\"inherited_lb\":"<<inherited
    <<",\"cpu\":"<<r.cpu<<",\"wall\":"<<elapsed<<",\"save_wall\":"<<io
    <<",\"iterations\":"<<r.iterations<<",\"certificates\":"<<r.certificates
    <<",\"cert_failures\":"<<r.cert_failures<<",\"refinements\":"<<r.refinements
    <<",\"margin\":"<<r.margin<<",\"status\":\""<<dnn_status_name(r.status)
    <<"\",\"state_bytes\":"<<(s.warm.Y.size()+s.warm.Z.size())*sizeof(double)<<"}"<<std::endl;
}
static double optimum(const QAP&q) {
  int n=q.A->m;std::vector<int> p(n);std::iota(p.begin(),p.end(),0);double best=HUGE_VAL;
  do {double v=q.shift;for(int i=0;i<n;i++) {v+=q.C->me[i][p[i]];
    for(int j=0;j<n;j++)v+=q.A->me[i][j]*q.B->me[p[i]][p[j]];}
    best=std::min(best,v);
  } while(std::next_permutation(p.begin(),p.end()));return best;
}
// reduceQAP expects the original zero-shift instance, not an already reduced
// QAP. Preserve a nonzero source shift explicitly for the synthetic fixtures.
static void reduce(QAP *q,QAP *r,const QAPAssignment &a) {
  reduceQAP(q,r,a);LowerSum s;s.add(r->shift);s.add(q->shift);r->shift=s.lower();
}
static double difference(const std::vector<double>&a,const std::vector<double>&b) {
  need(a.size()==b.size(),"state shape mismatch");double d=0;
  for(size_t i=0;i<a.size();i++)d=std::max(d,fabs(a[i]-b[i]));return d;
}
static void selftest(const std::string &path) {
  std::mt19937 rng(20261103);int count=0,rejections=0;double err=0,lb_err=0;
  for(int n=3;n<=7;n++)for(int seed=0;seed<2;seed++) {
    QAP q(n);q.shift=seed-4;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) {
      q.A->me[i][j]=int(rng()%13)-4;q.B->me[i][j]=int(rng()%11)-3;
      q.C->me[i][j]=int(rng()%9)-5;
    }
    double opt=optimum(q);QAPAssignment node(n);
    auto p=params(20,0);p.child_bounds=1;Snapshot parent;parent.id=identity(q,node);
    auto pr=dnn_bound(q.A,q.B,q.C,q.shift,HUGE_VAL,p,nullptr,&parent.warm);valid(pr);
    parent.lb=pr.lb;parent.child_lb=pr.child_lb;need(pr.lb<=opt+1e-7,"unsafe parent");
    need(parent.child_lb.size()==size_t(n*n),"missing child bounds");
    int i=n/2;double minopt=HUGE_VAL,minprobe=HUGE_VAL,minresume=HUGE_VAL;
    for(int j=0;j<n;j++) {
      QAPAssignment child=node;child.fix(i,j);QAP qc;reduce(&q,&qc,child);
      double co=optimum(qc);minopt=std::min(minopt,co);
      double inherited=std::max(parent.lb,parent.child_lb[i*n+j]);
      need(inherited<=co+1e-7,"unsafe inherited bound");
      DNNWarm start;need(dnn_warm_child(parent.warm,i,j,start),"child mapping");
      p=params(20,0);Snapshot probe;probe.id=identity(qc,child);
      auto a=dnn_bound(qc.A,qc.B,qc.C,qc.shift,HUGE_VAL,p,&start,&probe.warm);valid(a);
      probe.lb=std::max(inherited,a.lb);save(path,probe);
      Snapshot disk=load(path,probe.id,0);
      need(difference(disk.warm.Y,probe.warm.Y)==0&&difference(disk.warm.Z,probe.warm.Z)==0,"snapshot roundtrip");
      p.warm_primal=1;DNNWarm resumed,full;
      auto b=dnn_bound(qc.A,qc.B,qc.C,qc.shift,HUGE_VAL,p,&disk.warm,&resumed);valid(b);
      p.warm_primal=0;p.max_iter=40;
      auto c=dnn_bound(qc.A,qc.B,qc.C,qc.shift,HUGE_VAL,p,&start,&full);valid(c);
      double e=std::max(difference(resumed.Y,full.Y),difference(resumed.Z,full.Z));
      err=std::max(err,e);need(e<=1e-12,"full-spectrum continuation state mismatch");
      double blb=std::max(probe.lb,b.lb),clb=std::max(inherited,c.lb);
      lb_err=std::max(lb_err,fabs(blb-clb));need(fabs(blb-clb)<=1e-7,"continuation certificate mismatch");
      need(probe.lb<=co+1e-7&&blb<=co+1e-7&&clb<=co+1e-7,"unsafe child bound");
      minprobe=std::min(minprobe,probe.lb);minresume=std::min(minresume,blb);count++;
      auto wrong=probe.id;wrong[1]+=1;
      try{load(path,wrong,0);need(false,"wrong-node snapshot accepted");}
      catch(const std::runtime_error&e){need(std::string(e.what())=="snapshot node/objective mismatch","unexpected error");rejections++;}
    }
    need(minopt==opt,"assignment split is not exhaustive");
    need(minprobe>=parent.lb&&minprobe<=opt+1e-7&&minresume>=minprobe&&minresume<=opt+1e-7,"unsafe disjunctive bound");
  }
  remove(path.c_str());std::cout<<"PASS children="<<count<<" exact_families=10 wrong_node_rejections="<<rejections
    <<" max_state_error="<<std::setprecision(17)<<err<<" max_lb_error="<<lb_err<<std::endl;
}
int main(int argc,char **argv) {
 try {
  setvbuf(stdout,nullptr,_IOLBF,0);need(dnn_set_blas_threads(1),"cannot enforce one BLAS thread");
  if(argc==3&&std::string(argv[1])=="selftest"){selftest(argv[2]);return 0;}
  need(argc==8||argc==11,"usage: probe parent INSTANCE FIXES ITERS EIG CUTOFF STATE | probe {probe,resume,fresh} INSTANCE FIXES FAC LOC INPUT_STATE ITERS EIG CUTOFF OUTPUT_STATE");
  double t0=seconds();std::string mode=argv[1];bool parent_mode=mode=="parent"||mode=="control";
  need((parent_mode&&argc==8)||(!parent_mode&&argc==11),"wrong argument count");
  need(parent_mode||mode=="probe"||mode=="resume"||mode=="fresh","invalid mode");
  std::unique_ptr<QAP> q(readQAPLIB(argv[2]));need(bool(q),"cannot read instance");symmetrize_QAP(q.get());
  int n=q->A->m;QAPAssignment node(n);fixes(argv[3],node);QAP reduced;reduce(q.get(),&reduced,node);
  int it=integer(argv[parent_mode?4:7]),eig=integer(argv[parent_mode?5:8]);
  double cutoff=std::stod(argv[parent_mode?6:9]);need(std::isfinite(cutoff),"finite cutoff required");
  auto p=params(it,eig);if(mode=="control"){need(cutoff>0,"positive CPU budget required");p.time_limit=cutoff;}
  Snapshot out;out.eig=eig;DNNResult r;double inherited=0;int i=-1,j=-1;
  if(parent_mode) {
    p.child_bounds=1;out.id=identity(reduced,node);
    r=dnn_bound(reduced.A,reduced.B,reduced.C,reduced.shift,HUGE_VAL,p,nullptr,&out.warm);
    valid(r);out.lb=r.lb;out.child_lb=r.child_lb;inherited=r.lb;
    need(out.child_lb.size()==size_t(out.warm.m*out.warm.m),"missing parent child bounds");
  } else {
    i=integer(argv[4]);j=integer(argv[5]);need(i>=0&&i<n&&j>=0&&j<n&&node.p[i]==NONE&&node.pinv[j]==NONE,"invalid child");
    int li=0,lj=0;for(int k=0;k<i;k++)li+=node.p[k]==NONE;for(int k=0;k<j;k++)lj+=node.pinv[k]==NONE;
    QAPAssignment child=node;child.fix(i,j);QAP qr;reduce(q.get(),&qr,child);out.id=identity(qr,child);
    Snapshot input=load(argv[6],mode=="resume"?out.id:identity(reduced,node),eig);
    DNNWarm start;
    if(mode=="resume") {inherited=input.lb;start=std::move(input.warm);p.warm_primal=1;}
    else {
      need(input.warm.m==(int)reduced.A->m&&input.child_lb.size()==size_t(input.warm.m*input.warm.m),"parent state mismatch");
      inherited=std::max(input.lb,input.child_lb[li*input.warm.m+lj]);
      need(dnn_warm_child(input.warm,li,lj,start),"child mapping failed");
    }
    r=dnn_bound(qr.A,qr.B,qr.C,qr.shift,cutoff,p,&start,&out.warm);valid(r);
    out.lb=std::max(inherited,r.lb);
  }
  double ts=seconds();save(argv[parent_mode?7:10],out);double io=seconds()-ts;
  emit(mode,i,j,r,out,inherited,seconds()-t0,io);return 0;
 } catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<std::endl;return 1;}
}
