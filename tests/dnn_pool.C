// White-box objective/serialization checks: include the implementation once,
// then link without dnncut.o (see dnn_pool.mk). No test-only public solver API.
#include "../latest/dnncut.C"
#include <numeric>
#include <unistd.h>

static void check(bool ok, const char *why) { if (!ok) throw std::runtime_error(why); }
struct Matrix {
  MAT a{}; std::vector<double> v; std::vector<double*> rows;
  Matrix(int n): v(n*n), rows(n) { a.m=a.n=n; for(int i=0;i<n;i++) rows[i]=v.data()+i*n; a.me=rows.data(); }
};
static void maps(const std::vector<int>& assignment,int side,std::vector<int>& orig,std::vector<int>& other) {
  int n=assignment.size();std::vector<int> inv(n,-1),row,col;
  for(int i=0;i<n;i++)if(assignment[i]>=0)inv[assignment[i]]=i;
  for(int i=0;i<n;i++){if(assignment[i]<0)row.push_back(i);if(inv[i]<0)col.push_back(i);}
  orig=side?col:row;other=side?row:col;
}
static std::vector<double> objective(const std::vector<Region>&regions,const std::vector<Cut>&cuts,int m,int side,double&constant) {
  int N=m*m+1;std::vector<double>L(N*N,0);constant=0;Lift lift{m,side};
  for(const Cut&t:cuts){add_rcut(t,regions[t.s],L,lift,N,t.mu);constant-=t.mu*cut_rhs(t);}return L;
}
static double value(const std::vector<double>&L,double constant,const std::vector<int>&permutation) {
  int m=permutation.size(),N=m*m+1;std::vector<int>x{0};for(int i=0;i<m;i++)x.push_back(1+i*m+permutation[i]);
  double v=constant;for(int i:x)for(int j:x)v+=L[i*N+j];return v;
}
int main() {
 try {
  dnn_set_blas_threads(1);std::string file="/tmp/qap-pool-check-"+std::to_string((long long)getpid());
  int roundtrips=0,children=0,permutations=0,invalid=0,integrations=0;
  for(int side=0;side<2;side++) {
    DNNCutParams cp;cp.side=side;cp.pool_assignment=std::vector<int>(8,-1);cp.pool_assignment[1]=5;cp.pool_assignment[6]=2;
    std::vector<int>orig;maps(cp.pool_assignment,side,orig,cp.pool_other_orig);pool_maps(8,side,cp.pool_assignment,orig,cp.pool_other_orig);
    std::vector<Region>regions={{"a",{0,2,4}},{"b",{1,3}},{"c",{0,1,2,3,5}}};std::vector<Cut>cuts;
    for(int s=0;s<3;s++)for(int a=0;a<4;a++)for(int b=a+1;b<5;b++)for(int c=b+1;c<6;c++)for(int type=0;type<4;type++)
      cuts.push_back({type,a,b,c,s,(1+(a+b+c+type)%8)/7.0});
    write_pool(file,8,cp,orig,regions,cuts,7,true,-123.125);
    std::vector<Region>loaded;std::vector<Cut>restored;int kept,drop,sr;bool sq;double sv;
    read_pool(file,8,cp,orig,loaded,restored,kept,drop,sr,sq,sv);
    check(kept==(int)cuts.size()&&drop==0&&sr==7&&sq&&sv==-123.125,"roundtrip metadata/counts");
    double c0,c1;auto L=objective(regions,cuts,6,side,c0),R=objective(loaded,restored,6,side,c1);
    check(L==R&&c0==c1,"roundtrip Lagrangian coefficients/constant differ");roundtrips++;
    Matrix A(6),B(6),C(6);for(int i=0;i<6;i++)for(int j=0;j<6;j++){A.a.me[i][j]=abs(i-j);B.a.me[i][j]=(i+j)%5;C.a.me[i][j]=(i-j)*(i-j);}
    DNNParams par;par.eig_mode=DNN_EIG_DOUBLE;par.max_iter=50;par.tol=0;par.stall_window=0;
    par.L_add=L.data();for(double x:L)par.L_add_abs+=fabs(x);auto v0=dnn_bound(&A.a,&B.a,&C.a,c0,HUGE_VAL,par);
    par.L_add=R.data();auto v1=dnn_bound(&A.a,&B.a,&C.a,c1,HUGE_VAL,par);
    check(v0.status!=DNN_ERROR&&v0.lb==v1.lb&&v0.primal==v1.primal,"roundtrip DNN Lagrangian solve differs");
    std::vector<int>pfree_rows,pfree_cols;maps(cp.pool_assignment,0,pfree_rows,pfree_cols);
    for(int ri=0;ri<6;ri++)for(int cj=0;cj<6;cj++) {
      DNNCutParams child=cp;child.pool_assignment[pfree_rows[ri]]=pfree_cols[cj];std::vector<int>co;
      maps(child.pool_assignment,side,co,child.pool_other_orig);std::vector<Region>cr;std::vector<Cut>cc;
      read_pool(file,8,child,co,cr,cc,kept,drop,sr,sq,sv);
      int fixed_other=side?ri:cj,expected=0;for(const Cut&t:cuts)expected+=(t.a!=fixed_other&&t.b!=fixed_other&&t.c!=fixed_other);
      check(kept==expected&&kept+drop==(int)cuts.size(),"wrong-side child kept/dropped count");
      double cc0;auto CL=objective(cr,cc,5,side,cc0);std::vector<int>p(5);std::iota(p.begin(),p.end(),0);
      do {
        std::vector<int>pp(6);pp[ri]=cj;for(int i=0,ii=0;i<6;i++)if(i!=ri){int j=p[ii++];pp[i]=j+(j>=cj);}
        double full=value(L,c0,pp),inherited=value(CL,cc0,p);
        check(inherited>=full-1e-10,"dropped triangle increased negative Lagrangian penalty");
        // Untouched cuts retain exactly the same value under the child lift.
        std::vector<Cut>untouched;for(const Cut&t:cuts)if(t.a!=fixed_other&&t.b!=fixed_other&&t.c!=fixed_other)untouched.push_back(t);
        double uc;auto UL=objective(regions,untouched,6,side,uc);
        check(fabs(value(UL,uc,pp)-inherited)<1e-10,"child original-index remap changes retained objective");permutations++;
      } while(std::next_permutation(p.begin(),p.end()));children++;
    }
    // Corrupt valid symbolic records; no malformed multiplier/index may enter a certificate.
    Cut one={0,0,1,2,0,0.25};write_pool(file,8,cp,orig,regions,{one},1,false,-123.125);
    std::ifstream good_in(file.c_str());std::string good((std::istreambuf_iterator<char>(good_in)),{});
    const std::string indices=std::to_string(cp.pool_other_orig[0])+" "+std::to_string(cp.pool_other_orig[1])+" "+std::to_string(cp.pool_other_orig[2]);
    const std::string record="0 0 "+indices+" 0.25";
    for(const std::string &bad_record:std::vector<std::string>{"6 0 "+indices+" 0.25", "0 0 "+indices+" -1",
        "0 0 "+indices+" nan", "0 0 "+indices+" 1e308", "0 0 -1 1 2 0.25", "0 0 2 1 0 0.25", "0 999 "+indices+" 0.25"}) {
      std::string bad_text=good;size_t at=bad_text.find(record);check(at!=std::string::npos,"corruption fixture record");
      bad_text.replace(at,record.size(),bad_record);std::ofstream(file.c_str())<<bad_text;
      std::vector<Region>r;std::vector<Cut>c;bool rejected=false;
      try{read_pool(file,8,cp,orig,r,c,kept,drop,sr,sq,sv);}catch(const std::exception&){rejected=true;}
      check(rejected,"malformed symbolic record accepted");invalid++;
    }
    // The real library save and frozen-load path, with exact enumeration safety.
    DNNCutParams actual;actual.side=side;actual.prefix=true;actual.pool_assignment=cp.pool_assignment;
    actual.pool_other_orig=cp.pool_other_orig;actual.save_pool=file;actual.first_iters=20;actual.rounds=3;actual.iters=20;
    actual.kcuts=20;actual.freeze_iters=30;actual.delta=1;actual.cap_at_cutoff=0;
    par.L_add=nullptr;par.L_add_abs=0;double opt=HUGE_VAL;std::vector<int>p(6);std::iota(p.begin(),p.end(),0);
    do {double v=0;for(int i=0;i<6;i++){v+=C.a.me[i][p[i]];for(int j=0;j<6;j++)v+=A.a.me[i][j]*B.a.me[p[i]][p[j]];}opt=std::min(opt,v);}while(std::next_permutation(p.begin(),p.end()));
    DNNWarm w;auto source=dnn_cut_rounds(&A.a,&B.a,&C.a,0,opt,par,orig,8,actual,&w);
    check(source.status!=DNN_ERROR&&source.best_lb<=opt+1e-7,"saved source certificate invalid");
    actual.save_pool.clear();actual.load_pool=file;actual.rounds=0;actual.first_iters=0;actual.freeze_iters=25;
    auto frozen=dnn_cut_rounds(&A.a,&B.a,&C.a,0,opt,par,orig,8,actual,&w);
    check(frozen.status!=DNN_ERROR&&frozen.status!=DNN_FATHOMED&&frozen.best_lb<=opt+1e-7,"loaded frozen certificate invalid");
    check(frozen.iterations==25&&frozen.rounds==1&&std::isfinite(frozen.first_lb),"frozen mode must run exactly one solve");integrations++;
    // Failed ancestor and unsupported-family requests must not return a bound.
    actual.pool_assignment[1]=4;maps(actual.pool_assignment,side,orig,actual.pool_other_orig);
    auto bad=dnn_cut_rounds(&A.a,&B.a,&C.a,0,HUGE_VAL,par,orig,8,actual,&w);
    check(bad.status==DNN_ERROR&&!std::isfinite(bad.best_lb),"wrong ancestor accepted");invalid++;
    actual.pool_assignment=cp.pool_assignment;actual.pool_other_orig=cp.pool_other_orig;maps(cp.pool_assignment,side,orig,actual.pool_other_orig);
    actual.kcuts5=1;bad=dnn_cut_rounds(&A.a,&B.a,&C.a,0,HUGE_VAL,par,orig,8,actual,&w);
    check(bad.status==DNN_ERROR&&!std::isfinite(bad.best_lb),"unsupported family accepted");invalid++;
    std::ofstream(file.c_str())<<"QAP_REGIONAL_TRIANGLE_POOL 99\n";actual.kcuts5=0;
    bad=dnn_cut_rounds(&A.a,&B.a,&C.a,0,HUGE_VAL,par,orig,8,actual,&w);
    check(bad.status==DNN_ERROR&&!std::isfinite(bad.best_lb),"malformed file accepted");invalid++;
  }
  std::remove(file.c_str());printf("PASS roundtrips=%d child_maps=%d permutation_lifts=%d frozen_integrations=%d invalid_rejections=%d\n",roundtrips,children,permutations,integrations,invalid);
  return 0;
 }catch(const std::exception&e){fprintf(stderr,"FAIL %s\n",e.what());return 1;}
}
