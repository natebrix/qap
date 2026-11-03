#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include "dnn.h"
struct Matrix{MAT mat{};std::vector<double>v;std::vector<double*>rows;Matrix(int n):v(n*n),rows(n){mat.m=mat.n=n;for(int i=0;i<n;i++)rows[i]=v.data()+i*n;mat.me=rows.data();}};
void report(const char*label,DNNResult r){printf("%s lb=%.12f primal=%.9f iter=%d cpu=%.5f status=%s\n",label,r.lb,r.primal,r.iterations,r.cpu,dnn_status_name(r.status));fflush(stdout);}
int main(int argc,char**argv){std::ifstream in(argv[1]);int n;in>>n;Matrix A(n),B(n),C(n);for(auto&x:A.v)in>>x;for(auto&x:B.v)in>>x;dnn_set_blas_threads(1);DNNParams p;p.hopeless_after=0;p.stall_window=0;p.tol=0;p.gamma=p.gamma1=1;p.cheb_degree=3;p.tight_lam=2;p.eig_mode=argc>2?atoi(argv[2]):2;p.verbose=0;p.max_iter=200;DNNWarm warm;report("first200",dnn_bound(&A.mat,&B.mat,&C.mat,0,HUGE_VAL,p,nullptr,&warm));p.max_iter=100;unsetenv("SPEED_USE_Y");report("Zonly+100",dnn_bound(&A.mat,&B.mat,&C.mat,0,HUGE_VAL,p,&warm));setenv("SPEED_USE_Y","1",1);report("YZ+100",dnn_bound(&A.mat,&B.mat,&C.mat,0,HUGE_VAL,p,&warm));p.max_iter=300;report("continuous300",dnn_bound(&A.mat,&B.mat,&C.mat,0,HUGE_VAL,p));}
