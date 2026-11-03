#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
static double cpu(){timespec t;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);return t.tv_sec+t.tv_nsec*1e-9;}
inline double dot4(const double*x,const double*y,size_t n){double a0=0,a1=0,a2=0,a3=0;size_t i=0;for(;i+4<=n;i+=4){a0+=x[i]*y[i];a1+=x[i+1]*y[i+1];a2+=x[i+2]*y[i+2];a3+=x[i+3]*y[i+3];}for(;i<n;i++)a0+=x[i]*y[i];return(a0+a1)+(a2+a3);}
struct Data{int m=47,n1=2210;std::vector<double>L,Y,Z,V,Zc,tmp;Data():L(n1*n1),Y(n1*n1),Z(n1*n1),V(n1*n1),Zc(n1*n1),tmp(4*n1){for(size_t t=0;t<L.size();t++){L[t]=(int(t%131)-65)*.01;Y[t]=(t%127)*.001;Z[t]=(int(t%113)-56)*.001;V[t]=(t%117)*.001;}}};
template<bool SKIP> double step(Data&w,int it){int n1=w.n1,m=w.m;double beta=4.7,ibeta=1/beta,iscale=.01;bool check=it%10==0;long double pres2=0,dy2=0;double *yn=w.tmp.data();
for(int r=0;r<n1;r++){const double *lr=w.L.data()+size_t(r)*n1,*vr=w.V.data()+size_t(r)*n1;double*yr=w.Y.data()+size_t(r)*n1,*zr=w.Z.data()+size_t(r)*n1,*cr=w.Zc.data()+size_t(r)*n1;
if(check)for(int c=0;c<n1;c++)cr[c]=zr[c]+beta*(yr[c]-vr[c]);
for(int c=0;c<n1;c++)zr[c]+=beta*(yr[c]-vr[c]);
for(int c=0;c<n1;c++){double y=vr[c]-(lr[c]*iscale+zr[c])*ibeta;yn[c]=std::min(1.0,std::max(0.0,y));}
if(r==0)yn[0]=1;else{int i=(r-1)/m,j=(r-1)%m;for(int l=0;l<m;l++)if(l!=j)yn[1+i*m+l]=0;for(int k=0;k<m;k++)if(k!=i)yn[1+k*m+j]=0;}
double *rs=yn+n1,*ds=rs+n1;
if(!SKIP||check){for(int c=0;c<n1;c++){double y=yn[c],res_=y-vr[c];zr[c]+=beta*res_;rs[c]=res_;ds[c]=y-yr[c];yr[c]=y;}pres2+=dot4(rs,rs,n1);dy2+=dot4(ds,ds,n1);}
else for(int c=0;c<n1;c++){double y=yn[c],res_=y-vr[c];zr[c]+=beta*res_;yr[c]=y;}
}return double(pres2+dy2);}
int main(){Data a,b;for(int trial=0;trial<6;trial++){double t=cpu(),ca=0,cb=0;for(int i=1;i<=30;i++){double r=step<false>(a,i);if(i%10==0)ca+=r;}double ta=cpu()-t;t=cpu();for(int i=1;i<=30;i++){double r=step<true>(b,i);if(i%10==0)cb+=r;}double tb=cpu()-t;bool equal=!memcmp(a.Y.data(),b.Y.data(),a.Y.size()*8)&&!memcmp(a.Z.data(),b.Z.data(),a.Z.size()*8)&&!memcmp(a.Zc.data(),b.Zc.data(),a.Zc.size()*8)&&ca==cb;printf("trial%d old_ms=%.3f skip_ms=%.3f equal=%d\n",trial,ta/30*1e3,tb/30*1e3,int(equal));}}
