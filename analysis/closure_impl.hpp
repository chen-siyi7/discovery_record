#include "closure_numerics.h"
#include <numeric>
using std::vector;
using ledger::cpp_rational;

static bool feasible(const vector<double>& raw,const vector<unsigned char>& mask,
                     double alpha,int mode,const vector<double>& h) {
  int m=raw.size(),k=std::accumulate(mask.begin(),mask.end(),0);
  if(!k)return true;
  vector<double>s,u;
  for(int j=0;j<m;++j)(mask[j]?s:u).push_back(raw[j]);
  int q=u.size();
  if(mode!=1) {
    auto normalized=[&](double z){return mode==2?(z?m/z:0.):std::min(alpha*z,double(m));};
    auto exact=[&](double z)->cpp_rational {
      if(mode==2)return z?cpp_rational(m)/int(z):cpp_rational(0);
      cpp_rational x=cpp_rational(alpha)*cpp_rational(z);
      return x>m?cpp_rational(m):x;
    };
    vector<double>cs(k+1,0),cu(q+1,0);
    for(int j=0;j<k;++j)cs[j+1]=cs[j]+normalized(s[j]);
    for(int j=0;j<q;++j)cu[j+1]=cu[j]+normalized(u[j]);
    vector<cpp_rational>xs,xu;
    auto check=[&](int a) {
      int lo=0,hi=q;
      while(lo<hi) {
        int mid=(lo+hi)/2;double v=normalized(u[mid])*k;
        double err=ledger::roundoff(3,std::abs(v)+a);
        bool below=v-a< -err?true:(v-a>err?false:exact(u[mid])*k<a);
        if(below)lo=mid+1;else hi=mid;
      }
      int b=lo;double lhs=k*(cs[a]+cu[b]),rhs=double(a)*(a+b);
      double err=ledger::roundoff(a+b+3,std::abs(lhs)+rhs);
      if(lhs-rhs>err)return true;
      if(lhs-rhs< -err)return false;
      if(xs.empty()) {
        xs.resize(k+1);xu.resize(q+1);
        for(int j=0;j<k;++j)xs[j+1]=xs[j]+exact(s[j]);
        for(int j=0;j<q;++j)xu[j+1]=xu[j]+exact(u[j]);
      }
      return k*(xs[a]+xu[b])>=cpp_rational(a)*(a+b);
    };
    if(!check((k+1)/2))return false;
    for(int a=1;a<=k;++a)if(!check(a))return false;
  } else {
    for(int n=1;n<=m;++n) {
      vector<int>bs(k),bu(q);vector<double>cs(k+1,0),cu(q+1,0);
      for(int j=0;j<k;++j){bs[j]=ledger::by_bin(s[j],n,alpha,h[n]);cs[j+1]=cs[j]+(bs[j]?1./bs[j]:0.);}
      for(int j=0;j<q;++j){bu[j]=ledger::by_bin(u[j],n,alpha,h[n]);cu[j+1]=cu[j]+(bu[j]?1./bu[j]:0.);}
      vector<cpp_rational>xs,xu;
      for(int a=std::max(1,n-q);a<=std::min(k,n);++a) {
        int b=n-a;double lhs=k*(cs[a]+cu[b]);
        double err=ledger::roundoff(n+3,std::abs(lhs)+a);
        if(lhs-a< -err)return false;
        if(lhs-a<=err) {
          if(xs.empty()) {
            xs.resize(k+1);xu.resize(q+1);
            for(int j=0;j<k;++j)xs[j+1]=xs[j]+(bs[j]?cpp_rational(1)/bs[j]:cpp_rational(0));
            for(int j=0;j<q;++j)xu[j+1]=xu[j]+(bu[j]?cpp_rational(1)/bu[j]:cpp_rational(0));
          }
          if(k*(xs[a]+xu[b])<a)return false;
        }
      }
    }
  }
  return true;
}
extern "C" int closure_rows(const double* input,const unsigned char* force,
                            unsigned char* output,int rows,int m,double alpha,int mode) {
  auto h=ledger::harmonic(m);
  for(int row=0;row<rows;++row) {
    vector<int>order(m);std::iota(order.begin(),order.end(),0);
    auto val=[&](int j){return input[row*m+j];};
    auto strength=[&](int j){return mode==2?(val(j)?m/val(j):0.):val(j);};
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return mode==1?val(a)>val(b):strength(a)<strength(b);});
    vector<double>v(m);vector<unsigned char>mandatory(m),mask(m);
    for(int j=0;j<m;++j){v[j]=val(order[j]);mandatory[j]=force[row*m+order[j]];}
    int min_k=std::accumulate(mandatory.begin(),mandatory.end(),0);bool found=false;
    for(int k=m;k>=min_k;--k) {
      mask=mandatory;int extra=k-min_k;
      for(int j=m-1;j>=0&&extra;--j)if(!mask[j]){mask[j]=1;--extra;}
      if(feasible(v,mask,alpha,mode,h)){found=true;break;}
    }
    if(!found)return row+1;
    for(int j=0;j<m;++j)output[row*m+order[j]]=mask[j];
  }
  return 0;
}
