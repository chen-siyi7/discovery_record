#include <Rcpp.h>
#include "../analysis/closure_impl.hpp"

Rcpp::LogicalMatrix closure_rcpp(Rcpp::NumericMatrix x,
  Rcpp::LogicalMatrix mandatory, double alpha=0.05, int mode=0) {
  int nr=x.nrow(),m=x.ncol();
  if(mandatory.nrow()!=nr || mandatory.ncol()!=m)Rcpp::stop("Dimension mismatch");
  std::vector<double> data(nr*m);
  std::vector<unsigned char> force(nr*m),out(nr*m);
  for(int r=0;r<nr;++r)for(int j=0;j<m;++j){
    double z=x(r,j);if(!std::isfinite(z)||z<0||(mode==1&&z>1))Rcpp::stop("Invalid closure input");
    data[r*m+j]=z;force[r*m+j]=mandatory(r,j);
  }
  int rc=closure_rows(data.data(),force.data(),out.data(),nr,m,alpha,mode);
  if(rc)Rcpp::stop("No certified superset at row %i",rc);
  Rcpp::LogicalMatrix ans(nr,m);
  for(int r=0;r<nr;++r)for(int j=0;j<m;++j)ans(r,j)=out[r*m+j];
  return ans;
}

double e_cap_rcpp(int m,double alpha=.05){return ledger::e_cap(m,alpha);}

Rcpp::NumericVector by_cutoffs_rcpp(int m,double alpha=.05) {
  auto h=ledger::harmonic(m);Rcpp::NumericVector ans(m);
  for(int b=1;b<=m;++b)ans[b-1]=ledger::by_cut_down(m,b,alpha,h[m]);
  return ans;
}

Rcpp::IntegerMatrix by_bins_rcpp(Rcpp::NumericMatrix p,double alpha=.05) {
  int nr=p.nrow(),m=p.ncol();auto h=ledger::harmonic(m);
  Rcpp::IntegerMatrix ans(nr,m);
  for(int r=0;r<nr;++r)for(int j=0;j<m;++j)ans(r,j)=ledger::by_bin(p(r,j),m,alpha,h[m]);
  return ans;
}

Rcpp::LogicalMatrix reject_e_rcpp(Rcpp::NumericMatrix e,double alpha=.05) {
  int nr=e.nrow(),m=e.ncol();Rcpp::LogicalMatrix ans(nr,m);std::vector<double>v(m);
  for(int r=0;r<nr;++r) {
    for(int j=0;j<m;++j)v[j]=e(r,j);
    std::sort(v.begin(),v.end(),std::greater<double>());int k=0;
    for(int j=0;j<m;++j)if(ledger::e_pass(v[j],j+1,m,alpha))k=j+1;
    if(k)for(int j=0;j<m;++j)ans(r,j)=ledger::e_pass(e(r,j),k,m,alpha);
  }
  return ans;
}

Rcpp::LogicalMatrix step_up_rcpp(Rcpp::NumericMatrix x,
  Rcpp::NumericVector thresholds, bool decreasing=true) {
  int nr=x.nrow(),m=x.ncol();if(thresholds.size()!=m)Rcpp::stop("Threshold mismatch");
  Rcpp::LogicalMatrix out(nr,m);std::vector<double> v(m);
  for(int r=0;r<nr;++r){
    for(int j=0;j<m;++j)v[j]=x(r,j);
    if(decreasing)std::sort(v.begin(),v.end(),std::greater<double>());
    else std::sort(v.begin(),v.end());
    int k=0;for(int j=0;j<m;++j)if(decreasing?v[j]>=thresholds[j]:v[j]<=thresholds[j])k=j+1;
    if(k)for(int j=0;j<m;++j)out(r,j)=decreasing?x(r,j)>=thresholds[k-1]:x(r,j)<=thresholds[k-1];
  }
  return out;
}

RCPP_MODULE(discovery_functions) {
  Rcpp::function("closure_native", &closure_rcpp);
  Rcpp::function("e_cap_native", &e_cap_rcpp);
  Rcpp::function("by_cutoffs_native", &by_cutoffs_rcpp);
  Rcpp::function("by_bins_native", &by_bins_rcpp);
  Rcpp::function("reject_e_native", &reject_e_rcpp);
  Rcpp::function("step_up_native", &step_up_rcpp);
}
