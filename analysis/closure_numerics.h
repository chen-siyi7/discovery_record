#ifndef CHEN_CLOSURE_NUMERICS_H
#define CHEN_CLOSURE_NUMERICS_H
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace ledger {
using boost::multiprecision::cpp_rational;
inline std::vector<double> harmonic(int m) {
  std::vector<double> h(m+1,0.);
  for(int j=1;j<=m;++j) h[j]=h[j-1]+1./j;
  return h;
}

inline double by_cut(int n,int b,double alpha,double h){return (alpha*b)/(n*h);}
inline const cpp_rational& exact_harmonic(int n) {
  static std::vector<cpp_rational> hs(1,cpp_rational(0));
  while(int(hs.size())<=n){int j=hs.size();hs.push_back(hs.back()+cpp_rational(1)/j);}
  return hs[n];
}
inline bool by_leq(double p,int n,int b,double alpha,double h) {
  double lhs=(p*n)*h,rhs=alpha*b;
  double u=std::numeric_limits<double>::epsilon()/2,nu=(n+16.)*u;
  double err=(nu/(1-nu))*std::max(1.,std::abs(lhs)+std::abs(rhs));
  if(lhs-rhs< -err)return true;
  if(lhs-rhs>err)return false;
  return cpp_rational(p)*n*exact_harmonic(n)<=cpp_rational(alpha)*b;
}
inline double by_cut_down(int n,int b,double alpha,double h) {
  cpp_rational cut=cpp_rational(alpha)*b/(n*exact_harmonic(n));
  double x=by_cut(n,b,alpha,h);
  while(cpp_rational(x)>cut)x=std::nextafter(x,0.);
  while(cpp_rational(std::nextafter(x,std::numeric_limits<double>::infinity()))<=cut)
    x=std::nextafter(x,std::numeric_limits<double>::infinity());
  return x;
}
inline int by_bin(double p,int n,double alpha,double h) {
  if(!by_leq(p,n,n,alpha,h))return 0;
  double hint=(n*h*p)/alpha;
  int b=hint>=n?n:std::max(1,int(std::ceil(hint)));
  while(b>1 && by_leq(p,n,b-1,alpha,h))--b;
  while(b<=n && !by_leq(p,n,b,alpha,h))++b;
  return b<=n?b:0;
}

inline double roundoff(int terms,double magnitude) {
  double u=std::numeric_limits<double>::epsilon()/2;
  double nu=(terms+16.)*u;
  return (nu/(1-nu))*std::max(1.,magnitude);
}
inline bool e_pass(double e,int k,int m,double alpha) {
  if(std::isinf(e))return true;
  double lhs=(alpha*e)*k,err=roundoff(3,std::abs(lhs)+m);
  if(lhs-m>err)return true;
  if(lhs-m< -err)return false;
  return cpp_rational(alpha)*cpp_rational(e)*k>=m;
}
inline double e_cap(int m,double alpha) {
  double cap=m/alpha;
  if(cpp_rational(alpha)*cpp_rational(cap)<m)
    cap=std::nextafter(cap,std::numeric_limits<double>::infinity());
  return cap;
}
}
#endif
