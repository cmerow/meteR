#include <Rcpp.h>
// #include <Rmath.h> 
// #include <stdio.h>
// #include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gsne_la(int G0, int S0, int N0, double LA1, double BETA) {
  int n0=N0, s0=S0, g0=G0, i, j;
  double la1=LA1, beta=BETA;
  double num1=0, num2=0, denom=0;
  for ( i=1; i<=s0; i++ ) {
    for ( j=1; j<=n0; j++ ) {
      num1 = num1 + 1/((double)j) * exp(-i*(la1+beta*j));
      denom = denom + 1/((double)i*(double)j) * exp(-i*(la1+beta*j));
    }
    num2 = num2 + exp(-i*(la1+beta))/(1-exp(-beta*i));
  }
  NumericVector out(2);
  out[0] = num1/denom-(double)s0/(double)g0;
  out[1] = num2/denom-(double)n0/(double)g0;
  return out;
}

// [[Rcpp::export]]
double gsne_Z(double LA1, double LA2, double LA3, int S0, int N0) {
  int n0=N0, s0=S0, i, j;
  double la1=LA1, la2=LA2, la3=LA3;
  double beta=la2+la3;
  double Z=0;
  
  for(i=1; i<=s0; i++) {
    for(j=1; j<=n0; j++) {
      Z = Z + 1/((double)i*(double)j) * exp(-i*(la1+beta*j));
    } 
  }
  
  Z = Z/la3;
  return Z;
}
