#include <cstdio>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

const double EPSILON = std::numeric_limits<double>::epsilon();

// [[Rcpp::export]]
double CppCaseZeroZero(double alpha, double beta, double kappa, double B1, double B2, double B3) {
  double A(B1+B3);
  double b1(-alpha*B1/A), b2(-alpha*B2/A), b3(-alpha*B3/A);
  double tmp(1.0 - 2.0*pow(1+kappa/beta, -alpha) + pow(1+kappa/beta, b1+b3)*pow(1+2*kappa/beta, b2));

  return tmp;
}

// [[Rcpp::export]]
double CppCaseOneZero(NumericVector xs, double alpha, double beta, double kappa, double B1, double B2, double B3) {
  double A(B1+B3);
  double b1(-alpha*B1/A), b2(-alpha*B2/A);
  double tmp(alpha/beta * pow(1.0+kappa/beta, -alpha-1.0));

  if (std::abs(xs[0]) < EPSILON){
    xs[0] = xs[1];
    xs[1] = 0.0;
  }
  double x(xs[0]);
  
  tmp += 1.0/beta*pow(1.0+(kappa+x)/beta, b1-1.0)*pow(1.0+(2.0*kappa+x)/beta, b2-1.0)*pow(1.0+kappa/beta, b1)*(-alpha*(1.0+(kappa+x)/beta)+b1*kappa/beta);

  return tmp;
}

// [[Rcpp::export]]
double CppCaseOneOne(NumericVector xs, double alpha, double beta, double kappa, double B1, double B2, double B3) {
  double A(B1+B3);
  double b1(-alpha*B1/A), b2(-alpha*B2/A), b3(-alpha*B3/A);
  double x1(xs[1]), x2(xs[2]);
  double tmp_1(1.0/pow(beta, 2.0)*pow(1.0+(2.0*kappa+x1+x2), b2-2.0)*pow(1.0+(kappa+x1), b1-1.0)*pow(1+(kappa+x2), b3-1.0));
  double tmp_2(b2*(b2-1.0)*(1.0+(kappa+x1)/beta)*(1.0+(kappa+x2)/beta)+pow(b1, 2.0)*pow(1.0+(2.0*kappa+x1+x2)/beta, 2.0));
  double tmp_3(b2*b1*(1.0+(2*kappa+x1+x2))*((1.0+(kappa+x1)/beta)+(1.0+(kappa+x2)/beta)));
  
  return tmp_1 * (tmp_2 + tmp_3);
}
