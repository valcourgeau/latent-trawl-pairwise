#include <cstdio>
#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP modEx2(SEXP ns, SEXP xs) {
    int n = as<int>(ns);
    double x = as<double>(xs);
    for (int i=0; i<n; i++)
        x=1/(1+x);
    return wrap(x);
}


int main(void) {
    printf("Hello, World!\n");
}
