#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double min     (double a, double b);
void   copycol (const arma::mat& X, int i, arma::vec& y);

#endif
