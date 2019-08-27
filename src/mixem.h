#ifndef INCLUDE_MIXEM
#define INCLUDE_MIXEM

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void mixem (const arma::mat& L, const arma::vec& w, arma::vec& x,
	    const arma::vec& e, uint numiter);

#endif
