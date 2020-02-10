#ifndef INCLUDE_MIXEM
#define INCLUDE_MIXEM

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec mixem (const arma::mat& L, const arma::vec& w, const arma::vec& x0,
		 unsigned int numiter);

void mixem (const arma::mat& L1, const arma::vec& w, arma::vec& x,
	    arma::mat& P, unsigned int numiter);

void mixture_one_nonzero (const arma::mat& L1, unsigned int i, arma::vec& x);

#endif
