#ifndef INCLUDE_MIXEM
#define INCLUDE_MIXEM

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec mixem (const arma::mat& L, const arma::vec& w, const arma::vec& x0,
		 uint numiter);

void mixem (const arma::mat& L, const arma::vec& w, arma::vec& x, arma::mat& P,
	    uint numiter);

void mixem_update (const arma::mat& L, const arma::vec& w, arma::vec& x,
		   arma::mat& P);
  
#endif
