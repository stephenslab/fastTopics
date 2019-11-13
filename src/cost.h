#ifndef INCLUDE_COST
#define INCLUDE_COST

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec cost        (const arma::mat& X, const arma::mat& A,
		       const arma::mat& B, double e, bool poisson);
arma::vec cost_sparse (const arma::sp_mat& X, const arma::mat& A,
	   	       const arma::mat& B, double e, bool poisson);

#endif
