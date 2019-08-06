#ifndef INCLUDE_COST
#define INCLUDE_COST

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double cost        (const arma::mat& X, const arma::mat& A,
		    const arma::mat& B, double e);
double cost_sparse (const arma::sp_mat& X, const arma::mat& A,
		    const arma::mat& B, double e);

#endif
