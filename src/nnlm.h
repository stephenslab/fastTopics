//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R.h>

#define TINY_NUM 1e-16

using namespace arma;

void update (mat& H, const mat& Wt, const mat& A, uint max_iter);

void scd_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		    const vec& sumW, uint max_iter);
