//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R.h>

#define TINY_NUM 1e-16
#define MAX_ITER 500
#define NNMF_INNER_MAX_ITER 10
#define TRACE_STEP 10

using namespace arma;

int update (mat& H, const mat& Wt, const mat& A, const umat& mask,
	    uint max_iter);

int update_with_missing (mat& H, const mat& Wt, const mat& A,
			 const umat& mask, uint max_iter);

int scd_kl_update(subview_col<double> Hj, const mat& Wt, const vec& Aj,
		  const vec& sumW, const subview_col<uword> mask,
		  uint max_iter);
