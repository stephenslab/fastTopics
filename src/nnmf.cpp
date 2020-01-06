#include "nnlm.h"

// Using version 0.4-3 of the NNLM package.

//[[Rcpp::export]]
Rcpp::List c_nnmf (const arma::mat& A, const unsigned int k, arma::mat W,
		   arma::mat H, arma::umat Wm, arma::umat Hm,
		   uint max_iter, uint inner_max_iter) {
  
	/******************************************************************************************************
	 *              Non-negative Matrix Factorization(NNMF) using alternating scheme
	 *              ----------------------------------------------------------------
	 * Description:
	 * 	Decompose matrix A such that
	 * 		A = W H
	 * Arguments:
	 * 	A              : Matrix to be decomposed
	 * 	W, H           : Initial matrices of W and H, where ncol(W) = nrow(H) = k. # of rows/columns of W/H could be 0
	 * 	Wm, Hm         : Masks of W and H, s.t. masked entries are no-updated and fixed to initial values
	 * 	max_iter       : Maximum number of iteration
	 * 	inner_max_iter : Maximum number of iterations passed to each inner W or H matrix updating loop

	 * Return:
	 * 	A list (Rcpp::List) of
	 * 		W, H          : resulting W and H matrices
	 ******************************************************************************************************/

	unsigned int n = A.n_rows;
	unsigned int m = A.n_cols;
	unsigned int N_non_missing = n*m;
	uvec non_missing;
	bool any_missing = !A.is_finite();
	if (any_missing)
	{
		non_missing = find_finite(A);
		N_non_missing = non_missing.n_elem;
	}

	if (Wm.empty())
		Wm.resize(0, n);
	else
		inplace_trans(Wm);
	if (Hm.empty())
		Hm.resize(0, m);

	if (W.empty())
	{
		W.randu(k, n);
		W *= 0.01;
		if (!Wm.empty())
			W.elem(find(Wm > 0)).fill(0.0);
	}
	else
		inplace_trans(W);

	if (H.empty())
	{
		H.randu(k, m);
		H *= 0.01;
		if (!Hm.empty())
			H.elem(find(Hm > 0)).fill(0.0);
	}

	int total_raw_iter = 0;
	unsigned int i = 0;
	for(; i < max_iter; i++) {
		if (any_missing)
		{
			// update W
			total_raw_iter += update_with_missing(W, H, A.t(), Wm, inner_max_iter);
			// update H
			total_raw_iter += update_with_missing(H, W, A, Hm, inner_max_iter);
		}
		else
		{
			// update W
			total_raw_iter += update(W, H, A.t(), Wm, inner_max_iter);
			// update H
			total_raw_iter += update(H, W, A, Hm, inner_max_iter);

		}

	}

	return Rcpp::List::create(
		Rcpp::Named("W") = W.t(),
		Rcpp::Named("H") = H);
}
