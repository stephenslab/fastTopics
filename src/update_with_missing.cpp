#include "nnlm.h"

int update (mat& H, const mat & Wt, const mat& A, const umat& mask,
	    uint max_iter) {
	// A = W H, solve H
	// No missing in A, Wt = W^T

	unsigned int m = A.n_cols;
	int total_raw_iter = 0;

	bool is_masked = !mask.empty();
	mat WtW;
	vec mu, sumW;
	sumW = sum(Wt, 1);

	for (uint j = 0; j < m; j++) // by columns of H
	{
		// break if all entries of col_j are masked
		if (is_masked && arma::all(mask.col(j)))
			continue;

		int iter = 0;
		iter = scd_kl_update(H.col(j), Wt, A.col(j), sumW, mask.col(j), max_iter);
		total_raw_iter += iter;
	}
	return total_raw_iter;
}

int update_with_missing (mat& H, const mat& Wt, const mat& A, const umat& mask,
			 uint max_iter) {
	// A = W H, solve H
	// With missings in A, Wt = W^T

	uint m = A.n_cols;
	uint total_raw_iter = 0;

	bool is_masked = mask.n_elem > 0;
	mat WtW;
	vec mu;

	for (uint j = 0; j < m; j++) // by columns of H
	{
		// break if all entries of col_j are masked
		if (is_masked && arma::all(mask.col(j)))
			continue;

		bool any_missing = !is_finite(A.col(j));
		uvec non_missing;
		if (any_missing)
			non_missing = find_finite(A.col(j));
	
		int iter = 0;

		total_raw_iter += iter;
	}
	return total_raw_iter;
}
