#include "nnlm.h"

int scd_kl_update(subview_col<double> Hj, const mat& Wt, const vec& Aj,
		  const vec& sumW, const subview_col<uword> mask,
		  uint max_iter) {
  
	// Problem:  Aj = W * Hj
	// Method: Sequentially minimize KL distance using quadratic approximation
	// Wt = W^T
	// sumW = column sum of W
	// mask: skip updating

	double sumHj = sum(Hj);
	vec Ajt = Wt.t()*Hj;
	vec mu;
	double a; // 2nd-order-derivative
	double b; // 1st-order-derivative
	double tmp;
	bool is_masked = mask.n_elem > 0;

	uint t = 0; 
	for (; t < max_iter; t++)
	{
		for (uint k = 0; k < Wt.n_rows; k++)
		{
			if (is_masked && mask(k) > 0) continue;
			mu = Wt.row(k).t()/(Ajt + TINY_NUM);
			a = dot(Aj, square(mu));
			b = dot(Aj, mu) - sumW(k); // 0.5*ax^2 - bx
			b += a*Hj(k);
			tmp = b/(a+TINY_NUM); 
			if (tmp < 0) tmp = 0;
			if (tmp != Hj(k))
			{
				Ajt += (tmp - Hj(k)) * Wt.row(k).t();
				sumHj += tmp - Hj(k);
				Hj(k) = tmp;
			}
		}
	}
	return int(t);
}
