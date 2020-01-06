#include "nnlm.h"

void scd_kl_update(subview_col<double> Hj, const mat& Wt, const vec& Aj,
		   const vec& sumW, uint max_iter) {
  
  // Problem:  Aj = W * Hj
  // Method: Sequentially minimize KL distance using quadratic approximation
  // Wt = W^T
  // sumW = column sum of W

  double sumHj = sum(Hj);
  vec Ajt = Wt.t()*Hj;
  vec mu;
  double a; // 2nd-order-derivative
  double b; // 1st-order-derivative
  double tmp;
  
  for (uint t = 0; t < max_iter; t++) {
    for (uint k = 0; k < Wt.n_rows; k++) {
      mu = Wt.row(k).t()/(Ajt + TINY_NUM);
      a  = dot(Aj, square(mu));
      b  = dot(Aj, mu) - sumW(k);
      b  += a*Hj(k);
      tmp = b/(a+TINY_NUM); 
      if (tmp < 0)
	tmp = 0;
      
      if (tmp != Hj(k)) {
	  Ajt += (tmp - Hj(k)) * Wt.row(k).t();
	  sumHj += tmp - Hj(k);
	  Hj(k) = tmp;
	}
    }
  }
}
