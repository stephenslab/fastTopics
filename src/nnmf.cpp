#include <RcppArmadillo.h>

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void nnmf_update (mat& H, const mat& Wt, const mat& A, uint max_iter, 
		  double e);

void nnmf_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		     const vec& sumW, uint max_iter, double e);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List nnmf (const arma::mat& A, const uint k, arma::mat W, arma::mat H,
		 uint max_iter, uint inner_max_iter, double e) {
  for (uint i = 0; i < max_iter; i++) {
    
    // Update W.
    nnmf_update(W,H,A.t(),inner_max_iter,e);
    
    // Update H.
    nnmf_update(H,W,A,inner_max_iter,e);
  }
  return Rcpp::List::create(Rcpp::Named("W") = W,Rcpp::Named("H") = H);
}

void nnmf_update (mat& H, const mat & Wt, const mat& A, uint max_iter, 
		  double e) {
  uint m = A.n_cols;
  mat WtW;
  vec mu;
  const vec sumW = sum(Wt,1);
  for (uint j = 0; j < m; j++)
    nnmf_kl_update(H.col(j),Wt,A.col(j),sumW,max_iter,e);
}

void nnmf_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		     const vec& sumW, uint max_iter, double e) {
  
  double sumHj = sum(Hj);
  vec Ajt = Wt.t()*Hj;
  vec mu;
  double a;
  double b;
  double x;
  
  for (uint t = 0; t < max_iter; t++) {
    for (uint k = 0; k < Wt.n_rows; k++) {
      mu = Wt.row(k).t()/(Ajt + e);
      a  = dot(Aj, square(mu));
      b  = dot(Aj, mu) - sumW(k);
      b  += a*Hj(k);
      x  = b/(a + e); 
      if (x < 0)
	x = 0;
      
      if (x != Hj(k)) {
	  Ajt += (x - Hj(k)) * Wt.row(k).t();
	  sumHj += x - Hj(k);
	  Hj(k) = x;
	}
    }
  }
}
