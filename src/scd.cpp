#include <RcppParallel.h>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void scd_update (mat& H, const mat& Wt, const mat& A, uint numiter, double e);

void scd_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		    const vec& sumw, uint numiter, double e);

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// loadings updates in scd_update_factors_parallel_rcpp.
//
// [[Rcpp::depends(RcppParallel)]]
struct scd_factor_updater : public RcppParallel::Worker {
  const mat& A;
  const mat& Wt;
  mat&       H;
  const vec& sumw;
  uint       numiter;
  double     e;

  // This is used to create a scd_factor_updater object.
  scd_factor_updater (const mat& A, const mat& Wt, mat& H, const vec& sumw, 
		      uint numiter, double e) :
    A(A), Wt(Wt), H(H), sumw(sumw), numiter(numiter), e(e) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint j = begin; j < end; j++)
      scd_kl_update(H.col(j),Wt,A.col(j),sumw,numiter,e);
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform a sequential co-ordinate descent (SCD) update for the
// factors matrix (H).
//
// The inputs are: A, the m x n matrix to be decomposed, in which A is
// approximated as crossprod(H,Wt); Wt, the initial estimate of the k x
// n loadings matrix; H, the k x m factors matrix; numiter, the number
// of inner-loop iterations to perform; and e, a non-negative scalar
// included in the computations to prevent NaNs due to division by zero.
//
// This implementation is adapted from the R and C++ code developed by
// Xihui Lin and Paul Boutros, which is available for download at
// https://github.com/linxihui/NNLM.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_rcpp (const arma::mat& A, const arma::mat& Wt,
				   const arma::mat& H, uint numiter,
				   double e) {
  mat Hnew = H;
  scd_update(Hnew,Wt,A,numiter,e);
  return Hnew;
}

// This does the same thing as scd_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_parallel_rcpp (const arma::mat& A, 
					    const arma::mat& Wt,
					    const arma::mat& H, 
					    uint numiter, double e) {
  mat Hnew = H;
  vec sumw = sum(Wt,1);
  scd_factor_updater worker(A,Wt,Hnew,sumw,numiter,e);
  parallelFor(0,A.n_cols,worker);
  return Hnew;
}

// Iterate the SCD updates over all columns of H.
void scd_update (mat& H, const mat& Wt, const mat& A, uint numiter, double e) {
  vec sumw = sum(Wt,1);
  for (uint j = 0; j < A.n_cols; j++) 
    scd_kl_update(H.col(j),Wt,A.col(j),sumw,numiter,e);
}

// Implements the core part of the sequential co-ordinate descent
// (SCD) updates.
void scd_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		    const vec& sumw, uint numiter, double e) {
  uint   k   = Wt.n_rows;
  vec    Ajt = Wt.t()*Hj;
  vec    mu;
  double a;
  double b;
  double x;

  for (uint t = 0; t < numiter; t++)
    for (uint i = 0; i < k; i++) {
      mu = Wt.row(i).t()/(Ajt + e);
      a  = dot(Aj,square(mu));
      b  = dot(Aj,mu) - sumw(i) + a*Hj(i);
      x  = b/(a + e);
      x  = maximum(x,0);
      if (x != Hj(i)) {
	Ajt   += (x - Hj(i)) * Wt.row(i).t();
	Hj(i)  = x;
      }
    }
}
