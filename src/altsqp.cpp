#include <RcppArmadillo.h>
#include "misc.h"
#include "mixsqp.h"

// This is needed to tell R where to find the additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
double altsqp_main_loop_rcpp (const arma::mat& X, arma::mat& F, arma::mat& Fn,
			      arma::mat& Fy, arma::mat& Fbest, arma::mat& L,
			      arma::mat& Ln, arma::mat& Ly, arma::mat& Lbest,
			      double f, double fbest, const arma::vec& xsrow,
			      const arma::vec& xscol, double best,
			      double betamax, int numiter, int nc,
			      int extrapolate, double betainit,
			      double betaincrease, double betareduce,
			      double betamaxincrease, double e,
			      arma::mat& progress, bool verbose) {

  // Quantities used in the computations below.
  double f0;
  double beta, beta0;
  double d;
  
  if (verbose)
    Rprintf("Running altsqp_main_loop_cpp\n");

  // Iteratively apply the EM and SQP updates
  for (int iter = 0; iter < numiter; iter++) {
    
    // Store the value of the objective at the current iterate.
    f0 = f;

    // When the time is right, initiate the extrapolation scheme.
    if (beta == 0 && iter >= extrapolate) {
      beta  = betainit;
      beta0 = betainit;
    }

    // UPDATE FACTORS
    // --------------
    // Update the factors ("basis vectors").
    // TO DO.

    // Compute the extrapolated update for the factors. Note that
    // when beta = 0, Fy = Fn.
    // TO DO.

    // UPDATE LOADINGS
    // ---------------
    // Update the loadings ("activations").
    // TO DO.

    // Compute the extrapolated update for the loadings. Note that
    // when beta = 0, Ly = Ln.
    // TO DO.

    // Compute the value of the objective (cost) function at the
    // extrapolated solution for the factors (F) and the
    // non-extrapolated solution for the loadings (L).
    // TO DO.

    if (beta == 0) {

      // No extrapolation is used, so use the basic coordinate-wise
      // updates for the factors and loadings.
      F = Fn;
      L = Ln;
    } else {

      // Update the extrapolation parameters following Algorithm 3 of
      // Ang & Gillis (2019).
      if (f > f0) {

        // The solution did not improve, so restart the extrapolation
        // scheme.
        Fy       = F;
        Ly       = L;
        betamax  = beta0;
        beta    *= betareduce;
      } else {
        
        // The solution is improved; retain the basic co-ordinate ascent
        // update as well.
        F        = Fn;
        L        = Ln;
	beta    *= betaincrease;
        beta     = min(betamax,beta);
        beta0    = beta;
	betamax *= betamaxincrease;
        betamax  = min(0.99,betamax);
      }
    }
  
    // If the solution is improved, update the current best solution
    // using the extrapolated estimates of the factors (F) and the
    // non-extrapolated estimates of the loadings (L).
    if (f < fbest) {
      
      // TO DO.
      Fbest = Fy;
      Lbest = Ln;
      fbest = f;
    } else
      d = 0;

    // Record the algorithm's progress.
    // TO DO.
  }
  
  return fbest;
}
