#ifndef INCLUDE_POISMIXEM
#define INCLUDE_POISMIXEM

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec poismixem (const arma::mat& L, const arma::vec& w, 
		     const arma::vec& x0, unsigned int numiter);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w, 
		arma::vec& x, arma::mat& P, unsigned int numiter);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		const arma::uvec& i, arma::vec& x, unsigned int numiter);

arma::vec scd_kl_update (const arma::mat& L, const arma::vec& w,
			 const arma::vec& x0, unsigned int numiter, 
			 double e);

arma::vec scd_kl_update (const arma::mat& L, const arma::vec& u,
			 const arma::vec& w, const arma::vec& x0,
			 unsigned int numiter, double e);

arma::vec ccd_kl_update (const arma::mat& L, const arma::vec& w,
			 const arma::vec& x0, double e);

arma::vec ccd_kl_update (const arma::mat& L, const arma::vec& u,
			 const arma::vec& w, const arma::vec& x0,
			 double e);

#endif
