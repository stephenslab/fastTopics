#ifndef INCLUDE_POISMIXEM
#define INCLUDE_POISMIXEM

#include "mixsqp.h"

// FUNCTION DECLARATIONS
// ---------------------
arma::vec poismixem (const arma::mat& L, const arma::vec& w, 
		     const arma::vec& x0, uint numiter);

arma::vec poismixsqp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, uint numiter,
		      const mixsqp_control_params& control);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w, 
		arma::vec& x, arma::mat& P, uint numiter);

void poismixsqp (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		 arma::vec& x, arma::mat& Z, arma::mat& H, uint numiter,
		 const mixsqp_control_params& control);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		const arma::uvec& i, arma::vec& x, uint numiter);

void poismixsqp (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		 const arma::uvec& i, arma::vec& x, arma::mat& H, uint numiter,
		 const mixsqp_control_params& control);

arma::vec scd_kl_update (const arma::mat& L, const arma::vec& w,
			 const arma::vec& x0, uint numiter, double e);

void poismix_one_nonzero (const arma::mat& L1, const arma::vec& u, 
			  const arma::vec& w, uint i, arma::vec& x);

#endif
