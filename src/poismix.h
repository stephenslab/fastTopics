#ifndef INCLUDE_POISMIXEM
#define INCLUDE_POISMIXEM

#include "mixsqp.h"

// FUNCTION DECLARATIONS
// ---------------------
arma::vec poismixem (const arma::mat& L, const arma::vec& w, 
		     const arma::vec& x0, unsigned int numiter);

arma::vec poismixsqp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, unsigned int numiter,
		      const mixsqp_control_params& control);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w, 
		arma::vec& x, arma::mat& P, unsigned int numiter);

void poismixsqp (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		 arma::vec& x, arma::mat& Z, arma::mat& H, 
		 unsigned int numiter, const mixsqp_control_params& control);

void poismixem (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		const arma::uvec& i, arma::vec& x, unsigned int numiter);

void poismixsqp (const arma::mat& L1, const arma::vec& u, const arma::vec& w,
		 const arma::uvec& i, arma::vec& x, arma::mat& H, 
		 unsigned int numiter, const mixsqp_control_params& control);

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

void poismix_one_nonzero (const arma::mat& L1, const arma::vec& u, 
			  const arma::vec& w, unsigned int i, 
			  arma::vec& x);

#endif
