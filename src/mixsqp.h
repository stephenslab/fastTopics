#ifndef INCLUDE_MIXSQP
#define INCLUDE_MIXSQP

#include <RcppArmadillo.h>

// TYPE DEFINITIONS
// ----------------
// A list of parameters controlling behaviour of the mix-SQP algorithm.
typedef struct {
  double convtolactiveset;
  double zerothresholdsolution;
  double zerothresholdsearchdir;
  double suffdecr;
  double stepsizereduce;
  double minstepsize;
  double identitycontribincrease;
  uint   maxiteractiveset;
  double e;
} mixsqp_control_params;

// FUNCTION DECLARATIONS
// ---------------------
arma::vec mixsqp (const arma::mat& L, const arma::vec& w, const arma::vec& x0,
		  uint numiter, const mixsqp_control_params& control,
		  arma::vec& obj);

void mixsqp (const arma::mat& L, const arma::vec& w, arma::vec& x,
	     arma::mat& Z, arma::mat& H, uint numiter,
	     const mixsqp_control_params& control, arma::vec& obj);

#endif
