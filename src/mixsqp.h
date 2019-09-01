#ifndef INCLUDE_MIXSQP
#define INCLUDE_MIXSQP

#include <RcppArmadillo.h>

// TYPE DEFINITIONS
// ----------------
// A list of parameters controlling behaviour of the SQP optimization
// algorithm.
typedef struct {
  double activesetconvtol;
  double zerothreshold;
  double zerosearchdir;
  double suffdecr;
  double stepsizereduce;
  double minstepsize;
  uint   maxiteractiveset;
} mixsqp_control_params;

// FUNCTION DECLARATIONS
// ---------------------
mixsqp_control_params get_mixsqp_control_params (Rcpp::List control);

void mixsqp (const arma::mat& L, const arma::vec& w, arma::vec& x,
	     const arma::vec& e, uint numiter,
	     const mixsqp_control_params& control, bool verbose);

// FUNCTION DEFINITIONS
// --------------------
// Get the parameters controlling the behaviour of the SQP
// optimization algorithm from the named elements in a list.
inline mixsqp_control_params get_mixsqp_control_params (Rcpp::List control) {
  mixsqp_control_params x;
  x.activesetconvtol = control["activesetconvtol"];
  x.zerothreshold    = control["zerothreshold"];
  x.zerosearchdir    = control["zerosearchdir"];
  x.suffdecr         = control["suffdecr"];
  x.stepsizereduce   = control["stepsizereduce"];
  x.minstepsize      = control["minstepsize"];
  x.maxiteractiveset = control["maxiteractiveset"];
  return x;
}

#endif
