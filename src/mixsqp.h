#ifndef INCLUDE_MIXSQP
#define INCLUDE_MIXSQP

#include <RcppArmadillo.h>

// TYPE DEFINITIONS
// ----------------

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
  int    maxiteractiveset;
} mixsqp_control_params;

// FUNCTION DECLARATIONS
// ---------------------
void mixsqp (const arma::mat& L, const arma::vec& w, arma::vec& x,
	     const arma::vec& e, int numiter, mixsqp_control_params control,
	     bool verbose);

#endif
