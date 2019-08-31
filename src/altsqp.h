#ifndef INCLUDE_ALTSQP
#define INCLUDE_ALTSQP

#include "mixsqp.h"
#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec get_modified_problem_params (arma::mat& B, arma::vec& w,
				       const arma::vec& bs, double ws,
				       arma::vec& x);
void altsqp_update_em_sqp (arma::mat& B, arma::vec& w, const arma::vec& bs,
			   double ws, arma::vec& x, double e, uint numem,
			   uint numsqp, mixsqp_control_params control);

#endif
