#ifndef INCLUDE_ALTSQP
#define INCLUDE_ALTSQP

#include "mixem.h"
#include "mixsqp.h"
#include <RcppArmadillo.h>

// FUNCTION DEFINITIONS
// --------------------
// Set up the modified problem solved using either EM or the mix-SQP
// algorithm.
inline arma::vec get_modified_problem_params (arma::mat& B, arma::vec& w,
					      const arma::vec& bs, double ws,
					      arma::vec& x) {
  uint m = B.n_cols;
  arma::vec y(m);
  y.fill(ws);
  y /= bs;
  w /= ws;
  x /= y;
  scalecols(B,y);
  return y;
}

// Run one EM update and SQP update on the modified problem, then
// recover the updated solution to the unmodified problem.
inline void altsqp_update_em_sqp (arma::mat& B, arma::vec& w,
				  const arma::vec& bs, double ws,
				  arma::vec& x, double e, uint numem,
				  uint numsqp,
				  const mixsqp_control_params& control) {

  // Update the solution to the modified problem.
  uint n = B.n_rows;
  arma::vec y = get_modified_problem_params(B,w,bs,ws,x);
  arma::vec ev(n);
  ev.fill(e);
  // if (numem > 0)
  //   mixem(B,w,x,ev,numem);
  if (numsqp > 0)
    mixsqp(B,w,x,ev,numsqp,control,false);

  // Recover the updated solution to the unmodified problem.
  x %= y;
}

#endif
