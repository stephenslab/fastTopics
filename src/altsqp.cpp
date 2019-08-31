#include "misc.h"
#include "mixem.h"
#include "altsqp.h"

using namespace arma;

// Run one EM update and SQP update on the modified problem, then
// recover the updated solution to the unmodified problem.
void altsqp_update_em_sqp (mat& B, vec& w, const vec& bs, double ws, vec& x,
			   double e, uint numem, uint numsqp,
			   mixsqp_control_params control) {

  // Update the solution to the modified problem.
  uint n = B.n_rows;
  vec  y = get_modified_problem_params(B,w,bs,ws,x);
  vec  ev(n);
  ev.fill(e);
  if (numem > 0)
    mixem(B,w,x,ev,numem);
  if (numsqp > 0)
    mixsqp(B,w,x,ev,numsqp,control,false);

  // Recover the updated solution to the unmodified problem.
  x %= y;
}

// Set up the modified problem solved using either EM or the mix-SQP
// algorithm.
vec get_modified_problem_params (mat& B, vec& w, const vec& bs,
				 double ws, vec& x) {
  uint m = B.n_cols;
  vec y(m);
  y.fill(ws);
  y /= bs;
  w /= ws;
  x /= y;
  scalecols(B,y);
  return y;
}
