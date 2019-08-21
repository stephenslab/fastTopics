#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform a single EM update.
void mixem_update (const mat& L, const vec& w, vec& x, const vec& e) {
  uint n = L.n_rows;
  uint m = L.n_cols;
  mat  P(n,m);
  P = L;
  
  // Compute the n x m matrix of posterior mixture assignment
  // probabilities (L is an n x m matrix). This is the "E step".
  scalecols(P,x);
  addtorows(P,e);
  normalizerows(P);
  
  // Update the mixture weights. This is the "M step".
  x = trans(P) * w;
}
