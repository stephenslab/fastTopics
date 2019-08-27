#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform several EM updates.
void mixem (const mat& L, const vec& w, vec& x, const vec& e, uint numiter) {
  uint n = L.n_rows;
  uint m = L.n_cols;
  mat  P(n,m);

  // Iterate the E and M steps.
  for (uint iter = 0; iter < numiter; iter++) {
  
    // Compute the n x m matrix of posterior mixture assignment
    // probabilities (L is an n x m matrix). This is the "E step".
    P = L;
    scalecols(P,x);
    addtorows(P,e);
    normalizerows(P);
  
    // Update the mixture weights. This is the "M step".
    x = trans(P) * w;
  }
}
