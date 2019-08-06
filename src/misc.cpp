#include "misc.h"

// FUNCTION DEFINITIONS
// --------------------
// Return a or b, which ever is smaller.
double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}
