#include "StandardRecord/SRCovMx.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRCovMx::SRCovMx() :
    xx(kNaN), xy(kNaN), xz(kNaN), yy(kNaN), yz(kNaN), zz(kNaN)
  {
  }
}
