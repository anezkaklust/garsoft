#include "StandardRecord/SRDigit.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRDigit::SRDigit()
    : x(kNaN), y(kNaN), z(kNaN), t(kNaN), adc(-1), cellid(-1)
  {
  }
}
