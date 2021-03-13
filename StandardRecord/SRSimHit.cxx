#include "garsoft/StandardRecord/SRSimHit.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRSimHit::SRSimHit() :
    x(kNaN), y(kNaN), z(kNaN), t(kNaN), E(kNaN), trkid(-1), cellid(-1)
  {
  }
}
