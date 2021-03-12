#include "garsoft/StandardRecord/SRTrackPID.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRTrackPID::SRTrackPID() :
    id(-1), prob(kNaN)
  {
  }
}
