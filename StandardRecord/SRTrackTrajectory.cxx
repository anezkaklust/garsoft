#include "garsoft/StandardRecord/SRTrackTrajectory.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRTrackTrajectory::SRTrackTrajectory() :
    len(kNaN), chi2(kNaN), avgion(kNaN)
  {
  }
}
