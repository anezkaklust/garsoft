#include "garsoft/StandardRecord/SRParticleTrajectoryPoint.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRParticleTrajectoryPoint::SRParticleTrajectoryPoint() :
    x(kNaN), y(kNaN), z(kNaN), t(kNaN), E(kNaN), px(kNaN), py(kNaN), pz(kNaN)
  {
  }
}
