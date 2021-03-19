#include "StandardRecord/SRParticle.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRParticle::SRParticle() :
    trkid(-1), pdg(0), momidx(-1), momtrkid(-1), mompdg(-1),
    time(kNaN), vtxidx(-1), ntraj(0)
  {
  }
}
