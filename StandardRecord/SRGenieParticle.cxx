#include "garsoft/StandardRecord/SRGenieParticle.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRGenieParticle::SRGenieParticle() :
    intidx(-1), idx(-1), pdg(0), status(-1),
    firstmom(-1), lastmom(-1), firstdaugh(-1), lastdaugh(-1),
    E(kNaN), mass(kNaN)
  {
  }
}
