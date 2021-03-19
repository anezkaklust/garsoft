#include "StandardRecord/SRTPCRecoHit.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRTPCRecoHit::SRTPCRecoHit()
    : x(kNaN), y(kNaN), z(kNaN), sig(kNaN), rms(kNaN), chan(-1)
  {
  }
}
