#include "StandardRecord/SRTrack.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRTrack::SRTrack() :
    id(-1), startq(0), endq(0), NTPCClustersOnTrack(0), mcidx(-1), mcfrac(kNaN)
  {
  }
}
