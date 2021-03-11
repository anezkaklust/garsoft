#include "garsoft/StandardRecord/SRTruth.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRTruth::SRTruth() :
    SimnHits(0),
    SimEnergySum(kNaN),
    SimnHits_MuID(0),
    SimEnergySum_MuID(kNaN)
  {
  }
}
