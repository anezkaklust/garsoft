#include "garsoft/StandardRecord/StandardRecord.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  StandardRecord::StandardRecord()
    : event(-1), run(-1), subrun(-1),
      TPC_X(kNaN), TPC_Y(kNaN), TPC_Z(kNaN),
      SimnHits(0),
      SimEnergySum(kNaN),
      SimnHits_MuID(0),
      SimEnergySum_MuID(kNaN),
      DiginHits(0),
      DiginHits_MuID(0),
      ReconHits(0),
      RecoEnergySum(kNaN),
      ReconHits_MuID(0),
      RecoEnergySum_MuID(kNaN),
      nCluster(0)
  {
  }
}
