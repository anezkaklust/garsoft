#include "StandardRecord/SRHeader.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRHeader::SRHeader()
    : event(-1), run(-1), subrun(-1),
      TPC_X(kNaN), TPC_Y(kNaN), TPC_Z(kNaN)
  {
  }
}
