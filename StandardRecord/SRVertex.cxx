#include "garsoft/StandardRecord/SRVertex.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRVertex::SRVertex()
    : id(-1), x(kNaN), y(kNaN), z(kNaN), t(-1), ntrks(0), Q(0)
  {
  }
}
