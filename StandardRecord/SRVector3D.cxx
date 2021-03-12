#include "garsoft/StandardRecord/SRVector3D.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRVector3D::SRVector3D() :
    x(kNaN), y(kNaN), z(kNaN)
  {
  }

  SRVector3D::SRVector3D(float _x, float _y, float _z) :
    x(_x), y(_y), z(_z)
  {
  }
}
