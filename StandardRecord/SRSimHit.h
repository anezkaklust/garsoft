#ifndef CAFSRSIMHIT_H
#define CAFSRSIMHIT_H

#include <cstddef>

namespace caf
{
  class SRSimHit
  {
  public:
    SRSimHit();

    float x, y, z;
    float t;
    float E;
    int trkid;
    size_t cellid;
  };
}

#endif
