#ifndef CAFSRRECOHIT_H
#define CAFSRRECOHIT_H

#include <cstddef>

namespace caf
{
  class SRRecoHit
  {
  public:
    SRRecoHit();

    size_t id;
    float  x, y, z;
    float  t;
    float  E;
    size_t cellid;
  };
}

#endif
