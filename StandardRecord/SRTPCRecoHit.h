#ifndef CAFSRTPCRECOHIT_H
#define CAFSRTPCRECOHIT_H

#include <vector>

namespace caf
{
  class SRTPCRecoHit
  {
  public:
    SRTPCRecoHit();

    float x, y, z;
    float sig;
    float rms;
    unsigned int chan;
  };
}

#endif
