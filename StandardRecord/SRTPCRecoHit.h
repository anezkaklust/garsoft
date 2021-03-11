#ifndef CAFSRTPCRECOHIT_H
#define CAFSRTPCRECOHIT_H

#include <vector>

namespace caf
{
  class SRTPCRecoHit
  {
  public:
    SRTPCRecoHit();

    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<float>            sig;
    std::vector<float>            rms;
    std::vector<unsigned int>     chan;
  };
}

#endif
