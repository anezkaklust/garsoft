#ifndef CAFSRRECOHIT_H
#define CAFSRRECOHIT_H

#include <vector>

namespace caf
{
  class SRRecoHit
  {
  public:
    SRRecoHit();

    unsigned int                  nhits;
    std::vector<size_t>           id;
    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<float>            t;
    std::vector<float>            E;
    std::vector<size_t>           cellid;
    float                         totE;
  };
}

#endif
