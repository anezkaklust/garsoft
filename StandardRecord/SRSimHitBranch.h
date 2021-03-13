#ifndef CAFSRSIMHITBRANCH_H
#define CAFSRSIMHITBRANCH_H

#include "garsoft/StandardRecord/SRSimHit.h"

#include <vector>

namespace caf
{
  class SRSimHitBranch
  {
  public:
    SRSimHitBranch();

    unsigned int necal;
    std::vector<SRSimHit> ecal;
    float ecaltotE;

    unsigned int nmuid;
    std::vector<SRSimHit> muid;
    float muidtotE;
  };
}

#endif
