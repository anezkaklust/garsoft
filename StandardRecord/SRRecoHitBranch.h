#ifndef CAFSRRECOHITBRANCH_H
#define CAFSRRECOHITBRANCH_H

#include "StandardRecord/SRTPCRecoHit.h"
#include "StandardRecord/SRRecoHit.h"

#include <vector>

namespace caf
{
  class SRRecoHitBranch
  {
  public:
    SRRecoHitBranch();

    SRTPCRecoHit tpc;
    SRRecoHit ecal;
    SRRecoHit muid;
  };
}

#endif
