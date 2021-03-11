#ifndef CAFSRCLUSTERBRANCH_H
#define CAFSRCLUSTERBRANCH_H

#include "garsoft/StandardRecord/SRTPCCluster.h"
#include "garsoft/StandardRecord/SRECalCluster.h"

#include <vector>

namespace caf
{
  class SRClusterBranch
  {
  public:
    SRClusterBranch();

    SRTPCCluster tpc;
    SRECalCluster ecal;
  };
}

#endif
