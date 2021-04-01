#ifndef CAFSRCLUSTERBRANCH_H
#define CAFSRCLUSTERBRANCH_H

#include "StandardRecord/SRTPCCluster.h"
#include "StandardRecord/SRECalCluster.h"

#include <vector>

namespace caf
{
  class SRClusterBranch
  {
  public:
    SRClusterBranch();

    unsigned int ntpc;
    std::vector<SRTPCCluster> tpc;

    unsigned int necal;
    std::vector<SRECalCluster> ecal;
  };
}

#endif
