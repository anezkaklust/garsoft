#ifndef CAFSRTPCCLUSTER_H
#define CAFSRTPCCLUSTER_H

#include "garsoft/StandardRecord/SRCovMx.h"

#include <vector>

namespace caf
{
  class SRTPCCluster
  {
  public:
    SRTPCCluster();

    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<float>            sig;
    std::vector<float>            rms;
    std::vector<size_t>           trkid;
    std::vector<SRCovMx>          cov;
  };
}

#endif
