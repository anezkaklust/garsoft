#ifndef CAFSRTPCCLUSTER_H
#define CAFSRTPCCLUSTER_H

#include "StandardRecord/SRCovMx.h"

#include <cstddef>

namespace caf
{
  class SRTPCCluster
  {
  public:
    SRTPCCluster();

    float   x, y, z;
    float   sig;
    float   rms;
    size_t  trkid;
    SRCovMx cov;
  };
}

#endif
