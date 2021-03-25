#ifndef CAFSRVEEHYPOTHESIS_H
#define CAFSRVEEHYPOTHESIS_H

#include "StandardRecord/SRVector3D.h"

namespace caf
{
  class SRVeeHypothesis
  {
  public:
    SRVeeHypothesis();

    SRVector3D p;
    float E;
    float m;
  };
}

#endif
