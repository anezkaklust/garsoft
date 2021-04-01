#ifndef CAFSRVEE_H
#define CAFSRVEE_H

#include "StandardRecord/SRVeeHypothesis.h"

#include <cstddef>

namespace caf
{
  class SRVee
  {
  public:
    SRVee();

    size_t id;
    float  x, y, z;
    size_t t;

    SRVeeHypothesis Kpipi, Lppi, Lpip;
  };
}

#endif
