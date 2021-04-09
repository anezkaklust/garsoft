#ifndef CAFSRDIGIT_H
#define CAFSRDIGIT_H

#include <cstddef>

namespace caf
{
  class SRDigit
  {
  public:
    SRDigit();

    float x, y, z;
    float t;
    unsigned int adc;
    size_t cellid;
  };
}

#endif
