#ifndef CAFSRDIGIT_H
#define CAFSRDIGIT_H

#include <vector>

namespace caf
{
  class SRDigit
  {
  public:
    SRDigit();

    unsigned int                  nhits;
    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<float>            t;
    std::vector<unsigned int>     adc;
    std::vector<size_t>           cellid;
  };
}

#endif
