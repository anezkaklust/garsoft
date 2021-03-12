#ifndef CAFSRDIGITBRANCH_H
#define CAFSRDIGITBRANCH_H

#include "garsoft/StandardRecord/SRDigit.h"

#include <vector>

namespace caf
{
  class SRDigitBranch
  {
  public:
    SRDigitBranch();

    unsigned int necal;
    std::vector<SRDigit> ecal;

    unsigned int nmuid;
    std::vector<SRDigit> muid;
  };
}

#endif
