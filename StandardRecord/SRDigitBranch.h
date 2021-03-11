#ifndef CAFSRDIGITBRANCH_H
#define CAFSRDIGITBRANCH_H

#include "garsoft/StandardRecord/SRDigit.h"

namespace caf
{
  class SRDigitBranch
  {
  public:
    SRDigitBranch();

    SRDigit ecal;
    SRDigit muid;
  };
}

#endif
