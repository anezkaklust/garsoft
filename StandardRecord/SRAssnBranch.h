#ifndef CAFSRASSNS_H
#define CAFSRASSNS_H

#include "garsoft/StandardRecord/SRVertexAssn.h"
#include "garsoft/StandardRecord/SRVeeAssn.h"

#include <vector>

namespace caf
{
  class SRAssnBranch
  {
  public:
    SRAssnBranch();

    unsigned int nvtx;
    std::vector<SRVertexAssn> vtx;

    unsigned int nvee;
    std::vector<SRVeeAssn> vee;
  };
}

#endif
