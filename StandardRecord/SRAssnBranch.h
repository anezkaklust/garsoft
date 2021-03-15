#ifndef CAFSRASSNS_H
#define CAFSRASSNS_H

#include "garsoft/StandardRecord/SRVertexAssn.h"
#include "garsoft/StandardRecord/SRVeeAssn.h"
#include "garsoft/StandardRecord/SRClusterAssn.h"

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

    unsigned int necalclust;
    std::vector<SRClusterAssn> ecalclust;
  };
}

#endif
