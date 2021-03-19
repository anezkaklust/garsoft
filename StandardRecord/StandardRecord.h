#ifndef CAFSTANDARDRECORD_H
#define CAFSTANDARDRECORD_H

#include "StandardRecord/SRHeader.h"
#include "StandardRecord/SRTruth.h"
#include "StandardRecord/SRTrack.h"
#include "StandardRecord/SRVertex.h"
#include "StandardRecord/SRVee.h"
#include "StandardRecord/SRDigitBranch.h"
#include "StandardRecord/SRRecoHitBranch.h"
#include "StandardRecord/SRClusterBranch.h"
#include "StandardRecord/SRAssnBranch.h"

#include <string>
#include <vector>

namespace caf
{
  class StandardRecord
  {
  public:
    StandardRecord();

    SRHeader hdr; ///< global event info

    SRTruth mc;

    unsigned int ntrk;
    std::vector<SRTrack> trk;

    unsigned int nvtx;
    std::vector<SRVertex> vtx;

    unsigned int nvee;
    std::vector<SRVee> vee;

    SRDigitBranch dig;

    SRRecoHitBranch hit;

    SRClusterBranch clust;

    SRAssnBranch assn;
  };
}

#endif
