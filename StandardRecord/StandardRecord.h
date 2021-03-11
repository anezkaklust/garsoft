#ifndef CAFSTANDARDRECORD_H
#define CAFSTANDARDRECORD_H

#include "garsoft/StandardRecord/SRHeader.h"
#include "garsoft/StandardRecord/SRTruth.h"
#include "garsoft/StandardRecord/SRTrack.h"
#include "garsoft/StandardRecord/SRVertex.h"
#include "garsoft/StandardRecord/SRVee.h"
#include "garsoft/StandardRecord/SRDigitBranch.h"
#include "garsoft/StandardRecord/SRRecoHitBranch.h"
#include "garsoft/StandardRecord/SRClusterBranch.h"

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

    SRTrack trk;

    SRVertex vtx;

    SRVee vee;

    SRDigitBranch dig;

    SRRecoHitBranch hit;

    SRClusterBranch clust;
  };
}

#endif
