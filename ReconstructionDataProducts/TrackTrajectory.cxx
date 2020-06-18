//
//  TrackTrajectory.cxx
//  garsoft
//
//  Created by Tom Junk, June 17, 2020

#include "ReconstructionDataProducts/TrackTrajectory.h"

namespace gar {
  namespace rec {

    //--------------------------------------------------------------------------
    TrackTrajectory::TrackTrajectory() {}

    //--------------------------------------------------------------------------
    void TrackTrajectory::setData(std::vector<TVector3> ForwardTrajectory,
                 std::vector<TVector3> BackwardTrajectory) {
        fFWDTrajectory = ForwardTrajectory;
        fBAKTrajectory = BackwardTrajectory;
    }



  } // rec
} // gar
