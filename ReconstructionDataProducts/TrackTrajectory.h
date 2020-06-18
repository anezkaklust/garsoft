//
//  TrackIoniz.h
//
//  Created by Tom Junk June 17, 2020.
//  Track trajectory points.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_TrackTrajectory_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_TrackTrajectory_h

#include <vector>
#include <TVector3.h>

namespace gar {
  namespace rec {

    class TrackTrajectory {

    public:
      TrackTrajectory();
      inline const std::vector<TVector3> getFWDTrajectory() const {return fFWDTrajectory;};
      inline const std::vector<TVector3> getBAKTrajectory() const {return fBAKTrajectory;};

      // let the compiler provide the dtor

    private:
      std::vector<TVector3> fFWDTrajectory;  ///< Segment values ordered by forward fit
      std::vector<TVector3> fBAKTrajectory;  ///< Segment values ordered by backward fit

#ifndef __GCCXML__

    public:

      void setData(std::vector<TVector3> ForwardTrajectory,std::vector<TVector3> BackwardTrajectory);

#endif

    };

  } // rec
} // gar


#endif
