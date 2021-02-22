#ifndef GAR_MCDISPLAYTRACK_H
#define GAR_MCDISPLAYTRACK_H


#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include <vector>
#include <TLorentzVector.h>

namespace gar{
  namespace adp{

    class MCDisplayTrack {

      private:
        std::vector<TLorentzVector> fPositions;
        void Sparsify(simb::MCTrajectory traj);

      public:
        MCDisplayTrack() {}
        MCDisplayTrack(simb::MCTrajectory traj){
            Sparsify(traj);
        }
        MCDisplayTrack(const simb::MCParticle& m){
            simb::MCTrajectory traj = m.Trajectory();
            Sparsify(traj);
        }


        size_t NPoints() const;
        float X(const size_t i) const;
        float Y(const size_t i) const;
        float Z(const size_t i) const;
        float T(const size_t i) const;

    };

  }//ns adp
}//ns gar

#endif
