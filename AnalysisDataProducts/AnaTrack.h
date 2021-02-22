#ifndef GAR_ANATRACK_H
#define GAR_ANATRACK_H

#include "ReconstructionDataProducts/Track.h"
#include <vector>

using std::pair;
using std::vector;
using gar::rec::Track;

namespace gar {
  namespace adp {
    class AnaTrack {

      private:
        Track           fTrk;  // copy all reco::Track info (for now...)
        vector<pair<int,float>> fPidF; ///< PID: (PDG code, probability) front fit
        vector<pair<int,float>> fPidB; ///<  " backfit
        float           fIonF; ///< average ionization along foward fit
        float           fIonB; ///< average ionization along backward fit

      public:

        // constructors
        AnaTrack() {}
        AnaTrack(const Track& trk, const vector<pair<int,float>>& pidf, 
                 const vector<pair<int,float>>& pidb, float ionf, float ionb) :
          fTrk(trk),
          fPidF(pidf),
          fPidB(pidb),
          fIonF(ionf),
          fIonB(ionb)
          {}

        // getters
        const Track* GetTrack() const { return &fTrk;      }
        const gar::rec::IDNumber TrackID() const { return fTrk.getIDNumber(); }
        const vector<pair<int,float>>*  PidF() const { return &fPidF;  }
        const vector<pair<int,float>>*  PidB() const { return &fPidB;  }
        //const float  PidProbF() const { return fPidF.second; }
        //const float  PidProbB() const { return fPidB.second; }
        const float  IonizF()   const { return fIonF;        }
        const float  IonizB()   const { return fIonB;        }
    };//class
  }//adp
}//gar

#endif
