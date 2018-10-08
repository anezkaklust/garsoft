//
//  CaloDeposit.h
//
//  Created by Eldwan Brianne on 8/23/17.
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h
#define GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h

#include <list>
#include <vector>

namespace gar {
  namespace sdp {


    class CaloDeposit{
    public:

      CaloDeposit();

#ifndef __GCCXML__
      CaloDeposit(int trackID,
                    double t,
                    float  e,
                    float  xpos,
                    float  ypos,
                    float  zpos,
                    float  CaloID,
                    unsigned long long int CellID,
                    unsigned int Layer)
      : fTrackID  (trackID)
      , fTime     (t)
      , fEnergy   (e)
      , fX        (xpos)
      , fY        (ypos)
      , fZ        (zpos)
      , fCaloID   (CaloID)
      , fCellID   (CellID)
      , fLayer    (Layer)
      {}

      int    const& TrackID()   const { return fTrackID;   }
      double const& Time()      const { return fTime;      }
      float  const& Energy()    const { return fEnergy;    }
      float  const& X()         const { return fX;         }
      float  const& Y()         const { return fY;         }
      float  const& Z()         const { return fZ;         }
      int  const& CaloID()      const { return fCaloID;    }
      unsigned long long int  const& CellID()      const { return fCellID;    }
      unsigned int  const& Layer()      const { return fLayer;    }

      bool operator  <(gar::sdp::CaloDeposit const& b) const;

#endif

    private:

      int    fTrackID;   ///< g4 track ID of particle making the deposit
      double fTime;      ///< time of the energy deposit
      float  fEnergy;    ///< energy deposited
      float  fX;         ///< x position of the energy deposit
      float  fY;         ///< y position of the energy deposit
      float  fZ;         ///< z position of the energy deposit
      int  fCaloID;        ///< the size of the step for this energy deposit
      unsigned long long int fCellID;
      unsigned int fLayer;
    };


  } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h */