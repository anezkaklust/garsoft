//
//  EnergyDeposit.h
//
//  Created by Brian Rebel on 1/30/17.
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h
#define GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h

#include <list>

namespace gar {
  namespace sdp {
    
    
    class EnergyDeposit{
    public:
      
      EnergyDeposit();
      
#ifndef __GCCXML__
      EnergyDeposit(double t,
                    float  e,
                    float  xpos,
                    float  ypos,
                    float  zpos,
                    float  length)
      : fTime   (t)
      , fEnergy (e)
      , fX      (xpos)
      , fY      (ypos)
      , fZ      (zpos)
      , fdX     (length)
      {}
      
      double const& Time()   const { return fTime;   }
      float  const& Energy() const { return fEnergy; }
      float  const& X()      const { return fX;      }
      float  const& Y()      const { return fY;      }
      float  const& Z()      const { return fZ;      }
      float  const& dX()     const { return fdX;     }
      
      bool operator  <(gar::sdp::EnergyDeposit const& b) const;
      
#endif
      
    private:
      
      double fTime;    ///< time of the energy deposit
      float  fEnergy;  ///< energy deposited
      float  fX;       ///< x position of the energy deposit
      float  fY;       ///< y position of the energy deposit
      float  fZ;       ///< z position of the energy deposit
      float  fdX;      ///< the size of the step for this energy deposit
    };
    

    class EnergyDeposits{
      
    public:
      
      EnergyDeposits();
      
#ifndef __GCCXML__
      // Constructor specifying the track ID
      EnergyDeposits(int trackID);
      
      // Copy constructor with an offset from the G4 ID
      EnergyDeposits(EnergyDeposits const& deps,
                     int                   offset);
      
      void AddEnergyDeposit(EnergyDeposit const& dep);
      
      int                      const& TrackID()  const { return fTrackID;          }
      std::list<EnergyDeposit> const& Deposits() const { return fDeposits;         }
      void                            Sort()           { fDeposits.sort(); return; }
      
      bool operator  <(gar::sdp::EnergyDeposits const& b) const;
      
#endif
      
    private:
      
      int                      fTrackID;  ///< the G4 track ID for the particle making the deposits
      std::list<EnergyDeposit> fDeposits; ///< the collection of energy deposits, ordered by time
      
      
    };
  } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h */
