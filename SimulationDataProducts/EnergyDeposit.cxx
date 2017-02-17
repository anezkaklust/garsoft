//
//  EnergyDeposit.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 1/30/17.
//

#include "SimulationDataProducts/EnergyDeposit.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace sdp{
    
    //-------------------------------------------------
    EnergyDeposit::EnergyDeposit()
    : fTime   (std::numeric_limits<double>::max())
    , fEnergy (std::numeric_limits<float >::max())
    , fX      (std::numeric_limits<float >::max())
    , fY		  (std::numeric_limits<float >::max())
    , fZ		  (std::numeric_limits<float >::max())
    , fdX		  (std::numeric_limits<float >::max())
    {}
    
    //-------------------------------------------------
    // order the energy deposits by id, then time, then z position
    bool gar::sdp::EnergyDeposit::operator<(gar::sdp::EnergyDeposit const& b) const
    {
      if( std::abs(fTrackID) <  std::abs(b.TrackID()) )
        return true;
      else if( std::abs(fTrackID) == std::abs(b.TrackID()) ){
        if( fTime < b.Time() )
          return true;
        else if(fTime == b.Time() ){
          if(fZ < b.Z()) return true;
        }
      }

      return false;
    }
    
  } //sdp
} // gar