//
//  EnergyDeposit.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 1/30/17.
//

#include "SimulationDataProducts/EnergyDeposit.h"

#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/sim.h"

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
    // order the energy deposits by time.  If two deposits happen at the same
    // time, order them by z position in the detector
    bool gar::sdp::EnergyDeposit::operator<(gar::sdp::EnergyDeposit const& b) const
    {
      if( fTrackID < b.TrackID() ) return true;
      
      if( fTrackID == b.TrackID() ){
        if( fTime < b.Time() ) return true;
      }

      return false;
    }
    
  } //sdp
} // gar