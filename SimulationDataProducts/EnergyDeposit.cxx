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
      if(fTime < b.Time() ) return true;
      else if(fTime == b.Time()){
        if(fZ < b.Z()) return true;
      }
      
      return false;
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //-------------------------------------------------
    EnergyDeposits::EnergyDeposits()
    : fTrackID(std::numeric_limits<int>::min())
    {
      
    }

    //-------------------------------------------------
    EnergyDeposits::EnergyDeposits(int trackID)
    : fTrackID(trackID)
    {
      
    }

    //-------------------------------------------------
    EnergyDeposits::EnergyDeposits(EnergyDeposits const& deps,
                                   int                   offset)
    : fTrackID(offset + deps.TrackID())
    {
      for(auto const& dep : deps.Deposits()) fDeposits.push_back(dep);
      
      return;
    }

    //-------------------------------------------------
    void EnergyDeposits::AddEnergyDeposit(EnergyDeposit const& dep)
    {
      
      // if these objects end up being too big, we will have to find a
      // way to consolidate the energy deposit information, ie grouping
      // deposits that are close to each other in time and space.
      
      fDeposits.push_back(dep);
      
      return;
    }
    
    //-------------------------------------------------
    // order the energy deposit collections by track ID
    bool gar::sdp::EnergyDeposits::operator<(gar::sdp::EnergyDeposits const& b) const
    {
      return fTrackID < b.TrackID();
    }

  } //sdp
} // gar