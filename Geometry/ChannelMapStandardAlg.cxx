////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapStandardAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace gar {
  namespace geo{
    
    //----------------------------------------------------------------------------
    ChannelMapStandardAlg::ChannelMapStandardAlg(fhicl::ParameterSet const& p)
    : fNChannels(std::numeric_limits<unsigned int>::max())
    {
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Initialize()
    {
      // start over:
      Uninitialize();
      
      return;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Uninitialize()
    {
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::Nchannels() const
    {
      return fNChannels;
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::NearestChannel(float const xyz[3]) const
    {
      unsigned int channel = std::numeric_limits<unsigned int>::max();
      
      // insert code here to calculate the nearest channel based on the
      // supplied position.  If the position is outside the bounds for channels
      // throw an exception.
      
      return channel;
    }
    
  } // namespace
} // namespace gar