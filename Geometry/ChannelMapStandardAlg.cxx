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
    
    
  } // namespace
} // namespace gar