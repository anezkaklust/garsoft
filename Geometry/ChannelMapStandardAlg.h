////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELSTANDARDMAPALG_H
#define GEO_CHANNELSTANDARDMAPALG_H

#include <vector>
#include <set>
#include <iostream>

#include "Geometry/ChannelMapAlg.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar{
  namespace geo{
    
    class ChannelMapStandardAlg : public ChannelMapAlg{
      
    public:
      
      ChannelMapStandardAlg(fhicl::ParameterSet const& p);
      
      void         Initialize()                             override;
      void         Uninitialize()                           override;
      unsigned int Nchannels()                        const override;
      unsigned int NearestChannel(float const xyz[3]) const override;
      
    private:
      
      unsigned int fNChannels; ///< number of channels in the detector
      
    };
    
    
  } // end namespace geo
} // end namespace gar
#endif // GEO_CHANNELMAPSTANDARDALG_H

