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
#include "Geometry/GeoObjectSorterStandard.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar{
  namespace geo{
    
    class ChannelMapStandardAlg : public ChannelMapAlg{
      
    public:
      
      ChannelMapStandardAlg(fhicl::ParameterSet const& p);
      
      void         Initialize( GeometryData_t& geodata ) override;
      void         Uninitialize();
      unsigned int NChannels() const;
      
    private:
      
      unsigned int                  fNChannels;      ///< number of channels in the detector
      
    };
    
    
  } // end namespace geo
} // end namespace gar
#endif // GEO_CHANNELMAPSTANDARDALG_H

