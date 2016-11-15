////////////////////////////////////////////////////////////////////////////////
/// \file StandardGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "Geometry/StandardGeometryHelper.h"

// GArSoft libraries
#include "Geometry/ChannelMapStandardAlg.h"
#include "Geometry/GeometryCore.h"

// C/C++ libraries
#include <string>

namespace gar
{
  namespace geo
  {
    
    //----------------------------------------------------------------------------
    StandardGeometryHelper::StandardGeometryHelper(fhicl::ParameterSet const& pset,
                                                   ::art::ActivityRegistry    &)
    : fPset( pset )
    , fChannelMap()
    {}
    
    //----------------------------------------------------------------------------
    void StandardGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                                                          geo::GeometryCore        * geom)
    {
      fChannelMap.reset();
      
      std::string const detectorName = geom->DetectorName();
      
      // Migration note:
      // Should just create ChannelMapStandardAlg with no decision-making after transition
      // detector names in this code must be all lower case
      if((detectorName.find("argoneut")   != std::string::npos) ||
         (detectorName.find("microboone") != std::string::npos) ||
         (detectorName.find("bo")         != std::string::npos)
        ){
        fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>(sortingParameters);
      }
      
      if ( fChannelMap ){
        geom->ApplyChannelMap(fChannelMap);
      }
    
      return;
      
    } // StandardGeometryHelper::doConfigureChannelMapAlg()
    
    
    //----------------------------------------------------------------------------
    StandardGeometryHelper::ChannelMapAlgPtr_t StandardGeometryHelper::doGetChannelMapAlg() const
    {
      return fChannelMap;
    }
    
    //----------------------------------------------------------------------------
    
  } // namespace geo
} // namespace gar

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::geo::StandardGeometryHelper, gar::geo::ExptGeoHelperInterface)
