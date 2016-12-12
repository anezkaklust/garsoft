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
      
      fChannelMap = std::make_shared<geo::ChannelMapStandardAlg>(sortingParameters);
      
      if(fChannelMap) geom->ApplyChannelMap(fChannelMap);
    
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
