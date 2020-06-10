////////////////////////////////////////////////////////////////////////////////
/// \file StandardGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "Geometry/StandardGeometryHelper.h"

// GArSoft libraries
#include "Geometry/ChannelMapAlgs/ChannelMapStandardAlg.h"
#include "Geometry/ChannelMapAlgs/ECALSegmentationGridXYAlg.h"
#include "Geometry/ChannelMapAlgs/ECALSegmentationMultiGridStripXYAlg.h"
#include "Geometry/ChannelMapAlgs/MinervaSegmentationAlg.h"
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
    , fECALSegmentationAlg()
    {}

    //----------------------------------------------------------------------------
    void StandardGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                                                          geo::GeometryCore        * geom)
    {
      fChannelMap.reset();

      fChannelMap = std::make_shared<gar::geo::seg::ChannelMapStandardAlg>(sortingParameters);

      if(fChannelMap) geom->ApplyChannelMap(fChannelMap);

      return;

    } // StandardGeometryHelper::doConfigureChannelMapAlg()


    //----------------------------------------------------------------------------
    StandardGeometryHelper::ChannelMapAlgPtr_t StandardGeometryHelper::doGetChannelMapAlg() const
    {
      return fChannelMap;
    }

    //----------------------------------------------------------------------------
    void StandardGeometryHelper::doConfigureECALSegmentationAlg(fhicl::ParameterSet const& segParameters,
                                                          geo::GeometryCore        * geom)
    {
      fECALSegmentationAlg.reset();

      auto SegmentationAlgName = segParameters.get<std::string>("SegmentationAlgName");

      if(SegmentationAlgName.compare("GridXY") == 0)
      fECALSegmentationAlg = std::make_shared<gar::geo::seg::ECALSegmentationGridXYAlg>(segParameters);
      else if(SegmentationAlgName.compare("MultiGridStripXY") == 0)
      fECALSegmentationAlg = std::make_shared<gar::geo::seg::ECALSegmentationMultiGridStripXYAlg>(segParameters);
      else if(SegmentationAlgName.compare("Minerva") == 0)
      fECALSegmentationAlg = std::make_shared<gar::geo::seg::MinervaSegmentationAlg>(segParameters);
      else{
          throw cet::exception("StandardGeometryHelper::doConfigureECALSegmentationAlg")
          << "Unable to determine which ECAL Segmentation algorithm to use, bail";
      }

      if(fECALSegmentationAlg && SegmentationAlgName.compare("Minerva") != 0) geom->ApplyECALSegmentationAlg(fECALSegmentationAlg);
      //Case for Minerva segmentation algorithm
      if(fECALSegmentationAlg && SegmentationAlgName.compare("Minerva") == 0) geom->ApplyMinervaSegmentationAlg(fECALSegmentationAlg);

      return;

    } // StandardGeometryHelper::doConfigureChannelMapAlg()


    //----------------------------------------------------------------------------
    StandardGeometryHelper::ECALSegmentationAlgPtr_t StandardGeometryHelper::doGetECALSegmentationAlg() const
    {
      return fECALSegmentationAlg;
    }

    //----------------------------------------------------------------------------

  } // namespace geo
} // namespace gar

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::geo::StandardGeometryHelper, gar::geo::ExptGeoHelperInterface)
