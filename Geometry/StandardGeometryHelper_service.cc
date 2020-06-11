////////////////////////////////////////////////////////////////////////////////
/// \file StandardGeometryHelper_service.cc
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// class header
#include "Geometry/StandardGeometryHelper.h"

// GArSoft libraries
#include "Geometry/ChannelMapAlgs/ChannelMapStandardAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationGridXYAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationMultiGridStripXYAlg.h"
#include "Geometry/ChannelMapAlgs/MinervaSegmentationAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationStripXAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationStripYAlg.h"
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
    , fMinervaSegmentationAlg()
    , fMuIDSegmentationAlg()
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
      fECALSegmentationAlg = std::make_shared<gar::geo::seg::SegmentationGridXYAlg>(segParameters);
      else if(SegmentationAlgName.compare("MultiGridStripXY") == 0)
      fECALSegmentationAlg = std::make_shared<gar::geo::seg::SegmentationMultiGridStripXYAlg>(segParameters);
      else{
          throw cet::exception("StandardGeometryHelper::doConfigureECALSegmentationAlg")
          << "Unable to determine which ECAL Segmentation algorithm to use, bail";
      }

      if(fECALSegmentationAlg) geom->ApplyECALSegmentationAlg(fECALSegmentationAlg);

      return;

    } // StandardGeometryHelper::doConfigureECALSegmentationAlg()

    //----------------------------------------------------------------------------
    void StandardGeometryHelper::doConfigureMinervaSegmentationAlg(fhicl::ParameterSet const& segParameters,
                                                          geo::GeometryCore        * geom)
    {
      fMinervaSegmentationAlg.reset();

      auto SegmentationAlgName = segParameters.get<std::string>("SegmentationAlgName");

      if(SegmentationAlgName.compare("Minerva") == 0)
      fMinervaSegmentationAlg = std::make_shared<gar::geo::seg::MinervaSegmentationAlg>(segParameters);
      else{
          throw cet::exception("StandardGeometryHelper::doConfigureMinervaSegmentationAlg")
          << "Unable to determine which Tracker Sc Segmentation algorithm to use, bail";
      }

      if(fMinervaSegmentationAlg) geom->ApplyMinervaSegmentationAlg(fMinervaSegmentationAlg);

      return;

    } // StandardGeometryHelper::doConfigureMinervaSegmentationAlg()

    //----------------------------------------------------------------------------
    void StandardGeometryHelper::doConfigureMuIDSegmentationAlg(fhicl::ParameterSet const& segParameters,
                                                          geo::GeometryCore        * geom)
    {
      fMuIDSegmentationAlg.reset();

      auto SegmentationAlgName = segParameters.get<std::string>("SegmentationAlgName");

      if(SegmentationAlgName.compare("SegmentationStripXAlg") == 0)
      fMuIDSegmentationAlg = std::make_shared<gar::geo::seg::SegmentationStripXAlg>(segParameters);
      else if(SegmentationAlgName.compare("SegmentationStripYAlg") == 0)
      fMuIDSegmentationAlg = std::make_shared<gar::geo::seg::SegmentationStripYAlg>(segParameters);
      else{
          throw cet::exception("StandardGeometryHelper::doConfigureMuIDSegmentationAlg")
          << "Unable to determine which MuID Segmentation algorithm to use, bail";
      }

      if(fMuIDSegmentationAlg) geom->ApplyMuIDSegmentationAlg(fMuIDSegmentationAlg);

      return;

    } // StandardGeometryHelper::doConfigureMuIDSegmentationAlg()

    //----------------------------------------------------------------------------
    StandardGeometryHelper::SegmentationAlgPtr_t StandardGeometryHelper::doGetECALSegmentationAlg() const
    {
      return fECALSegmentationAlg;
    }

    //----------------------------------------------------------------------------
    StandardGeometryHelper::SegmentationAlgPtr_t StandardGeometryHelper::doGetMinervaSegmentationAlg() const
    {
      return fMinervaSegmentationAlg;
    }

    //----------------------------------------------------------------------------
    StandardGeometryHelper::SegmentationAlgPtr_t StandardGeometryHelper::doGetMuIDSegmentationAlg() const
    {
      return fMuIDSegmentationAlg;
    }

    //----------------------------------------------------------------------------

  } // namespace geo
} // namespace gar

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::geo::StandardGeometryHelper, gar::geo::ExptGeoHelperInterface)
