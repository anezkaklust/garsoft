////////////////////////////////////////////////////////////////////////
//
//  \file DetectorProperties_service.cc
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// GArSoft includes
#include "DetectorInfo/DetectorPropertiesServiceStandardGAr.h"
//#include "DetectorInfo/ECALProperties.h"
#include "DetectorInfo/GArProperties.h"
#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/ServicePack.h" // gar::extractProviders()
#include "DetectorInfo/GArPropertiesService.h"
//#include "DetectorInfo/ECALPropertiesService.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"

#include "TFile.h"

// Art includes
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace detinfo{

    //--------------------------------------------------------------------
    DetectorPropertiesServiceStandardGAr::DetectorPropertiesServiceStandardGAr(fhicl::ParameterSet const& pset,
                                                                         ::art::ActivityRegistry    &reg)
    : fInheritNumberTimeSamples(pset.get<bool>("InheritNumberTimeSamples", false))
    {
      // Register for callbacks.

      reg.sPostOpenFile.watch    (this, &DetectorPropertiesServiceStandardGAr::postOpenFile);
      reg.sPreProcessEvent.watch (this, &DetectorPropertiesServiceStandardGAr::preProcessEvent);
      fProp = std::make_unique<detinfo::DetectorPropertiesStandard>(pset,
                                                                    gar::extractProviders<geo::GeometryGAr,
                                                                                          detinfo::GArPropertiesService,
                                                                                          //detinfo::ECALPropertiesService,
                                                                                          detinfo::DetectorClocksServiceGAr>(),
                                                                    std::set<std::string>({ "InheritNumberTimeSamples" })
                                                                    );

      // at this point we need and expect the provider to be fully configured
      fProp->CheckIfConfigured();

      // Save the parameter set.
      fPS = pset;

    }

    //--------------------------------------------------------------------
    void DetectorPropertiesServiceStandardGAr::reconfigure(fhicl::ParameterSet const& p)
    {
      fProp->ValidateAndConfigure(p, { "InheritNumberTimeSamples" });

      // Save the parameter set.
      fPS = p;

      return;
    }

    //-------------------------------------------------------------
    void DetectorPropertiesServiceStandardGAr::preProcessEvent(const ::art::Event& , art::ScheduleContext)
    {
      // Make sure TPC Clock is updated with TimeService (though in principle it shouldn't change
      fProp->UpdateClocks(gar::providerFrom<detinfo::DetectorClocksServiceGAr>());
    }


    //--------------------------------------------------------------------
    //  Callback called after input file is opened.
    void DetectorPropertiesServiceStandardGAr::postOpenFile(const std::string& filename)
    {
      // Use this method to figure out whether to inherit configuration
      // parameters from previous jobs.

      // There is no way currently to correlate parameter sets saved in
      // sqlite RootFileDB with process history (from MetaData tree).
      // Therefore, we use the approach of scanning every historical
      // parameter set in RootFileDB, and finding all parameter sets
      // that appear to be DetectorPropertiesService configurations.  If all
      // historical parameter sets are in agreement about the value of
      // an inherited parameter, then we accept the historical value,
      // print a message, and override the configuration parameter.  In
      // cases where the historical configurations are not in agreement
      // about the value of an inherited parameter, we ignore any
      // historical parameter values that are the same as the current
      // configured value of the parameter (that is, we resolve the
      // conflict in favor of parameters values that are different than
      // the current configuration).  If two or more historical values
      // differ from the current configuration, throw an exception.
      // Note that it is possible to give precendence to the current
      // configuration by disabling inheritance for that configuration
      // parameter.

      // Don't do anything if no parameters are supposed to be inherited.

      if(!fInheritNumberTimeSamples) return;

      // The only way to access art service metadata from the input file
      // is to open it as a separate TFile object.  Do that now.

      if(filename.size() != 0) {

        TFile* file = (TFile*) new TFile(filename.c_str(), "READ");
        if(file != 0 && !file->IsZombie() && file->IsOpen()) {

          // Open the sqlite datatabase.

          ::art::SQLite3Wrapper sqliteDB(file, "RootFileDB");

          // Loop over all stored ParameterSets.

          unsigned int iNumberTimeSamples = 0;  // Combined value of NumberTimeSamples.
          unsigned int nNumberTimeSamples = 0;  // Number of NumberTimeSamples parameters seen.

          sqlite3_stmt * stmt = 0;
          sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, NULL);
          while (sqlite3_step(stmt) == SQLITE_ROW) {
            //fhicl::ParameterSet ps;
            //fhicl::make_ParameterSet(reinterpret_cast<char const *>(sqlite3_column_text(stmt, 0)), ps);
	    auto ps = fhicl::ParameterSet::make(reinterpret_cast<char const *>(sqlite3_column_text(stmt, 0)));

            // Is this a DetectorPropertiesService parameter set?

            bool psok = isDetectorPropertiesServiceStandardGAr(ps);
            if(psok) {

              // Check NumberTimeSamples

              //	    if(fInheritNumberTimeSamples) {
              unsigned int newNumberTimeSamples = ps.get<unsigned int>("NumberTimeSamples");

              // Ignore parameter values that match the current configuration.

              if(newNumberTimeSamples != fPS.get<unsigned int>("NumberTimeSamples")) {
                if(nNumberTimeSamples == 0)
                  iNumberTimeSamples = newNumberTimeSamples;
                else if(newNumberTimeSamples != iNumberTimeSamples) {
                  throw cet::exception(__FUNCTION__)
                  << "Historical values of NumberTimeSamples do not agree: "
                  << iNumberTimeSamples << " " << newNumberTimeSamples << "\n" ;
                }
                ++nNumberTimeSamples;
                //	    }
              }
            }
          }

          // Done looping over parameter sets.
          // Now decide which parameters we will actually override.

          if(// fInheritNumberTimeSamples &&
             nNumberTimeSamples != 0 &&
             iNumberTimeSamples != fProp->NumberTimeSamples()) {
            MF_LOG_INFO("DetectorPropertiesServiceStandardGAr")
            << "Overriding configuration parameter NumberTimeSamples using historical value.\n"
            << "  Configured value:        " << fProp->NumberTimeSamples() << "\n"
            << "  Historical (used) value: " << iNumberTimeSamples;

            fProp->SetNumberTimeSamples(iNumberTimeSamples);
          }
        }

        // Close file.
        if(file != 0) {
          if(file->IsOpen())
            file->Close();
          delete file;
        }
      }

    }

    //--------------------------------------------------------------------
    //  Determine whether a parameter set is a DetectorPropertiesService configuration.

    bool DetectorPropertiesServiceStandardGAr::isDetectorPropertiesServiceStandardGAr(const fhicl::ParameterSet& ps) const
    {
      // This method uses heuristics to determine whether the parameter
      // set passed as argument is a DetectorPropertiesService configuration
      // parameter set.

      return
      (ps.get<std::string>("service_type", "") == "DetectorPropertiesService")
      && (ps.get<std::string>("service_provider", "") == "DetectorPropertiesServiceStandardGAr")
      ;
#if 0
      // old heuristics here:
      std::string s;
      double d;
      int i;
      unsigned int u;

      bool result = !ps.get_if_present("module_label", s);
      result = result && ps.get_if_present("TriggerOffset", i);
      result = result && ps.get_if_present("SamplingRate", d);
      result = result && ps.get_if_present("NumberTimeSamples", u);

      return result;
#endif // 0
    }

  } // namespace detinfo

} // gar

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::DetectorPropertiesServiceStandardGAr,
                                  gar::detinfo::DetectorPropertiesService)
