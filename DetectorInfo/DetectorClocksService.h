////////////////////////////////////////////////////////////////////////
// DetectorClocksService.h
//
// Pure virtual service interface for DetectorClocks functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORCLOCKSSERVICE_H
#define DETECTORCLOCKSSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "DetectorInfo/DetectorClocks.h"
#include "CoreUtils/ServiceUtil.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    class DetectorClocksService {
      
    public:
      typedef detinfo::DetectorClocks provider_type;
      
    public:
      virtual ~DetectorClocksService() = default;
      
      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const  detinfo::DetectorClocks* provider() const = 0;
      
    }; // class DetectorClocksService
  } //namespace detinfo
} //gar

DECLARE_ART_SERVICE_INTERFACE(gar::detinfo::DetectorClocksService, LEGACY)

#endif // DETECTORCLOCKSSERVICE_H
