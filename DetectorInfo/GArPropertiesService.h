////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef GArPropertiesSERVICE_H
#define GArPropertiesSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "DetectorInfo/GArProperties.h"
#include "CoreUtils/ServiceUtil.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    class GArPropertiesService {
      
    public:
      typedef detinfo::GArProperties provider_type;
      
    public:
      virtual ~GArPropertiesService() = default;
      
      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const  detinfo::GArProperties* provider() const = 0;
      
    }; // class GArPropertiesService
  } //namespace detinfo
  DECLARE_ART_SERVICE_INTERFACE(detinfo::GArPropertiesService, LEGACY)
} // gar
#endif // GArPropertiesSERVICE_H
