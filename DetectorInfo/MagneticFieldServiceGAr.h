#ifndef MAGNETICFIELDSERVICEGAR_H
#define MAGNETICFIELDSERVICEGAR_H

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"
#include "DetectorInfo/GArMagneticField.h"

namespace mag {
  class MagneticFieldServiceGAr: public MagneticFieldService {
  public:
    MagneticFieldServiceGAr(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

  private:
    void reconfigure(fhicl::ParameterSet const& pset);
    void preBeginRun(const art::Run& run);

    const provider_type* provider() const override { return &fProp; }

    mag::GArMagneticField fProp;
  };
}

DECLARE_ART_SERVICE_INTERFACE_IMPL(mag::MagneticFieldServiceGAr,
                                   mag::MagneticFieldService,
                                   SHARED)

#endif