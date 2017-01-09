//
//  ElectronDriftAlg.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/18/16.
//

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "GArG4/ElectronDriftAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar {
  namespace garg4 {
   
    //--------------------------------------------------------------------------
    // The parameter set is in case a derived algorithm needs to be configured
    // beyond what is available from either DetectorPropertiesService or
    // GArPropertiesService.
    ElectronDriftAlg::ElectronDriftAlg(CLHEP::HepRandomEngine      & engine,
                                       fhicl::ParameterSet    const& /*pset*/)
    : fEngine(engine)
    {
      // get the corrections and constants from the necessary places
      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      auto garProp = gar::providerFrom<detinfo::GArPropertiesService>();
      
      fDriftVelocity         = detProp->DriftVelocity(detProp->Efield(),
                                                      detProp->Temperature());
      fLifetimeCorrection    = -1000. * detProp->ElectronLifetime();
      fLongitudinalDiffusion = garProp->LongitudinalDiffusion();
      fTransverseDiffusion   = garProp->TransverseDiffusion();
      fFanoFactor            = garProp->FanoFactor();
      
      fInverseVelocity       = 1. / fDriftVelocity;
      fLongDiffConst         = std::sqrt(2. * fLongitudinalDiffusion);
      fTransDiffConst        = std::sqrt(2. * fTransverseDiffusion);
      
      return;
    }
    
    //--------------------------------------------------------------------------
    ElectronDriftAlg::~ElectronDriftAlg()
    {
    }

    
    
  }
}