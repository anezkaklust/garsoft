//
//  ElectronDriftAlg.h
//
//  Algorithm base class for drifting electrons in an electric field
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef ElectronDriftAlg_h
#define ElectronDriftAlg_h

#include <vector>

#include "CLHEP/Random/RandGauss.h"

class G4Step;

namespace fhicl {
  class ParameterSet;
}

namespace gar {
  
  namespace sdp{
    class EnergyDeposit;
  }
  
  namespace rosim{
    
    class ElectronDriftAlg{
      
    public:
      
      ElectronDriftAlg(CLHEP::HepRandomEngine      & engine,
                       fhicl::ParameterSet    const& pset);
      virtual ~ElectronDriftAlg();
      
      // The user must call IonizationAndScintillation::Instance()->Reset(dep)
      // from within the implementation of this method in the derived class
      virtual void DriftElectronsToReadout(gar::sdp::EnergyDeposit const& dep) = 0;
      
    protected:
      
      double                         fDriftVelocity;         ///< electron drift velocity
      double                         fInverseVelocity;       ///< stored for computational convenience
      double                         fLifetimeCorrection;    ///< electron lifetime correction in negative ms
      double                         fLongitudinalDiffusion; ///< diffusion along the drift
      double                         fLongDiffConst;         ///< stored for computational convenience
      double                         fTransverseDiffusion;   ///< diffusion transverse to the drift
      double                         fTransDiffConst;        ///< stored for computational convenience
      double                         fFanoFactor;            ///< Fano factor
      CLHEP::HepRandomEngine&        fEngine;                ///< random number engine
    };
    
  }
} // gar

#endif /* ElectronDriftAlg_h */
