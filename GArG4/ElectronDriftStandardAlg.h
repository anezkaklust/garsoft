//
//  ElectronDriftStandardAlg.h
//
//  Default implementation of an electron drift algorithm
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef ElectronDriftStandardAlg_h
#define ElectronDriftStandardAlg_h

#include "GArG4/ElectronDriftAlg.h"

namespace gar {
  namespace garg4{
    
    class ElectronDriftStandardAlg : public ElectronDriftAlg {
      
    public:
      
      ElectronDriftStandardAlg(CLHEP::HepRandomEngine      & engine,
                               fhicl::ParameterSet    const& pset);
      virtual ~ElectronDriftStandardAlg();
      
      std::vector<sdp::TDCIDE> DriftElectronsToReadout(const G4Step* step);
      
    private:
      
      int fElectronsPerCluster;     ///< Number of electrons to drift in a cluster
      int fMinimumElectronClusters; ///< Minimum number of clusters to drift
      
    };
    
  }
} // gar



#endif /* ElectronDriftStandardAlg_h */
