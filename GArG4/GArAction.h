//
//  GArAction.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/12/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//
/// This class implements the G4Base::UserAction interface in order to
/// accumulate a list of hits in the gaseous argon modeled by Geant4.
//

#ifndef GArAction_hpp
#define GArAction_hpp

#include <stdio.h>

// NuTools includes
#include "nutools/G4Base/UserAction.h"

// GArSoft includes
#include "Geometry/Geometry.h"
#include "SimulationDataProducts/SimChannel.h"

//ART includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Forward declarations.
class G4Event;
class G4Track;
class G4Step;
class G4EnergyLossForExtrapolator;

namespace gar {
  namespace garg4 {
    
    ///list of energy deposits from Geant4
    class GArAction : public g4b::UserAction {
      
    public:
      
      // Standard constructors and destructors;
      explicit GArAction(fhicl::ParameterSet const& pset);
      virtual ~GArAction();
      
      void reconfigure(fhicl::ParameterSet const& pset);
      
      // UserActions method that we'll override, to obtain access to
      // Geant4's particle tracks and trajectories.
      void BeginOfEventAction(const G4Event*);
      void EndOfEventAction  (const G4Event*);
      void PreTrackingAction (const G4Track*);
      void PostTrackingAction(const G4Track*);
      void SteppingAction    (const G4Step*);
      
      //  Returns the FLSHitList accumulated during the current event.
      std::vector<gar::sdp::IDE> const& GArIDEList() const { return fIDEList; }
      
    private:
      std::vector<gar::sdp::IDE>        fIDEList;    ///< The accumulated information for hits in the event.
      double                            fEnergyCut;  ///< The minimum energy in GeV for a particle to
                                                     ///< be included in the list.
      art::ServiceHandle<geo::Geometry> fGeo;        ///< handle to geometry service
      
    };
    
  } // garg4
} // gar

#endif /* GArAction_hpp */
