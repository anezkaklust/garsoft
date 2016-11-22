//
//  AuxDetAction.cxx
//
//  Created by Brian Rebel on 11/22/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"

  // G4 includes
#include "Geant4/G4Event.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PrimaryParticle.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4VUserPrimaryParticleInformation.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4VProcess.hh"

#include "GArG4/AuxDetAction.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar {
  
  namespace garg4 {
    
    //-------------------------------------------------------------
    // Constructor.
    AuxDetAction::AuxDetAction(CLHEP::HepRandomEngine*    engine,
                               fhicl::ParameterSet const& pset)
    : fEngine(engine)
    {
      fGeo = gar::providerFrom<geo::Geometry>();

      this->reconfigure(pset);
    }
    
    //-------------------------------------------------------------
    // Destructor.
    AuxDetAction::~AuxDetAction()
    {
      // Delete anything that we created with "new'.
    }
    
    //-------------------------------------------------------------
    void AuxDetAction::reconfigure(fhicl::ParameterSet const& pset)
    {
      fEnergyCut = pset.get<double>("EnergyCut") * CLHEP::GeV;
      
      return;
    }
    
    //-------------------------------------------------------------
    void AuxDetAction::BeginOfEventAction(const G4Event*)
    {
      // Clear any previous information.
      fAuxDetSimChannels.clear();
    }
    
    //-------------------------------------------------------------
    void AuxDetAction::PreTrackingAction(const G4Track* track)
    {
    }
    
    //-------------------------------------------------------------
    void AuxDetAction::PostTrackingAction( const G4Track* /*track*/)
    {
    }
    
    //-------------------------------------------------------------
    // With every step, check and see if we are in an auxdet and save the energy
    void AuxDetAction::SteppingAction(const G4Step* step)
    {
      LOG_DEBUG("AuxDetAction")
      << "AuxDetAction::SteppingAction";
      
      // Get the pointer to the track
      G4Track *track = step->GetTrack();
      
      const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector &stop  = track->GetPosition();
      const CLHEP::Hep3Vector &mom   = track->GetMomentum();
      
      LOG_DEBUG("AuxDetAction")
      << "step momentum = "
      << mom.x()
      << " "
      << mom.y()
      << " "
      << mom.z()
      << " "
      << mom.mag();
      
      double tpos0[3] = {start.x()/CLHEP::cm, start.y()/CLHEP::cm, start.z()/CLHEP::cm};
      double tpos1[3] = {stop .x()/CLHEP::cm, stop .y()/CLHEP::cm, stop .z()/CLHEP::cm};
      
      // If it's a null step, don't use it.
      if(tpos0[0] == tpos1[0] &&
         tpos0[1] == tpos1[1] &&
         tpos0[2] == tpos1[2]  ) return;
      
      // check that we are in the correct material to record a hit
      std::string volname = track->GetVolume()->GetName();
      if(volname.find("AuxDet") == std::string::npos ) {
        return;
      }
      
      // only worry about energy depositions larger than the minimum required
      if(step->GetTotalEnergyDeposit() > fEnergyCut){
        
        // look for the nearest AuxDetChannel to the energy deposition
        
        // find that object in the fAuxDetSimChannels collection
        
        //  Add this step to the collected IDEs for the AuxDetSimChannel 
        
      } // end if enough energy to worry about this step
      
    }// end of AuxDetAction::SteppingAction
    
    
    //------------------------------------------------------------------------------
    void AuxDetAction::EndOfEventAction(const G4Event*)
    {
      
    }
    
    
  } // garg4
  
} // gar