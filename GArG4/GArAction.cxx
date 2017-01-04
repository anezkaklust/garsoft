//
//  GArAction.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/12/16.
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

#include "GArG4/GArAction.h"
#include "GArG4/ParticleListAction.h"
#include "GArG4/IonizationAndScintillation.h"
#include "GArG4/ElectronDriftStandardAlg.h"

namespace gar {

  namespace garg4 {
    
    //-------------------------------------------------------------
    // Constructor.
    GArAction::GArAction(CLHEP::HepRandomEngine*    engine,
                         fhicl::ParameterSet const& pset)
    : fDriftAlg(nullptr)
    , fEngine(engine)
    {
      this->reconfigure(pset);
    }
    
    //-------------------------------------------------------------
    // Destructor.
    GArAction::~GArAction()
    {
      // Delete anything that we created with "new'.
    }
    
    //-------------------------------------------------------------
    void GArAction::reconfigure(fhicl::ParameterSet const& pset)
    {
      fEnergyCut = pset.get<double>("EnergyCut") * CLHEP::GeV;
      
      auto driftAlgName = pset.get<std::string>("ElectronDriftAlg", "Standard");
      
      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::garg4::ElectronDriftStandardAlg>(*fEngine, pset);
      else
        throw cet::exception("GArAction")
        << "Unable to determine which electron drift algorithm to use, bail";

      return;
    }
    
    //-------------------------------------------------------------
    void GArAction::BeginOfEventAction(const G4Event*)
    {
      // Clear any previous information.
      fDriftAlg->Reset();
    }
    
    //-------------------------------------------------------------
    void GArAction::PreTrackingAction(const G4Track* track)
    {
    }
    
    //-------------------------------------------------------------
    void GArAction::PostTrackingAction( const G4Track* /*track*/)
    {
    }
    
    //-------------------------------------------------------------
    // With every step, handle energy deposition in gaseous argon
    void GArAction::SteppingAction(const G4Step* step)
    {
      LOG_DEBUG("GArAction")
      << "GArAction::SteppingAction";
      
      // Get the pointer to the track
      G4Track *track = step->GetTrack();
      
      const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector &stop  = track->GetPosition();
      
      // If it's a null step, don't use it.
      if(start[0] == stop[0] &&
         start[1] == stop[1] &&
         start[2] == stop[2]  ) return;
      
      // check that we are in the correct material to record a hit
      std::string material = track->GetMaterial()->GetName();
      if(material.compare("GAr") != 0 ) {
        return;
      }
      
      // only worry about energy depositions larger than the minimum required
      if(step->GetTotalEnergyDeposit() > fEnergyCut){
        
        // drift the ionization electrons to the readout
        fDriftAlg->DriftElectronsToReadout(step);
        
        // here is where we would also call an algorithm to figure out how much
        // scintillation light was produced
        
      } // end if enough energy to worry about this step
      
    }// end of GArAction::SteppingAction
    
    
    //------------------------------------------------------------------------------
    void GArAction::EndOfEventAction(const G4Event*)
    {
      
    }
    
    
  } // garg4
  
} // gar