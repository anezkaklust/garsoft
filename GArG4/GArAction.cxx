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

namespace gar {

  namespace garg4 {
    
    //-------------------------------------------------------------
    // Constructor.
    GArAction::GArAction(fhicl::ParameterSet const& pset)
    : fDriftAlg(nullptr)
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
    void GArAction::reconfigure(fhicl::ParameterSet const& pset )
    {
      fEnergyCut = pset.get<double>("EnergyCut") * CLHEP::GeV;
      
      auto driftAlgName = pset.get<std::string>("ElectronDriftAlg", "Standard");
      
      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::garg4::ElectronDriftStandardAlg>();
      else
        throw cet::exception("GArAction")
        << "Unable to determine which electron drift algorithm to use, bail"

      return;
    }
    
    //-------------------------------------------------------------
    void GArAction::BeginOfEventAction(const G4Event*)
    {
      // Clear any previous particle information.
      fIDEList.clear();
    }
    
    //-------------------------------------------------------------
    // Create our initial sim::FLSHit object and add it to the sim::FLSHitList.
    void GArAction::PreTrackingAction(const G4Track* track)
    {
    }
    
    //-------------------------------------------------------------
    void GArAction::PostTrackingAction( const G4Track* /*track*/)
    {
    }
    
    //-------------------------------------------------------------
    // With every step, add to the particle's trajectory.
    void GArAction::SteppingAction(const G4Step* step)
    {
      //mf::LogInfo("GArAction") << "GArAction::SteppingAction";
      
      // Get the pointer to the track
      G4Track *track = step->GetTrack();
      
      const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector &stop  = track->GetPosition();
      const CLHEP::Hep3Vector &mom   = track->GetMomentum();
      
      LOG_DEBUG("GArAction")
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
      std::string material = track->GetMaterial()->GetName();
      if(material.compare("GAr") != 0 ) {
        return;
      }
      
      
      // Getting Energy depositions
      const double edep = step->GetTotalEnergyDeposit()/CLHEP::GeV;
      
      // only worry about non-zero energy depositions.
      if(edep > 0){
        // reset the IonizationAndScintillation singleton
        gar::garg4::IonizationAndScintillation::Instance()->Reset(step);
        
        auto ides = fDriftAlg->DriftElectronsToReadout(step);
        
      
      }
      
//        fFLSHit->AddEdep( edep );
//
//        // Get the position and time the particle first enters
//          // the volume, as well as the pdg code.  that information is
//          // obtained from the G4Track object that corresponds to this step
//          // account for the fact that we use cm, ns, GeV rather than the G4 defaults
//        fFLSHit->SetPDG    ( track->GetDefinition()->GetPDGEncoding() );
//        fFLSHit->SetTrackId( ParticleListAction::GetCurrentTrackID()  );
      
      
    }// end of GArAction::SteppingAction
    
    
    //------------------------------------------------------------------------------
    void GArAction::EndOfEventAction(const G4Event*)
    {
      
    }
    
    
  } // garg4
  
} // gar