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
    : fEngine(engine)
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
      
      return;
    }
    
    //-------------------------------------------------------------
    void GArAction::BeginOfEventAction(const G4Event*)
    {
      // Clear any previous information.
      fDeposits.clear();
    }
    
    //-------------------------------------------------------------
    void GArAction::PreTrackingAction(const G4Track* /*track*/)
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
      if(step->GetTotalEnergyDeposit() / CLHEP::GeV > fEnergyCut){
        
        // save this deposition
        this->AddEnergyDeposition(step);
        
      } // end if enough energy to worry about this step
      
    }// end of GArAction::SteppingAction

    //------------------------------------------------------------------------------
    void GArAction::AddEnergyDeposition(const G4Step* step)
    {
      // get the track id for this step
      auto trackID = step->GetTrack()->GetTrackID();
      
      // the step mid point
      auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
                             step->GetPostStepPoint()->GetPosition() );
      
      // now figure out the energy deposit
      gar::sdp::EnergyDeposit dep(step->GetPreStepPoint()->GetGlobalTime(),
                                  step->GetTotalEnergyDeposit() / CLHEP::GeV,
                                  midPoint.x() / CLHEP::cm,
                                  midPoint.y() / CLHEP::cm,
                                  midPoint.z() / CLHEP::cm,
                                  midPoint.mag() / CLHEP::cm);
      
      // try inserting a new EnergyDeposits into the fDeposits set, the return
      // of the attemp is a pair whose first element is an iterator either to the
      // new EnergyDeposits object or to the one that was already there
      gar::sdp::EnergyDeposits edeps(trackID);
      
      std::set<gar::sdp::EnergyDeposits>::iterator itr = fDeposits.find(edeps);

      // if we already have an object with that track ID, copy it and then
      // erase the existing one. We are doing this because the EnergyDeposits
      // object is the key in the set, so the set doesn't allow us to mess with
      // it.  However, we can copy it out, erase the set member, and the do
      // a new insert after adding our current energy deposition no problem.
      if( itr != fDeposits.end() ){
        edeps = *itr;
        fDeposits.erase(itr);
      }

      edeps.AddEnergyDeposit(dep);
      fDeposits.insert(edeps);
      
      return;
    }
    
    //------------------------------------------------------------------------------
    void GArAction::EndOfEventAction(const G4Event*)
    {
      
      // sort the EnergyDeposit lists in each EnergyDeposits object
      for(auto deps : fDeposits) deps.Sort();
      
    }
    
    
  } // garg4
  
} // gar