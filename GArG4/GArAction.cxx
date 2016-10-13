//
//  GArAction.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/12/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "GArG4/GArAction.h"

namespace gar {

  namespace garg4 {
    
    //-------------------------------------------------------------
    // Constructor.
    GArAction::GArAction(fhicl::ParameterSet const& pset)
    {
      this->reconfigure();
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
      
        /// Get the pointer to the track
      G4Track *track = step->GetTrack();
      
      
      const CLHEP::Hep3Vector &pos0 = step->GetPreStepPoint()->GetPosition(); // Start of the step
      const CLHEP::Hep3Vector &pos  = track->GetPosition();                   // End of the step
      const CLHEP::Hep3Vector &mom  = track->GetMomentum();
      
      LOG_DEBUG("GArAction") << "flshits momentum = "
				  << mom.x() << " " << mom.y() << " " << mom.z()
				  << " " << mom.mag();
      
      double tpos0[3] = {pos0.x()/CLHEP::cm, pos0.y()/CLHEP::cm, pos0.z()/CLHEP::cm}; ///< Start of the step
      double tpos1[3] = {pos.x()/CLHEP::cm , pos.y()/CLHEP::cm , pos.z()/CLHEP::cm};  ///< End of the step
      double tpos2[3] = {0.};                                    ///< Some temporary vector that we'll make use of
      
      // If it's a null step, don't use it. Otherwise it may induce an additional
      // FLSHit, which is wrong
      if(tpos0[0] == tpos1[0] &&
         tpos0[1] == tpos1[1] &&
         tpos0[2] == tpos1[2]  ) return;
      
      //check that we are in the correct material to record a hit - ie scintillator
      std::string material = track->GetMaterial()->GetName();
      if(material.compare("Scintillator") != 0 ) {
        return;
      }
      
      //check if the current cell pointer is empty, that means
      //this is the first step for this track
      double dcos[3]  = {mom.x()/mom.mag(), mom.y()/mom.mag(), mom.z()/mom.mag()};
      
      /// Get density for Birks' Law calculation
      TGeoManager *geomanager = fGeo->ROOTGeoManager();
      if ( ! geomanager ) {
        throw cet::exception("NoTGeoManager")
        << "GArAction: no TGeoManager given by fGeo.  "
        << "I need this to get the scintillator density.  Failing.\n"
        << __FILE__ << ":" << __LINE__ << "\n";
        return;
      }
      TGeoMaterial *mat = geomanager->GetMaterial("Scintillator");
      if ( ! mat ) {
        throw cet::exception("NoTGeoMaterial")
        << "GArAction: no TGeoMaterial names 'Scintillator' found.  "
        << "I need this to get the scintillator density.  Failing.\n"
        << __FILE__ << ":" << __LINE__ << "\n";
        return;
      }
      
      /// Getting Energy depositions
      const double edep      = step->GetTotalEnergyDeposit()/CLHEP::GeV;
      fNovaG4EmSaturation.SetChouConstant(fChouConstant);
      const double edepBirks = fNovaG4EmSaturation.VisibleEnergyDeposition(step)/CLHEP::GeV;
      
      
        fFLSHit->AddEdep( edep );

        // Get the position and time the particle first enters
          // the volume, as well as the pdg code.  that information is
          // obtained from the G4Track object that corresponds to this step
          // account for the fact that we use cm, ns, GeV rather than the G4 defaults
        fFLSHit->SetPDG    ( track->GetDefinition()->GetPDGEncoding() );
        fFLSHit->SetTrackId( ParticleListAction::GetCurrentTrackID()  );
        
      
    }// end of GArAction::SteppingAction
    
    
    //------------------------------------------------------------------------------
    void GArAction::EndOfEventAction(const G4Event*)
    {
      
    }
    
    
  } // garg4
  
} // gar