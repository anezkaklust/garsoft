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
      double tpos2[3] = {0.};
      
      // If it's a null step, don't use it.
      if(tpos0[0] == tpos1[0] &&
         tpos0[1] == tpos1[1] &&
         tpos0[2] == tpos1[2]  ) return;
      
      //check that we are in the correct material to record a hit
      std::string material = track->GetMaterial()->GetName();
      if(material.compare("GAr") != 0 ) {
        return;
      }
      
      double dcos[3]  = {mom.x()/mom.mag(),
                         mom.y()/mom.mag(),
                         mom.z()/mom.mag()};
      
      /// Get density for Birks' Law calculation
      TGeoManager *geomanager = fGeo->ROOTGeoManager();
      if ( ! geomanager ) {
        throw cet::exception("NoTGeoManager")
        << "GArAction: no TGeoManager given by fGeo.  "
        << "I need this to get the density.  Failing.\n"
        << __FILE__ << ":" << __LINE__ << "\n";
        return;
      }
      TGeoMaterial *mat = geomanager->GetMaterial("GAr");
      if ( ! mat ) {
        throw cet::exception("NoTGeoMaterial")
        << "GArAction: no TGeoMaterial names 'GAr' found.  "
        << "I need this to get the density.  Failing.\n"
        << __FILE__ << ":" << __LINE__ << "\n";
        return;
      }
      
      /// Getting Energy depositions
      const double edep = step->GetTotalEnergyDeposit()/CLHEP::GeV;
      
      
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