////////////////////////////////////////////////////////////////////////
/// \file  GArG4.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GARG4GArG4ANA_H
#define GARG4GArG4ANA_H 

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// LArSoft Includes
#include "MCCheater/BackTracker.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "Geometry/Geometry.h"

// ROOT includes
#include <TTree.h>

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

typedef struct{
  int run;
  int subrun;
  int event;
} EventInfo;

typedef struct{
  int    trackID;
  int    pdg;
  float  x;
  float  y;
  float  z;
  float  e;
  float  dX;
  float  t;
} EnergyDep;

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

///Geant4 interface
namespace gar {
  namespace garg4 {
    
    class GArG4Ana : public ::art::EDAnalyzer{
    public:
      
      /// Standard constructor and destructor for an FMWK module.
      explicit GArG4Ana(fhicl::ParameterSet const& pset);
      virtual ~GArG4Ana();
      
      void analyze (const ::art::Event& evt);
      void beginJob();
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      std::string fG4ModuleLabel;    ///< module label for the Geant
      std::string fTruthModuleLabel; ///< module label for the Geant
      
      TTree*      fTree;             ///< Tree to keep track of sanity check info
      EventInfo   fEvt;              ///< Struct containing event identification
      EnergyDep   fEDep;             ///< Struct containing energy deposition info
    };
    
  } // namespace garg4
  
  namespace garg4 {
    
    //-----------------------------------------------------------------------
    // Constructor
    GArG4Ana::GArG4Ana(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    {
      this->reconfigure(pset);
    }
    
    //-----------------------------------------------------------------------
    // Destructor
    GArG4Ana::~GArG4Ana()
    {
    }
    
    //-----------------------------------------------------------------------
    void GArG4Ana::beginJob()
    {
      ::art::ServiceHandle<::art::TFileService> tfs;
      ::art::ServiceHandle<geo::Geometry> geo;
      
      fTree    = tfs->make<TTree>("G4Tree", "G4Tree");

      std::string description("run/I:subrun/I:event/I");
      
      fTree->Branch("info", &fEvt, description.c_str());
      
      description = "trackID/I:pdg/I:x/F:y/F:z/F:e/F:dX/F:t/F";
      
      fTree->Branch("edep", &fEDep, description.c_str());
      
      
    }
    
    //-----------------------------------------------------------------------
    void GArG4Ana::reconfigure(fhicl::ParameterSet const& p)
    {
      fG4ModuleLabel    = p.get< std::string >("GeantModuleLabel");
      fTruthModuleLabel = p.get< std::string >("TruthModuleLabel");
      
      return;
    }
    
    //-----------------------------------------------------------------------
    void GArG4Ana::analyze(const ::art::Event& evt)
    {
      // set the event id info
      fEvt.run = evt.run();
      fEvt.subrun = evt.subRun();
      fEvt.event  = evt.id().event();
      
      // get the list of particles from this event
      ::art::ServiceHandle<cheat::BackTracker> bt;

      // get the energy depositions for this event
      auto edepsCol = evt.getValidHandle<std::vector<gar::sdp::EnergyDeposit> >(fG4ModuleLabel);
      
      if( edepsCol.failedToGet() ){
        LOG_VERBATIM("GArG4Ana")
        << "failed to get energy depositions, bail";
        return;
      }

      // loop over the energy deposition collections to fill the tree
      for(auto edep : *edepsCol){
        
        // get the MCParticle for this track ID
        auto part = bt->TrackIDToParticle(edep.TrackID());
        
        if( !part ) continue;
        
        fEDep.trackID = edep.TrackID();
        fEDep.pdg     = part->PdgCode();
        fEDep.x       = edep.X();
        fEDep.y       = edep.Y();
        fEDep.z       = edep.Z();
        fEDep.dX      = edep.dX();
        fEDep.t       = edep.Time();
        fEDep.e       = edep.Energy();
        
        LOG_VERBATIM("GArG4Ana")
        << "pos: ("
        << fEDep.x
        << ", "
        << fEDep.y
        << ", "
        << fEDep.z
        << ") e: "
        << fEDep.e
        << " t: "
        << fEDep.t
        << " dX: "
        << fEDep.dX;
        
          // make the tree flat in terms of the energy depositions
        fTree->Fill();
        
      } // end loop over collection of EnergyDeposits
      
      return;
    } // analyze
    
  } // namespace garg4
  
  namespace garg4 {
    
    DEFINE_ART_MODULE(GArG4Ana)
    
  } // namespace garg4
  
} // gar

#endif // GARG4GARG4H

