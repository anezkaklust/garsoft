//
//  RawDigitAna_module.cc
//
//  Created by Brian Rebel on 3/2/17.
//

#ifndef GAR_ROSIM_RAWDIGITANA
#define GAR_ROSIM_RAWDIGITANA

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// GArSoft Includes
#include "MCCheater/BackTracker.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "Geometry/Geometry.h"
#include "RawDataProducts/RawDigit.h"

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
  int   trackID;
  int   pdg;
  float x;
  float y;
  float z;
  float e;
  float dX;
  float t;
} EnergyDep;

typedef struct{
  int   channel;
  float x;
  float y;
  float z;
} ChannelInfo;

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

namespace gar {
  namespace rosim {
    
    class RawDigitAna : public ::art::EDAnalyzer{
    public:
      
        /// Standard constructor and destructor for an FMWK module.
      explicit RawDigitAna(fhicl::ParameterSet const& pset);
      virtual ~RawDigitAna();
      
      void analyze (const ::art::Event& evt);
      void beginJob();
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      std::string  fReadoutModuleLabel; ///< module label for the Geant
      std::string  fG4ModuleLabel;      ///< module label for the Geant
      
      TTree*       fEDepTree;          ///< Tree to keep track of sanity check info
      EventInfo    fEvt;               ///< Struct containing event identification
      EnergyDep    fEDep;              ///< Struct containing energy deposition info
      ChannelInfo  fChannelInfo;       ///< Struct containing channel info
    };
    
  } // namespace rosim
  
  namespace rosim {
    
    //-----------------------------------------------------------------------
    // Constructor
    RawDigitAna::RawDigitAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    {
      this->reconfigure(pset);
    }
    
    //-----------------------------------------------------------------------
    // Destructor
    RawDigitAna::~RawDigitAna()
    {
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::beginJob()
    {
      ::art::ServiceHandle<::art::TFileService> tfs;
      ::art::ServiceHandle<geo::Geometry> geo;
      
      fEDepTree = tfs->make<TTree>("EDepTree", "EDepTree");
      
      std::string description("run/I:subrun/I:event/I");
      
      fEDepTree->Branch("info", &fEvt, description.c_str());
      
      description = "trackID/I:pdg/I:x/F:y/F:z/F:e/F:dX/F:t/F";
      
      fEDepTree->Branch("edep", &fEDep, description.c_str());
      
      description = "channel/I:x/F:y/F:z/F";
      
      fEDepTree->Branch("chan", &fChannelInfo, description.c_str());
      
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::reconfigure(fhicl::ParameterSet const& p)
    {
      fG4ModuleLabel      = p.get< std::string >("GeantModuleLabel"  );
      fReadoutModuleLabel = p.get< std::string >("ReadoutModuleLabel");
      
      return;
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::analyze(const ::art::Event& evt)
    {
      // set the event id info
      fEvt.run    = evt.run();
      fEvt.subrun = evt.subRun();
      fEvt.event  = evt.id().event();
      
      // get the list of particles from this event
      ::art::ServiceHandle<cheat::BackTracker> bt;
      
      // get the raw digits for this event
      auto digCol = evt.getValidHandle<std::vector<gar::raw::RawDigit> >(fReadoutModuleLabel);
      
      if( digCol.failedToGet() ){
        LOG_VERBATIM("RawDigitAna")
        << "failed to get raw digits, bail";
        return;
      }
 
      // get the find many for the energy depositions
      ::art::FindMany<gar::sdp::EnergyDeposit> fmEnergyDep(digCol, evt, fReadoutModuleLabel);
      
      if( !fmEnergyDep.isValid() ){
        LOG_WARNING("RawDigitAna")
        << "Unable to find valid association between RawDigits and "
        << "energy deposits, no analysis in this module is possible";
        
        return;
      }

      auto geo = gar::providerFrom<geo::Geometry>();

      std::vector<const gar::sdp::EnergyDeposit*> edepsCol;

      float xyz[3] = {0.};
      
      for(size_t d = 0; d < digCol->size(); ++d){
        fmEnergyDep.get(d, edepsCol);
        if(edepsCol.size() < 1) continue;
        
        // fill the channel information
        fChannelInfo.channel = (*digCol)[d].Channel();

        geo->ChannelToPosition((*digCol)[d].Channel(), xyz);
        
        fChannelInfo.x = xyz[0];
        fChannelInfo.y = xyz[1];
        fChannelInfo.z = xyz[2];

        // loop over the energy deposition collections to fill the tree
        for(auto edep : edepsCol){
        
          // get the MCParticle for this track ID
          auto part = bt->TrackIDToParticle(edep->TrackID());
          
          fEDep.trackID = edep->TrackID();
          fEDep.pdg     = part->PdgCode();
          fEDep.x       = edep->X();
          fEDep.y       = edep->Y();
          fEDep.z       = edep->Z();
          fEDep.dX      = edep->dX();
          fEDep.t       = edep->Time();
          fEDep.e       = edep->Energy();
          
          LOG_DEBUG("RawDigitAna")
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
          fEDepTree->Fill();
          
        } // end loop over collection of EnergyDeposits
      } // end loop over RawDigit collection
      
      return;
    } // analyze
    
  } // namespace rosim
  
  namespace rosim {
    
    DEFINE_ART_MODULE(RawDigitAna)
    
  } // namespace rosim
  
} // gar

#endif // GAR_ROSIM_RAWDIGITANA

