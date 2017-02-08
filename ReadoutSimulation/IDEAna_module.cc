////////////////////////////////////////////////////////////////////////
//
// \file IDEAna_module.cc
//
// \author dmckee@phys.ksu.edu
//
// \brief Build some histograms based on the som::SimChannels created
// \brief by LArVoxelReadout
////////////////////////////////////////////////////////////////////////
// C++ std library includes
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>
// POSIX includes
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// Root Includes
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1.h"
#include "TVector3.h"
// GArSoft includes
#include "Geometry/Geometry.h"
#include "SimulationDataProducts/SimChannel.h"

namespace gar {
  namespace geo { class Geometry; }
  
  namespace rosim {
    
    /// Base class for creation of raw signals on wires.
    class IDEAna : public ::art::EDAnalyzer {
      
    public:
      
      explicit IDEAna(fhicl::ParameterSet const& pset);
      virtual ~IDEAna();
      
      /// read/write access to event
      void analyze (const ::art::Event& evt);
      void beginJob();
      void endJob();
      void reconfigure(fhicl::ParameterSet const& p);
      
    private:
      
      std::string fIonizationModuleLabel;
      
      TH1D * fChargeXpos;     ///< position of the MC Truth charge deposition
      TH1D * fChargeYpos;     ///< position of the MC Truth charge deposition
      TH1D * fChargeZpos;     ///< position of the MC Truth charge deposition
      
      TH1D * fTDC;            ///< Which TDCs have activity
      
      TH1D * fTDCsPerChannel; ///< Number of TDCs with activity
      TH1D * fIDEsPerChannel;
      
      TH1D * fElectrons;      ///< Electrons in the whole channel entry
      TH1D * fEnergy;         ///< Energy in the whole channel entry
      
      TH1D * fElectronsPerTDC;
      TH1D * fEnergyPerTDC;
      
      TH1D * fElectronsPerIDE;
      TH1D * fEnergyPerIDE;
      
    }; // class IDEAna
    
    
    //-------------------------------------------------
    IDEAna::IDEAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fChargeXpos()
    , fChargeYpos()
    , fChargeZpos()
    , fTDC()
    , fTDCsPerChannel()
    , fIDEsPerChannel()
    , fElectrons()
    , fEnergy()
    , fElectronsPerTDC()
    , fEnergyPerTDC()
    , fElectronsPerIDE()
    , fEnergyPerIDE()
    {
      this->reconfigure(pset);
    }
    
    //-------------------------------------------------
    IDEAna::~IDEAna()
    {
    }
    
    void IDEAna::reconfigure(fhicl::ParameterSet const& p)
    {
      fIonizationModuleLabel = p.get< std::string >("IonizationModuleLabel");
      return;
    }
    
    //-------------------------------------------------
    void IDEAna::beginJob()
    {
      // get access to the TFile service
      ::art::ServiceHandle<::art::TFileService> tfs;

      int    timeSamples = 10000.;
      
      fTDC             = tfs->make<TH1D>("hTDC",
                                         "Active TDC;TDCs;Events;",
                                         timeSamples, 0, 1. * timeSamples);
      fTDCsPerChannel  = tfs->make<TH1D>("hTDCsPerChannel",
                                         "TDCs per channel entry;# TDCs;Events",
                                         128, 0, 1. * timeSamples);
      fIDEsPerChannel  = tfs->make<TH1D>("hIDEsPerChannel",
                                         "IDE per channel entry;# IDEs;Events",
                                         100, 0, 20000);
      fElectrons       = tfs->make<TH1D>("hElectrons",
                                         "Electrons per channel;Electrons;Events",
                                         100, 0, 2e7);
      fElectronsPerIDE = tfs->make<TH1D>("hElectronsPerIDE",
                                         "Electrons per IDE;Electrons;Events",
                                         100, 0, 10000);
      fElectronsPerTDC = tfs->make<TH1D>("hElectronsPerTDC",
                                         "Electrons per TDC;Electrons;Events",
                                         100, 0, 10000);
      return;
      
    }
    
    //-------------------------------------------------
    void IDEAna::endJob() {}
    
    //-------------------------------------------------
    void IDEAna::analyze(const ::art::Event& evt)
    {
      
      if (evt.isRealData()) {
        throw cet::exception("IDEAna")
        << "Not for use on Data\n";
      }
      
      ::art::ServiceHandle<geo::Geometry> geom;
      
      auto ideCol = evt.getValidHandle<std::vector<sdp::IDE> >(fIonizationModuleLabel);
      
      //++++++++++
      // Loop over the Chnnels and fill histograms
      //++++++++++
      unsigned int totalIDEs      = ideCol->size();
      double       totalElectrons = 0;
      double       tdcElectrons   = 0;
      
      for(const auto& ide : *ideCol){

        fTDC->Fill(ide.TDC());
          
        totalElectrons += ide.NumElectrons();
        tdcElectrons   += ide.NumElectrons();
        
        fElectronsPerIDE->Fill(ide.NumElectrons());
        
      } // end loop over IDEs
      
      fTDCsPerChannel->Fill(std::numeric_limits<float>::max());
      fIDEsPerChannel->Fill(totalIDEs);
      fElectrons->Fill(totalElectrons);
      fElectronsPerTDC->Fill(tdcElectrons);

      return;
    }//end analyze method
    
    DEFINE_ART_MODULE(IDEAna)
    
  } // end of rosim namespace
} // gar

