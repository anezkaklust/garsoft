////////////////////////////////////////////////////////////////////////
//
// \file GArSimChannelAna_module.cc
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
// LArSoft includes
#include "Geometry/Geometry.h"
#include "SimulationDataProducts/SimChannel.h"

namespace gar {
  namespace geo { class Geometry; }
  
  namespace garg4 {
    
    /// Base class for creation of raw signals on wires.
    class GArSimChannelAna : public ::art::EDAnalyzer {
      
    public:
      
      explicit GArSimChannelAna(fhicl::ParameterSet const& pset);
      virtual ~GArSimChannelAna();
      
      /// read/write access to event
      void analyze (const ::art::Event& evt);
      void beginJob();
      void endJob();
      void reconfigure(fhicl::ParameterSet const& p);
      
    private:
      
      std::string fGArG4ModuleLabel;
      
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
      
    }; // class GArSimChannelAna
    
    
    //-------------------------------------------------
    GArSimChannelAna::GArSimChannelAna(fhicl::ParameterSet const& pset)
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
    GArSimChannelAna::~GArSimChannelAna()
    {
    }
    
    void GArSimChannelAna::reconfigure(fhicl::ParameterSet const& p)
    {
      fGArG4ModuleLabel = p.get< std::string >("GArGeantModuleLabel");
      return;
    }
    
    //-------------------------------------------------
    void GArSimChannelAna::beginJob()
    {
      // get access to the TFile service
      ::art::ServiceHandle<::art::TFileService> tfs;

      // geometry data.
      ::art::ServiceHandle<gar::geo::Geometry> geom;
      
      // get the dimensions of the detector
      double width       = geom->DetHalfWidth() * 2.;
      double halfHeight  = geom->DetHalfHeight();
      double length      = geom->DetLength();
      int    timeSamples = 10000.;
      
      fChargeXpos      = tfs->make<TH1D>("hChargeXpos",
                                         "X charge depositions;X (cm);Events",
                                         101, 0.0, width);
      fChargeYpos      = tfs->make<TH1D>("hChargeYpos",
                                         "Y charge depositions;Y (cm);Events",
                                         101, -halfHeight, halfHeight);
      fChargeZpos      = tfs->make<TH1D>("hChargeZpos",
                                         "Z charge depositions;Z (cm);Events",
                                         101, 0.0, length);
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
      fEnergy          = tfs->make<TH1D>("hEnergy",
                                         "Energy per channel;energy;Events",
                                         100, 0, 2500);
      fElectronsPerIDE = tfs->make<TH1D>("hElectronsPerIDE",
                                         "Electrons per IDE;Electrons;Events",
                                         100, 0, 10000);
      fEnergyPerIDE    = tfs->make<TH1D>("hEnergyPerIDE",
                                         "Energy per IDE;energy;Events",
                                         100, 0, 50);
      fElectronsPerTDC = tfs->make<TH1D>("hElectronsPerTDC",
                                         "Electrons per TDC;Electrons;Events",
                                         100, 0, 10000);
      fEnergyPerTDC    = tfs->make<TH1D>("hEnergyPerTDC",
                                         "Energy per YDC;energy;Events",
                                         100, 0, 50);
      return;
      
    }
    
    //-------------------------------------------------
    void GArSimChannelAna::endJob() {}
    
    //-------------------------------------------------
    void GArSimChannelAna::analyze(const ::art::Event& evt)
    {
      
      if (evt.isRealData()) {
        throw cet::exception("GArSimChannelAna") << "Not for use on Data yet...\n";
      }
      
      ::art::ServiceHandle<geo::Geometry> geom;
      
      auto chanHandle = evt.getValidHandle<std::vector<sdp::SimChannel> >(fGArG4ModuleLabel);
      std::vector<sdp::SimChannel> const& scVec(*chanHandle);
      
      //++++++++++
      // Loop over the Chnnels and fill histograms
      //++++++++++
      unsigned int totalIDEs      = 0;
      double       totalElectrons = 0;
      double       totalEnergy    = 0;
      
      for(const auto& sc : scVec){

        const auto & tdcidemap = sc.TDCIDEs();
        fTDCsPerChannel->Fill(tdcidemap.size());
        
        for(auto const& tdcide : tdcidemap){
          unsigned int tdc          = tdcide.fTDC;
          auto const&  ideVec       = tdcide.fIDEs;
          double       tdcElectrons = 0.;
          double       tdcEnergy    = 0.;

          totalIDEs += ideVec.size();
          
          fTDC->Fill(tdc);
          
          for(auto const& ide : ideVec) {
            totalElectrons += ide.numElectrons;
            totalEnergy    += ide.energy;
            tdcElectrons   += ide.numElectrons;
            tdcEnergy      += ide.energy;
            
            fChargeXpos->Fill(ide.x);
            fChargeYpos->Fill(ide.y);
            fChargeZpos->Fill(ide.z);
            fElectronsPerIDE->Fill(ide.numElectrons);
            fEnergyPerIDE->Fill(ide.energy);
          }
          fElectronsPerTDC->Fill(tdcElectrons);
          fEnergyPerTDC->Fill(tdcEnergy);
          
        } // end loop over TDCIDEs
      } // end loop over sim channels
      
      fIDEsPerChannel->Fill(totalIDEs);
      fElectrons->Fill(totalElectrons);
      fEnergy->Fill(totalEnergy);

      return;
    }//end analyze method
    
    DEFINE_ART_MODULE(GArSimChannelAna)
    
  } // end of garg4 namespace
} // gar

