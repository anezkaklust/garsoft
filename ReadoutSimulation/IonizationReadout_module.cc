////////////////////////////////////////////////////////////////////////
/// \file  IonizationReadout_module.cc
/// \brief Create the readout from the ionization signal
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_IONIZATIONREADOUT
#define GAR_READOUTSIMULATION_IONIZATIONREADOUT

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "DetectorInfo/DetectorClocksService.h"
#include "GArG4/G4SimulationParameters.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ElectronDriftAlg.h"
#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {

    struct edepIDE {
      
      edepIDE();
      edepIDE(int            trackID,
              float          numE,
              unsigned int   chan,
              unsigned short tdc,
              size_t         edepIdx)
      : ide    (gar::sdp::IDE(trackID, numE, chan, tdc))
      {
        edepLocs.push_back(edepIdx);
      }

      edepIDE(gar::sdp::IDE       const& id,
              std::vector<size_t> const& edepIdx)
      : ide(gar::sdp::IDE(id))
      {
        for(auto const e : edepIdx) edepLocs.push_back(e);
      }

      void AddEDep(size_t edepLoc) { edepLocs.push_back(edepLoc); }
      
      bool operator <(edepIDE const& b) const
      {
        return (this->ide < b.ide);
      }

      void operator +=(edepIDE const& b) const
      {
        this->ide += b.ide;
        
        for(auto const e : b.edepLocs)
          edepLocs.push_back(e);

        return;
      }

      gar::sdp::IDE       ide;
      std::vector<size_t> edepLocs;
    };
    
    /**
     * @brief Runs Readout simulation including propagation of electrons and photons to readout
     *
     *
     * Randomness
     * -----------
     *
     * The random number generators used by this process are:
     * - 'GEANT' instance: used by Geant4
     * - 'propagation' instance: used in electron propagation
     * - 'radio' instance: used for radiological decay
     *
     *
     */
    class IonizationReadout : public ::art::EDProducer{
    public:
      
      /// Standard constructor and destructor for an FMWK module.
      explicit IonizationReadout(fhicl::ParameterSet const& pset);
      virtual ~IonizationReadout();
      
      /// The main routine of this module: Fetch the primary particles
      /// from the event, simulate their evolution in the detctor, and
      /// produce the detector response.
      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      
    private:
      
      void DriftElectronsToReadout(std::vector<sdp::EnergyDeposit const& edepCol,
                                   std::vector<edepIDE>                & edepIDEs);
      void CombineIDEs(std::vector<edepIDE> & edepIDEs)
      void IDEToRawDigits(std::vector<sdp::IDE>      const& ides,
                          std::vector<raw::RawDigit>      & rawDigits);
      
      std::string                         fG4Label;  ///< label of G4 module
      std::unique_ptr<ElectronDriftAlg>   fDriftAlg; ///< algorithm to drift ionization electrons
      gar::detinfo::ElecClock             fClock;    ///< electronics clock
      const gar::detinfo::DetectorClocks* fTime;     ///< electronics clock
      std::unique_ptr<TPCReadoutSimAlg>   fROSimAlg; ///< algorithm to simulate the electronics
      
    };
    
  } // namespace rosim
  
  namespace rosim {
    
    //----------------------------------------------------------------------
    // Constructor
    IonizationReadout::IonizationReadout(fhicl::ParameterSet const& pset)
    {
      fTime  = gar::providerFrom<detinfo::DetectorClocksService>();
      fClock = fTime->TPCClock();
      
      // initialize the GArSimulationParameters singleton
      garg4::G4SimulationParameters::CreateInstance(pset.get<fhicl::ParameterSet>("GArSimParsPSet"));

      LOG_DEBUG("IonizationReadout") << "Debug: IonizationReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;
      
      // setup the random number service for Geant4, the "G4Engine" label is a
      // special tag setting up a global engine for use by Geant4/CLHEP;
      // obtain the random seed from NuRandomService
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "ionization", pset, "IonizationSeed");
      
      fG4Label = pset.get<std::string>("G4ModuleLabel", "geant");
      
      auto driftAlgPars = pset.get<fhicl::ParameterSet>("ElectronDriftAlgPars");
      auto driftAlgName = driftAlgPars.get<std::string>("DriftAlgType");
      
      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::rosim::ElectronDriftStandardAlg>(rng->getEngine("ionization"),
                                                                           driftAlgPars);
      else
        throw cet::exception("IonizationReadout")
        << "Unable to determine which electron drift algorithm to use, bail";

      auto tpcROAlgPars = pset.get<fhicl::ParameterSet>("TPCReadoutSimAlgPars");
      auto tpcROAlgName = driftAlgPars.get<std::string>("TPCReadoutSimType");
      
      if(tpcROAlgName.compare("Standard") == 0)
        fROSimAlg = std::make_unique<gar::rosim::TPCReadoutSimStandardAlg>(rng->getEngine("ionization"),
                                                                           tpcROAlgPars);
      else
        throw cet::exception("IonizationReadout")
        << "Unable to determine which TPC readout simulation algorithm to use, bail";

      
      produces< std::vector<raw::RawDigit>                      >();
      produces< std::vector<sdp::IDE>                           >();
      produces< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit> >();
      produces< ::art::Assns<sdp::EnergyDeposit, sdp::IDE>      >();
      
      
      return;
    }
    
    //----------------------------------------------------------------------
    // Destructor
    IonizationReadout::~IonizationReadout()
    {
    }
    
    //----------------------------------------------------------------------
    void IonizationReadout::beginJob()
    {
      auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());
      
      // create the ionization and scintillation calculator;
      // this is a singleton (!) so it does not make sense
      // to create it in LArVoxelReadoutGeometry
      IonizationAndScintillation::CreateInstance(rng->getEngine("ionization"));
      
      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::beginRun(::art::Run& run)
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    void IonizationReadout::produce(::art::Event& evt)
    {
      LOG_DEBUG("IonizationReadout") << "produce()";
      
      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<raw::RawDigit>                      > rdCol (new std::vector<raw::RawDigit>                     );
      std::unique_ptr< std::vector<sdp::IDE>                           > ideCol(new std::vector<sdp::IDE>                          );
      std::unique_ptr< ::art::Assns<raw::RawDigit, sdp::EnergyDeposit> > erassn(new ::art::Assns<raw::RawDigit, sdp::EnergyDeposit>);
      std::unique_ptr< ::art::Assns<sdp::IDE,      sdp::EnergyDeposit> > eiassn(new ::art::Assns<sdp::IDE,      sdp::EnergyDeposit>);
      
      // first get the energy deposits from the event record
      auto eDepCol = evt.getValidHandle< std::vector<sdp::EnergyDeposit> >(fG4Label);
      std::vector<edepIDE> eDepIDEs;

      // drift the ionization electrons to the readout and create sdp::IDE objects
      this->DriftElectronsToReadout(*eDepCol, eDepIDEs);
      
      // make the IDE vector and set their associations to the energy deposits
      for(auto edide : eDepIDEs){
        ideCol->emplace_back(edide.ide);
        
        // loop over the locations in the eDepCol to make the associations
        for(auto ed : edide.edepLocs){
          util::CreateAssn(this, evt, *ideCol, art::Ptr<sdp::EnergyDeposit>(eDepCol, ed), *eiassn);
        }
        
      } // end loop to fill IDE vector and make IDE/EnergyDeposit associations
      
      // now create the RawDigits using the IDEs
      this->IDEToRawDigits(*ideCol, *rdCol);
      
      evt.put(std::move(rdCol));
      evt.put(std::move(ideCol));
      evt.put(std::move(eiassn));
      evt.put(std::move(erassn));
      
      return;
    } // IonizationReadout::produce()
    
    //--------------------------------------------------------------------------
    void IonizationReadout::DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                                    std::vector<edepIDE>                 & edepIDEs)
    {
      ::art::ServiceHandle<gar::geo::Geometry> geo;
      
      float          xyz[3] = {0.};
      float          numEl  = 0.;
      unsigned int   chan   = 0;
      unsigned short tdc    = 0;

      // now instantiate an ElectronDriftInfo object to keep track of the
      // drifted locations of each electron cluster from each energy deposit
      rosim::ElectronDriftInfo driftInfo;

      // loop over the energy deposits
      for(size_t e = 0; e < eDepCol->size(); ++e){
        
        // get the positions, arrival times, and electron cluster sizes
        // for this energy deposition
        fDriftAlg->DriftElectronsToReadout((*eDepCol)[e], driftInfo);
        
        auto clusterXPos = driftInfo.ClusterXPos();
        auto clusterYPos = driftInfo.ClusterYPos();
        auto clusterZPos = driftInfo.ClusterZPos();
        auto clusterTime = driftInfo.ClusterTime();
        auto clusterSize = driftInfo.ClusterSize();
        
        // the vectors should all have the same size by the time we get them
        // here (verified by the ElectronDriftInfo object when they are filled)
        for(size_t c = 0; c < clusterXPos.size(); ++c){
          
          xyz[0] = clusterXPos[c];
          xyz[1] = clusterYPos[c];
          xyz[2] = clusterZPos[c];
          numEl  = clusterSize[c];
          chan   = geo->NearestChannel(xyz);
          tdc    = fClock.Ticks(fTime->G4ToElecTime(clusterTime[c]));
          
          eDepIDEs.emplace_back((*eDepCol)[e].TrackID(), numEl, chan, tdc, e);
          
        }

      } // end loop over deposit collections

      this->CombineIDEs(eDepIDEs);
      
      
      return;
    }
    
    //--------------------------------------------------------------------------
    void IonizationReadout::CombineIDEs(std::vector<edepIDE> & edepIDEs)
    {
      std::vector<edepIDE> temp;
      
      // sort the edepIDE objects.  This is the sorting by IDEs because that is
      // all the < operator of edepIDEs does
      std::sort(eDepIDEs.begin(), eDepIDEs.end());
      
      unsigned int   chan     = 0;
      unsigned int   prevChan = edepIDEs.front().ide.Channel();
      unsigned short tdc      = 0;
      unsigned short prevTDC  = edepIDEs.front().ide.TDC();
      
      edepIDE sum(edepIDEs.front());
      
      // now loop over the sorted vector and combine any edepIDEs
      // with the same channel and tdc values for the IDEs
      for(size_t e = 1; e < edepIDEs.size(); ++e){
      
        chan = edepIDEs[e].ide.Channel();
        tdc  = edepIDEs[e].ide.TDC();
        
        if(chan != prevChan ||
           tdc  != prevTDC){
          
          // put the summed edepIDE into the temp vector
          temp.push_back(sum);
          
          // start over with a fresh sum
          sum      = edepIDE(edepIDEs[e]);
          prevChan = chan;
          prevTDC  = tdc;
          
        }
        else
          sum += edepIDEs[e];
        
      } // end loop to sum edepIDEs with the same channel and tdc values
      
      // now swap the input vector with the temp vector
      temp.swap(edepIDEs);
      
      return;
    }
    
    //--------------------------------------------------------------------------
    void IonizationReadout::IDEToRawDigit(std::vector<sdp::IDE>      const& ides,
                                          std::vector<raw::RawDigit>      & rawDigits)
    {
      auto temp = fROSimAlg->IDEToRawDigits(ids);
      
      temp.swap(rawDigits);
      
      return;
    }
    
  } // namespace rosim
  
  namespace rosim {
    
    DEFINE_ART_MODULE(IonizationReadout)
    
  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_IONIZATIONREADOUT
