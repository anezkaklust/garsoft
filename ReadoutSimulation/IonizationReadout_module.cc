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
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {
    
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
      
      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      void DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                   std::vector<edepIDE>                 & edepIDEs);
      void CombineIDEs(std::vector<edepIDE> & edepIDEs);
      
      std::string                         fG4Label;    ///< label of G4 module
      std::unique_ptr<ElectronDriftAlg>   fDriftAlg;   ///< algorithm to drift ionization electrons
      const gar::detinfo::DetectorClocks* fTime;       ///< electronics clock
      std::unique_ptr<TPCReadoutSimAlg>   fROSimAlg;   ///< algorithm to simulate the electronics
      fhicl::ParameterSet                 fISCalcPars; ///< parameter set for the IS calculator
      
    };
    
  } // namespace rosim
  
  namespace rosim {
    
    //----------------------------------------------------------------------
    // Constructor
    IonizationReadout::IonizationReadout(fhicl::ParameterSet const& pset)
    {
      fTime  = gar::providerFrom<detinfo::DetectorClocksService>();
      
      // setup the random number service
      // obtain the random seed from NuRandomService
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "ionization", pset, "IonizationSeed");
      
      this->reconfigure(pset);
      
      produces< std::vector<raw::RawDigit>                      >();
      produces< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit> >();
      
      return;
    }
    
    //----------------------------------------------------------------------
    // Destructor
    IonizationReadout::~IonizationReadout()
    {
    }

    //----------------------------------------------------------------------
    void IonizationReadout::reconfigure(fhicl::ParameterSet const& pset)
    {
      LOG_DEBUG("IonizationReadout") << "Debug: IonizationReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;
      
      fISCalcPars = pset.get<fhicl::ParameterSet>("ISCalcPars"            );
      fG4Label    = pset.get<std::string        >("G4ModuleLabel", "geant");
      
      auto driftAlgPars = pset.get<fhicl::ParameterSet>("ElectronDriftAlgPars");
      auto driftAlgName = driftAlgPars.get<std::string>("DriftAlgType");
      
      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::rosim::ElectronDriftStandardAlg>(rng->getEngine("ionization"),
                                                                           driftAlgPars);
      else
        throw cet::exception("IonizationReadout")
        << "Unable to determine which electron drift algorithm to use, bail";
      
      auto tpcROAlgPars = pset.get<fhicl::ParameterSet>("TPCReadoutSimAlgPars");
      auto tpcROAlgName = tpcROAlgPars.get<std::string>("TPCReadoutSimType");
      
      if(tpcROAlgName.compare("Standard") == 0)
        fROSimAlg = std::make_unique<gar::rosim::TPCReadoutSimStandardAlg>(rng->getEngine("ionization"),
                                                                           tpcROAlgPars);
      else
        throw cet::exception("IonizationReadout")
        << "Unable to determine which TPC readout simulation algorithm to use, bail";
      
      return;
    }
    
    //----------------------------------------------------------------------
    void IonizationReadout::beginJob()
    {
      auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());
      
      // create the ionization and scintillation calculator;
      // this is a singleton (!) so we just need to make the instance in one
      // location
      IonizationAndScintillation::CreateInstance(rng->getEngine("ionization"),
                                                 fISCalcPars);
      
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
      std::unique_ptr< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit> > erassn(new ::art::Assns<sdp::EnergyDeposit, raw::RawDigit>);
      
      // first get the energy deposits from the event record
      auto eDepCol = evt.getValidHandle< std::vector<sdp::EnergyDeposit> >(fG4Label);
      std::vector<edepIDE> eDepIDEs;

      // drift the ionization electrons to the readout and create edepIDE objects
      this->DriftElectronsToReadout(*eDepCol, eDepIDEs);
      
      // The IDEs should have been combined already so that there are no repeat
      // TDC values for any channel
      auto               numTicks = gar::providerFrom<detinfo::DetectorPropertiesService>()->NumberTimeSamples();
      unsigned int       prevChan = eDepIDEs.front().Channel;
      std::set<size_t>   digitEDepLocs;
      std::vector<float> electrons(numTicks, 0.);

      // make the signal raw digits and set their associations to the energy deposits
      for(auto edide : eDepIDEs){
        
        if(edide.Channel != prevChan){
          rdCol->emplace_back(fROSimAlg->CreateRawDigit(edide.Channel, electrons));
        
          // loop over the locations in the eDepCol to make the associations
          for(auto ed : digitEDepLocs){
            auto const ptr = art::Ptr<sdp::EnergyDeposit>(eDepCol, ed);
            util::CreateAssn(*this, evt, *rdCol, ptr, *erassn);
          }

          digitEDepLocs.clear();
          electrons.clear();
          electrons.resize(numTicks, 0);
          prevChan = edide.Channel;
        }

        electrons[edide.TDC] = edide.NumElect;
        
        for(auto loc : edide.edepLocs) digitEDepLocs.insert(loc);
        
      } // end loop to fill signal raw digit vector and make EnergyDeposit associations

      // now make the noise digits
      fROSimAlg->CreateNoiseDigits(*rdCol);
      
      evt.put(std::move(rdCol));
      evt.put(std::move(erassn));
      
      return;
    } // IonizationReadout::produce()
    
    //--------------------------------------------------------------------------
    void IonizationReadout::DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                                    std::vector<edepIDE>                 & edepIDEs)
    {
      ::art::ServiceHandle<gar::geo::Geometry> geo;
      
      float        xyz[3] = {0.};
      unsigned int chan   = 0;

      // now instantiate an ElectronDriftInfo object to keep track of the
      // drifted locations of each electron cluster from each energy deposit
      rosim::ElectronDriftInfo driftInfo;

      // loop over the energy deposits
      for(size_t e = 0; e < edepCol.size(); ++e){
        
        // get the positions, arrival times, and electron cluster sizes
        // for this energy deposition
        fDriftAlg->DriftElectronsToReadout(edepCol[e], driftInfo);
        
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
          
          // see if this cluster of electrons can be mapped to a channel.
          // if not, move on to the next one
          try{
            chan = geo->NearestChannel(xyz);
          }
          catch(cet::exception &e){
            continue;
          }
          
          edepIDEs.emplace_back(clusterSize[c],
                                chan,
                                fTime->TPCG4Time2TDC(clusterTime[c]),
                                e);
          
        }

      } // end loop over deposit collections

      this->CombineIDEs(edepIDEs);
      
      return;
    }
    
    //--------------------------------------------------------------------------
    void IonizationReadout::CombineIDEs(std::vector<edepIDE> & edepIDEs)
    {
      std::vector<edepIDE> temp;
      
      // sort the edepIDE objects.  This is the sorting by IDEs because that is
      // all the < operator of edepIDEs does
      std::sort(edepIDEs.begin(), edepIDEs.end());
      
      unsigned int   chan     = 0;
      unsigned int   prevChan = edepIDEs.front().Channel;
      unsigned short tdc      = 0;
      unsigned short prevTDC  = edepIDEs.front().TDC;
      
      edepIDE sum(edepIDEs.front());
      
      // now loop over the sorted vector and combine any edepIDEs
      // with the same channel and tdc values for the IDEs
      for(auto edepIDE : edepIDEs){
      
        chan = edepIDE.Channel;
        tdc  = edepIDE.TDC;
        
        if(chan != prevChan ||
           tdc  != prevTDC){
          
          // put the summed edepIDE into the temp vector
          temp.push_back(sum);
          
          // start over with a fresh sum
          sum      = edepIDE;
          prevChan = chan;
          prevTDC  = tdc;
          
        }
        else
          sum += edepIDE;
        
      } // end loop to sum edepIDEs with the same channel and tdc values
      
      // now swap the input vector with the temp vector
      temp.swap(edepIDEs);
      
      return;
    }
    
  } // namespace rosim
  
  namespace rosim {
    
    DEFINE_ART_MODULE(IonizationReadout)
    
  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_IONIZATIONREADOUT
