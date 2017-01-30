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
#include "GArG4/G4SimulationParameters.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ElectronDriftAlg.h"
#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {
    
    // Forward declarations within namespace.
    
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
      
      std::unique_ptr<ElectronDriftAlg> fDriftAlg;  ///< algorithm to drift ionization electrons
      
      
    };
    
  } // namespace rosim
  
  namespace rosim {
    
    //----------------------------------------------------------------------
    // Constructor
    IonizationReadout::IonizationReadout(fhicl::ParameterSet const& pset)
    {
      // initialize the GArSimulationParameters singleton
      garg4::G4SimulationParameters::CreateInstance(pset.get<fhicl::ParameterSet>("GArSimParsPSet"));

      LOG_DEBUG("IonizationReadout") << "Debug: IonizationReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;
      
      // setup the random number service for Geant4, the "G4Engine" label is a
      // special tag setting up a global engine for use by Geant4/CLHEP;
      // obtain the random seed from NuRandomService,
      // unless overridden in configuration with key "Seed" or "GEANTSeed"
      // same thing for the propagation engine:
      // and again for radio decay
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "propagation", pset, "PropagationSeed");
      
      
      auto driftAlgPars = pset.get<fhicl::ParameterSet>("ElectronDriftAlgPars");
      auto driftAlgName = driftAlgPars.get<std::string>("DriftAlgType");
      
      if(driftAlgName.compare("Standard") == 0)
        fDriftAlg = std::make_unique<gar::rosim::ElectronDriftStandardAlg>(rng->getEngine("propagation"),
                                                                           driftAlgPars);
      else
        throw cet::exception("IonizationReadout")
        << "Unable to determine which electron drift algorithm to use, bail";
      
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
    void IonizationReadout::beginJob()
    {
      ::art::ServiceHandle<geo::Geometry> geom;
      auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());
      
      // create the ionization and scintillation calculator;
      // this is a singleton (!) so it does not make sense
      // to create it in LArVoxelReadoutGeometry
      IonizationAndScintillation::CreateInstance(rng->getEngine("propagation"));
      
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
      std::unique_ptr< std::vector<raw::RawDigit>                       > rdCol (new std::vector<raw::RawDigit>                      );
      std::unique_ptr< ::art::Assns<sdp::EnergyDeposits, raw::RawDigit> > erassn(new ::art::Assns<sdp::EnergyDeposits, raw::RawDigit>);
      
      // drift the ionization electrons to the readout
      //fDriftAlg->DriftElectronsToReadout(dep);
      
      evt.put(std::move(rdCol));
      evt.put(std::move(erassn));
      
      return;
    } // IonizationReadout::produce()
    
  } // namespace rosim
  
  namespace rosim {
    
    DEFINE_ART_MODULE(IonizationReadout)
    
  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_IONIZATIONREADOUT
