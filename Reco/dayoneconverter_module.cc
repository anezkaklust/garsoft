////////////////////////////////////////////////////////////////////////
// Class:       dayoneconverter
// Plugin Type: producer (art v3_00_00)
// File:        dayoneconverter_module.cc
//
// Generated at Tue Jun 23 12:28:55 2020 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include <memory>

#include "TMath.h"
#include "TVector3.h"

#include "Geant4/G4ThreeVector.hh"
#include "nutools/MagneticField/MagneticField.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandGauss.h"

// GArSoft Includes

#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"

class dayoneconverter : public art::EDProducer {
public:
  explicit dayoneconverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  dayoneconverter(dayoneconverter const&) = delete;
  dayoneconverter(dayoneconverter&&) = delete;
  dayoneconverter& operator=(dayoneconverter const&) = delete;
  dayoneconverter& operator=(dayoneconverter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  std::string fInputEdepLabel;  ///<   Input label for edeps
  std::string fInputEdepInstanceTPC;  ///<   Input instance for TPC edeps
  std::string fInputEdepInstanceMuID;  ///<  Input instance for MuID edeps
  float       fSmearX;         ///< amount by which to smear X, in cm
  float       fSmearY;         ///< amount by which to smear Y, in cm
  float       fSmearZ;         ///< amount by which to smear Z, in cm
  cet::exempt_ptr<CLHEP::HepRandomEngine> fEngine;  //< random engine

};


dayoneconverter::dayoneconverter(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  fInputEdepLabel = p.get<std::string>("InputLabel","edepconvert");
  fInputEdepInstanceTPC = p.get<std::string>("InputInstanceTPC","TrackerSc");
  fInputEdepInstanceMuID = p.get<std::string>("InputInstanceMuID","MuID");
  fSmearX = p.get<float>("SmearX",0.3); // in cm
  fSmearY = p.get<float>("SmearY",0.3); // in cm
  fSmearZ = p.get<float>("SmearZ",0.3); // in cm

  art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
  art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
  consumes<std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
  consumes<std::vector<gar::sdp::CaloDeposit> >(muidedeptag);
  produces<std::vector<gar::rec::Track> >();
  produces<std::vector<gar::rec::TPCCluster> >();
  produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();

      // setup the random number service
      // obtain the random seed from NuRandomService
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "dayoneconverter", p, "RandomSeed");

      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;
      auto& engine = rng->getEngine(art::ScheduleID::first(),
				    p.get<std::string>("module_label"),"dayoneconverter");
      fEngine = cet::make_exempt_ptr(&engine);

}

void dayoneconverter::produce(art::Event& e)
{
  // Implementation of required member function here.

      std::unique_ptr< std::vector<gar::rec::TPCCluster> > TPCClusterCol(new std::vector<gar::rec::TPCCluster>);
      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
      //art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
      auto tccdHandle = e.getValidHandle< std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
      auto const& tccds = *tccdHandle;

// put these back when we have tracks and can associate them
//    auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
//    auto const tpcclusPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e);
  
      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      CLHEP::RandGauss GaussRand(*fEngine);

      // first stab, just make all the TPC clusters in the event, and then worry later
      // about pattern recognition

      for (auto const& cd : tccds)
	{
	  double vrand[3];
          GaussRand.fireArray(3, vrand, 0., 1.0);

	  float fcpos[3];
	  fcpos[0] = cd.X() + vrand[0]*fSmearX;
	  fcpos[1] = cd.Y() + vrand[1]*fSmearY;
	  fcpos[2] = cd.Z() + vrand[2]*fSmearZ;
	  
          TPCClusterCol->emplace_back(cd.Energy(),   // do we need a scale factor here?
                                      fcpos,
                                      cd.Time(),     // to check -- units
                                      cd.Time(),
                                      cd.Time(),
                                      fSmearX);	  
	}


      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterCol));
      e.put(std::move(TPCClusterTrkAssns));

}

DEFINE_ART_MODULE(dayoneconverter)
