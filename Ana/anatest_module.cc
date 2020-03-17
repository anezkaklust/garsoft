////////////////////////////////////////////////////////////////////////////////
// Class:       anatest
// Plugin Type: analyzer (art v2_11_02)
// File:        anatest_module.cc
//
// Generated at Tue Mar 3 11:11:11 2020 by Leo Bellantoni
//
// a plugin for testing reco output - derived from anatree_module.cc
//
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "RawDataProducts/CaloRawDigit.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <string>
#include <vector>
#include <unordered_map>

#include <chrono>
using namespace std::chrono;



namespace gar {

    class anatest : public art::EDAnalyzer {
    public:
        explicit anatest(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        anatest(anatest const &) = delete;
        anatest(anatest &&) = delete;
        anatest & operator = (anatest const &) = delete;
        anatest & operator = (anatest &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;



    private:

        void ClearVectors();
        void FillVectors(art::Event const & e);

        // Position of TPC from geometry service; 1 S Boston Ave.
        double ItsInTulsa[3];

        // Input data labels
        std::vector<std::string> fGeneratorLabels;
        std::string fGeantLabel;
        std::string fHitLabel;
        std::string fTPCClusterLabel;
        std::string fTrackLabel;
        std::string fVertexLabel;
        std::string fRawCaloHitLabel;
        std::string fCaloHitLabel;
        std::string fClusterLabel;
        std::string fECALAssnLabel;

        // the analysis output tree
        TTree *fTree;

        // global event info
        Int_t fEvent;        ///< number of the event being processed
        Int_t fRun;          ///< number of the run being processed
        Int_t fSubRun;       ///< number of the sub-run being processed

        // Some random fields for now
        std::vector<Int_t>              fNeutrinoType;
        std::vector<Int_t>              fInteractionType;
        std::vector<Float_t>            fQ2;
    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatest::anatest(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    bool usegenlabels =
        p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
    if (!usegenlabels) fGeneratorLabels.clear();

    fGeantLabel       = p.get<std::string>("GEANTLabel","geant");
    fHitLabel         = p.get<std::string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<std::string>("TrackLabel","track");
    fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
    fRawCaloHitLabel  = p.get<std::string>("RawCaloHitLabel","daqecal");
    fCaloHitLabel     = p.get<std::string>("CaloHitLabel","calohit");
    fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");
    fECALAssnLabel    = p.get<std::string>("ECALAssnLabel","trkecalassn");



    if (usegenlabels) {
        for (size_t i=0; i<fGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::MCTruth> >();
    }



    consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);

    consumes<std::vector<rec::Hit> >(fHitLabel);
    consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);    
    consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Vertex> >(fVertexLabel);
    consumes<art::Assns<rec::Track, rec::Vertex> >(fVertexLabel);

    consumes<std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
    consumes<std::vector<raw::CaloRawDigit> >(fRawCaloHitLabel);
    consumes<std::vector<rec::CaloHit> >(fCaloHitLabel);
    consumes<std::vector<rec::Cluster> >(fClusterLabel);
    consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);
    
    return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::beginJob() {

    art::ServiceHandle<geo::Geometry> euclid;
    ItsInTulsa[0] = euclid->TPCXCent();
    ItsInTulsa[1] = euclid->TPCYCent();
    ItsInTulsa[2] = euclid->TPCZCent();



    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("MPDtestTree","MPDtestTree");



    fTree->Branch("Run",           &fRun,         "Run/I");
    fTree->Branch("SubRun",        &fSubRun,      "SubRun/I");
    fTree->Branch("Event",         &fEvent,       "Event/I");

    fTree->Branch("NType",         &fNeutrinoType);
    fTree->Branch("InterT",        &fInteractionType);
    fTree->Branch("MC_Q2",         &fQ2);

    return;
}  // End of :anatest::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::analyze(art::Event const & e) {

    ClearVectors();
    FillVectors(e);
    fTree->Fill();
    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::ClearVectors() {

    // clear out all our vectors
    fNeutrinoType.clear();
    fInteractionType.clear();
    fQ2.clear();

    return;
} // end :anatest::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::FillVectors(art::Event const & e) {



    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory
    std::vector< art::Handle<std::vector<simb::MCTruth>> > mcthandlelist;
    if (fGeneratorLabels.size()<1) {
        e.getManyByType(mcthandlelist);    // get them all (even if there are none)
    } else {
        mcthandlelist.resize(fGeneratorLabels.size());
        for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
            if (!e.getByLabel(fGeneratorLabels.at(i),mcthandlelist.at(i)))
                throw cet::exception("anatest") << " No simb::MCTruth branch."
                     << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }



    art::Handle< std::vector<simb::MCParticle> > MCPHandle;
    if (!e.getByLabel(fGeantLabel, MCPHandle)) {
        throw cet::exception("anatest") << " No simb::MCParticle branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handles for MCCaloInfo
    art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;
    if (!e.getByLabel(fGeantLabel, SimHitHandle)) {
        throw cet::exception("anatest") << " No gar::sdp::CaloDeposit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for TPC hit data
    art::Handle< std::vector<rec::Hit> > HitHandle;
    if (!e.getByLabel(fHitLabel, HitHandle)) {
        throw cet::exception("anatest") << " No rec::Hit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for TPCClusters
    art::Handle< std::vector<rec::TPCCluster> > TPCClusterHandle;
    if (!e.getByLabel(fTPCClusterLabel, TPCClusterHandle)) {
        throw cet::exception("anatest") << " No rec::TPCCluster branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handles for Tracks and also Assn's to TPCClusters, TrackIoniz
    art::Handle< std::vector<rec::Track> > TrackHandle;
    //art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
    if (!e.getByLabel(fTrackLabel, TrackHandle)) {
        throw cet::exception("anatest") << " No rec::Track branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);



    // Get handle for Vertices; also Assn's to Tracks
    art::Handle< std::vector<rec::Vertex> > VertexHandle;
    //art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
    if (!e.getByLabel(fVertexLabel, VertexHandle)) {
        throw cet::exception("anatest") << " No rec::Vertex branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);



    // Get handle for CaloDigits
    art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
    if (!e.getByLabel(fRawCaloHitLabel, RawHitHandle)) {
        throw cet::exception("anatest") << " No :raw::CaloRawDigit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for CaloHits
    art::Handle< std::vector<rec::CaloHit> > RecoHitHandle;
    if (!e.getByLabel(fCaloHitLabel, RecoHitHandle)) {
        throw cet::exception("anatest") << " No rec::CaloHit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for CaloClusters; also Assn for matching tracks
    art::Handle< std::vector<rec::Cluster> > RecoClusterHandle;
    //art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
    if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {
        throw cet::exception("anatest") << " No rec::Cluster branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>
    //                                   (RecoClusterHandle,e,fECALAssnLabel);



    // =============  Pull art handles =========================================
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();

    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
        for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
            if (mct.NeutrinoSet()) {
                simb::MCNeutrino nuw = mct.GetNeutrino();
                fNeutrinoType.push_back(nuw.Nu().PdgCode());
                fInteractionType.push_back(nuw.InteractionType());
            }
        }
    }





    // Testing the backtracker
	// Here, a few lines to set up useful data and pick some good events.
    auto bt = gar::providerFrom<cheat::BackTracker>();
    typedef int TrkId;

    std::unordered_map<TrkId, int> localMap;    // From track ID to index in event
    int index = 0;
    for ( auto const& mcp : *MCPHandle ) {
        int TrackId = mcp.TrackId();
        localMap[TrackId] = index++;
    }

    auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, HitHandle.id());
    std::vector<art::Ptr<gar::rec::Hit>> allhits;
    size_t nHits = HitHandle->size();
    for (size_t iHit=0; iHit<nHits; ++iHit) {
        auto const hitPointer = hitPtrMaker(iHit);
        allhits.push_back(hitPointer);
    }

    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
		for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
			if (mct.NeutrinoSet()) {
				simb::MCNeutrino nuw = mct.GetNeutrino();
				double x = nuw.Nu().EndX();
				double r = std::hypot( nuw.Nu().EndZ(), nuw.Nu().EndY());
				bool inFiducial = ( r < 445.0 /2.0) && ( fabs(x) < 430.0 /2.0 );
				if (inFiducial) goto somebodyIn;
			}
		}
	}
	return;
	somebodyIn:


    for (simb::MCParticle mcpFromEvent : *MCPHandle) {
        TrkId traky = mcpFromEvent.TrackId();
        simb::MCParticle* mcpFromList = bt->TrackIDToParticle(traky);

		if ( mcpFromList->Mother() != 0) continue;
		double x = mcpFromList->Vx();
		double r = std::hypot( mcpFromList->Vy(), mcpFromList->Vz());
		bool inFiducial = ( r < 445.0 /2.0) && ( fabs(x) < 430.0 /2.0 );
		if (!inFiducial) continue;

		std::vector<art::Ptr<gar::rec::Hit>> somehits =
			bt->ParticleToHits(mcpFromList, allhits);

		std::pair<double,double> straight = bt->HitCollectionEfficiency(mcpFromList,
												somehits,allhits,false);
		std::pair<double,double> chaser   = bt->HitCollectionEfficiency(mcpFromList,
												somehits,allhits,true);
		std::cout << straight.first << " +/- " << straight.second << "\nno,\n" 
			<< chaser.first << " +/- " << chaser.second << std::endl << std::endl;
	}




    return;
}



DEFINE_ART_MODULE(gar::anatest)
