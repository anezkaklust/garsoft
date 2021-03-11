////////////////////////////////////////////////////////////////////////////////
// Class:       anatree
// Plugin Type: analyzer (art v2_11_02)
// File:        anatree_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Additions from Leo Bellantoni, 2019-20
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
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
#include "SimulationDataProducts/GenieParticle.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/Vee.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "RawDataProducts/CaloRawDigit.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH2.h"
#include "TFile.h"
#include "TVector3.h"

#include "CLHEP/Random/RandFlat.h"

#include <string>
#include <vector>
#include <unordered_map>

#include "StandardRecord/StandardRecord.h"

namespace gar {

    typedef std::pair<float, std::string> P;

    class anatree : public art::EDAnalyzer {
    public:
        explicit anatree(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        anatree(anatree const &) = delete;
        anatree(anatree &&) = delete;
        anatree & operator = (anatree const &) = delete;
        anatree & operator = (anatree &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;

    private:
        //Working Horse
        cheat::BackTrackerCore* BackTrack;
        void FillGeneratorMonteCarloInfo(art::Event const & e, caf::StandardRecord& rec);
        void FillRawInfo(art::Event const & e, caf::StandardRecord& rec);
        void FillRecoInfo(art::Event const & e, caf::StandardRecord& rec);
        void FillHighLevelRecoInfo(art::Event const & e, caf::StandardRecord& rec);

        //Helpers
        void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
        float& forwardIonVal, float& backwardIonVal);
        float processOneDirection(std::vector<std::pair<float,float>> SigData,
        float ionizeTruncate);
        // Helper method for processOneDirection
        static bool lessThan_byE(std::pair<float,float> a, std::pair<float,float> b)
        {return a.first < b.first;}

        // Compute T for coherent pion analysis
        float computeT( simb::MCTruth theMCTruth );

        // Calculate the track PID based on reco momentum and Tom's parametrization
        std::vector< std::pair<int, float> > processPIDInfo( float p );

        // Position of TPC from geometry service; 1 S Boston Ave.
        float ItsInTulsa[3];
        float xTPC;
        float rTPC;

        // Input data labels
        std::vector<std::string> fGeneratorLabels;
        std::vector<std::string> fGENIEGeneratorLabels;
        std::string fGeantLabel; ///< module label for geant4 simulated hits
        std::string fGeantInstanceCalo; ///< Instance name ECAL for sdp::CaloDeposit
        std::string fGeantInstanceMuID; ///< Instance name MuID for sdp::CaloDeposit

        std::string fHitLabel; ///< module label for reco TPC hits rec::Hit
        std::string fTPCClusterLabel; ///< module label for TPC Clusters rec::TPCCluster
        std::string fTrackLabel; ///< module label for TPC Tracks rec:Track
        std::string fTrackTrajectoryLabel; ///< module label for TPC Track Trajectories rec:TrackTrajectory
        std::string fVertexLabel; ///< module label for vertexes rec:Vertex
        std::string fVeeLabel; ///< module label for conversion/decay vertexes rec:Vee

        std::string fRawCaloHitLabel; ///< module label for digitized calo hits raw::CaloRawDigit
        std::string fRawCaloHitInstanceCalo; ///< Instance name ECAL for raw::CaloRawDigit
        std::string fRawCaloHitInstanceMuID; ///< Instance name MuID for raw::CaloRawDigit

        std::string fCaloHitLabel; ///< module label for reco calo hits rec::CaloHit
        std::string fCaloHitInstanceCalo; ///< Instance name ECAL for rec::CaloHit
        std::string fMuIDHitLabel;
        std::string fCaloHitInstanceMuID; ///< Instance name MuID for rec::CaloHit

        std::string fClusterLabel; ///< module label for calo clusters rec::Cluster
        std::string fPFLabel; ///< module label for reco particles rec::PFParticle
        std::string fECALAssnLabel; ///< module label for track-clusters associations

        // Optionally keep/drop parts of the analysis tree
        bool  fWriteMCinfo;        ///< Info from MCTruth, GTruth     Default=true
        bool  fWriteMCPTrajectory; ///< Write MCP Trajectory                Default=true
        bool  fWriteMCPTrajMomenta; ///< Write Momenta associated with MCP Trajectory  Default=false
        bool  fWriteMCCaloInfo;    ///< Write MC info for calorimeter Default=true
        float fMatchMCPtoVertDist; ///< MCParticle to MC vertex match Default=roundoff

        bool  fWriteHits;          ///< Write info about TPC Hits     Default=false
        bool  fWriteTPCClusters;   ///< Write TPCClusters info        Default=true
        bool  fWriteTracks;        ///< Start/end X, P for tracks     Default=true
        bool  fWriteTrackTrajectories;        ///< Point traj of reco tracks     Default=false
        bool  fWriteVertices;      ///< Reco vertexes & their tracks  Default=true
        bool  fWriteVees;          ///< Reco vees & their tracks      Default=true

        bool  fWriteCaloDigits;    ///< Raw digits for calorimetry.   Default=false
        bool  fWriteCaloHits;      ///< Write ECAL hits.              Default=true
        bool  fWriteCaloClusters;  ///< Write ECAL clusters.          Default=true
        bool  fWriteMatchedTracks; ///< Write ECAL-track Assns        Default=true

        // Truncation parameter for dE/dx (average this fraction of the lowest readings)
        float fIonizTruncate;      ///<                               Default=1.00;
        bool  fWriteCohInfo;       ///< MC level t for coherent pi+.  Default=false

        // the analysis tree
        TTree *fTree;

        Float_t fTPC_X;      ///< center of TPC stored as per-event & compressed by root
        Float_t fTPC_Y;
        Float_t fTPC_Z;


        //Geometry
        const geo::GeometryCore* fGeo; ///< pointer to the geometry

        typedef int TrkId;
        std::unordered_map<TrkId, Int_t> TrackIdToIndex;

        //map of PID TH2 per momentum value
        CLHEP::HepRandomEngine &fEngine;  ///< random engine
        std::unordered_map<int, TH2F*> m_pidinterp;
    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatree::anatree(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed")) {
    fGeo     = gar::providerFrom<geo::GeometryGAr>();

    bool usegenlabels =
    p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
    if (!usegenlabels) fGeneratorLabels.clear();

    bool usegeniegenlabels  =
    p.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
    if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

    //Sim Hits
    fGeantLabel        = p.get<std::string>("GEANTLabel","geant");
    fGeantInstanceCalo = p.get<std::string>("GEANTInstanceCalo","ECAL");
    fGeantInstanceMuID = p.get<std::string>("GEANTInstanceMuID","MuID");

    fHitLabel         = p.get<std::string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<std::string>("TrackLabel","track");
    fTrackTrajectoryLabel  = p.get<std::string>("TrackTrajectoryLabel","track");
    fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
    fVeeLabel         = p.get<std::string>("VeeLabel","veefinder1");

    //Calorimetric related ECAL/MuID
    fRawCaloHitLabel        = p.get<std::string>("RawCaloHitLabel","daqecal");
    fRawCaloHitInstanceCalo = p.get<std::string>("RawCaloHitInstanceCalo","ECAL");
    fRawCaloHitInstanceMuID = p.get<std::string>("RawCaloHitInstanceMuID","MuID");

    fCaloHitLabel        = p.get<std::string>("CaloHitLabel","sipmhit");
    fCaloHitInstanceCalo = p.get<std::string>("CaloHitInstanceCalo","ECAL");
    fMuIDHitLabel        = p.get<std::string>("MuIDHitLabel","sipmhit");
    fCaloHitInstanceMuID = p.get<std::string>("CaloHitInstanceMuID","MuID");

    fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");
    fPFLabel          = p.get<std::string>("PFLabel","pandora");
    fECALAssnLabel    = p.get<std::string>("ECALAssnLabel","trkecalassn");

    // What to write
    fWriteMCinfo              = p.get<bool>("WriteMCinfo",        true);
    fWriteMCPTrajectory       = p.get<bool>("WriteMCPTrajectory", true);
    fWriteMCPTrajMomenta      = p.get<bool>("WriteMCPTrajMomenta",false);
    fWriteMCCaloInfo          = p.get<bool>("WriteMCCaloInfo",    true);
    float MCPtoVertDefault    = 10.0*std::numeric_limits<Float_t>::epsilon();
    fMatchMCPtoVertDist       = p.get<float>("MatchMCPtoVertDist",MCPtoVertDefault);

    fWriteHits                = p.get<bool>("WriteHits",         false);
    fWriteTPCClusters         = p.get<bool>("WriteTPCClusters",  true);
    fWriteTracks              = p.get<bool>("WriteTracks",       true);
    fWriteTrackTrajectories   = p.get<bool>("WriteTrackTrajectories", false);
    fWriteVertices            = p.get<bool>("WriteVertices",     true);
    fWriteVees                = p.get<bool>("WriteVees",         true);

    fWriteCaloDigits          = p.get<bool>("WriteCaloDigits",   false);
    fWriteCaloHits            = p.get<bool>("WriteCaloHits",     true);
    fWriteCaloClusters        = p.get<bool>("WriteCaloClusters", true);
    fWriteMatchedTracks       = p.get<bool>("WriteMatchedTracks",true);

    fIonizTruncate            = p.get<float>("IonizTruncate",    0.70);
    fWriteCohInfo             = p.get<bool> ("WriteCohInfo",     false);

    if (usegenlabels) {
        for (size_t i=0; i<fGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::MCTruth> >();
    }

    if (usegeniegenlabels) {
        for (size_t i=0; i<fGENIEGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::GTruth> >(fGENIEGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::GTruth> >();
    }

    consumesMany<std::vector<sdp::GenieParticle> >();
    //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
    consumes<std::vector<simb::MCParticle> >(fGeantLabel);

    //TPC related
    consumes<std::vector<sdp::EnergyDeposit> >(fGeantLabel);
    consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Hit> >(fHitLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);
    consumes<std::vector<rec::TrackTrajectory> >(fTrackTrajectoryLabel);
    consumes<std::vector<rec::Vertex> >(fVertexLabel);
    consumes<art::Assns<rec::Track, rec::Vertex> >(fVertexLabel);
    consumes<std::vector<rec::Vee> >(fVeeLabel);
    consumes<art::Assns<rec::Track, rec::Vee> >(fVeeLabel);
    consumes<std::vector<rec::TrackIoniz>>(fTrackLabel);
    consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

    //Calorimetry related
    art::InputTag ecalgeanttag(fGeantLabel, fGeantInstanceCalo);
    consumes<std::vector<gar::sdp::CaloDeposit> >(ecalgeanttag);

    art::InputTag ecalrawtag(fRawCaloHitLabel, fRawCaloHitInstanceCalo);
    consumes<std::vector<raw::CaloRawDigit> >(ecalrawtag);

    art::InputTag ecalhittag(fCaloHitLabel, fCaloHitInstanceCalo);
    consumes<std::vector<rec::CaloHit> >(ecalhittag);

    //Muon system related
    if (fGeo->HasMuonDetector()) {
        art::InputTag muidgeanttag(fGeantLabel, fGeantInstanceMuID);
        consumes<std::vector<gar::sdp::CaloDeposit> >(muidgeanttag);

        art::InputTag muidrawtag(fRawCaloHitLabel, fRawCaloHitInstanceMuID);
        consumes<std::vector<raw::CaloRawDigit> >(muidrawtag);

        art::InputTag muidhittag(fCaloHitLabel, fCaloHitInstanceMuID);
        consumes<std::vector<rec::CaloHit> >(muidhittag);
    }

    consumes<std::vector<rec::Cluster> >(fClusterLabel);
    consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

    return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::beginJob() {

    fTPC_X = ItsInTulsa[0] = fGeo->TPCXCent();
    fTPC_Y = ItsInTulsa[1] = fGeo->TPCYCent();
    fTPC_Z = ItsInTulsa[2] = fGeo->TPCZCent();

    xTPC = fGeo->TPCLength() / 2.;
    rTPC = fGeo->TPCRadius();

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

    // Create the branch. We will update the address before we write the tree
    caf::StandardRecord* rec = 0;
    fTree->Branch("rec", &rec);

    if (fWriteMCPTrajectory) {
        // Checking for parameter file validity only once to simplify
        // later code in anatree::process etc.
        if (!fWriteMCinfo) {
            throw cet::exception("anatree")
            << " fWriteMCPTrajectory, but !fWriteMCinfo."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        } 
    }

    if (fWriteMCCaloInfo) {
        if (!fWriteMCinfo) {
            throw cet::exception("anatree")
            << " fWriteMCCaloInfo, but !fWriteMCinfo."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Reco'd verts & their track-ends
    if (fWriteVertices) {
        if (!fWriteTracks) {
            throw cet::exception("anatree")
            << " fWriteVertices, but !fWriteTracks."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Reco'd vees & their track-ends
    if (fWriteVees) {
        if (!fWriteTracks) {
            throw cet::exception("anatree")
            << " fWriteVees, but !fWriteTracks."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    if (fWriteMatchedTracks) {
        if (!fWriteTracks || !fWriteCaloClusters) {
            throw cet::exception("anatree")
            << " fWriteMatchedTracks, but (!fWriteTracks || !fWriteCaloClusters)."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    std::string filename = "${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
    TFile infile(filename.c_str(), "READ");

    m_pidinterp.clear();
    char str[11];
    for (int q = 0; q < 501; ++q) {
        sprintf(str, "%d", q);
        std::string s = "pidmatrix";
        s.append(str);
        // read the 500 histograms one by one; each histogram is a
        // 6 by 6 matrix of probabilities for a given momentum value
        m_pidinterp.insert( std::make_pair(q, (TH2F*) infile.Get(s.c_str())->Clone("pidinterp")) );
    }

    return;
}  // End of :anatree::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::analyze(art::Event const & e) {

    caf::StandardRecord rec;
    caf::StandardRecord* prec = &rec;
    fTree->SetBranchAddress("rec", &prec);

    rec.run    = e.run();
    rec.subrun = e.subRun();
    rec.event  = e.id().event();

    rec.TPC_X = fTPC_X;
    rec.TPC_Y = fTPC_Y;
    rec.TPC_Z = fTPC_Z;


    // Need a non-constant backtracker instance, for now, in anayze ot beginJob
    cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
    BackTrack = const_cast<cheat::BackTrackerCore*>(const_bt);

    //Fill generator and MC Information
    if (fWriteMCinfo) FillGeneratorMonteCarloInfo(e, rec);

    //Fill Raw Information
    if (fWriteCaloDigits) FillRawInfo(e, rec);

    //Fill Reco Information
    FillRecoInfo(e, rec);

    //Fill HL Reco Information
    FillHighLevelRecoInfo(e, rec);

    fTree->Fill();
    return;
}


//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillGeneratorMonteCarloInfo(art::Event const & e, caf::StandardRecord& rec) {

    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
    std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;

    if (fGeneratorLabels.size()<1) {
        e.getManyByType(mcthandlelist);    // get them all (even if there are none)
    } else {
        mcthandlelist.resize(fGeneratorLabels.size());
        for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
            // complain if we wanted a specific one but didn't find it
            if (!e.getByLabel(fGeneratorLabels.at(i),mcthandlelist.at(i))) {
                throw cet::exception("anatree") << " No simb::MCTruth branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }
    }

    if (fGENIEGeneratorLabels.size()<1) {
        e.getManyByType(gthandlelist);  // get them all (even if there are none)
    } else {
        gthandlelist.resize(fGENIEGeneratorLabels.size());
        for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i) {
            // complain if we wanted a specific one but didn't find it
            if (!e.getByLabel(fGENIEGeneratorLabels.at(i),gthandlelist.at(i))) {
                throw cet::exception("anatree") << " No simb::GTruth branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }
    }

    // save MCTruth info
    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
        for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
            if (mct.NeutrinoSet()) {
                simb::MCNeutrino nuw = mct.GetNeutrino();
                rec.neutrinoType.push_back(nuw.Nu().PdgCode());
                rec.ccnc.push_back(nuw.CCNC());
                rec.mode.push_back(nuw.Mode());
                rec.interactionType.push_back(nuw.InteractionType());
                rec.Q2.push_back(nuw.QSqr());
                rec.W.push_back(nuw.W());
                rec.X.push_back(nuw.X());
                rec.Y.push_back(nuw.Y());
                rec.theta.push_back(nuw.Theta());
                if (fWriteCohInfo) {
                  rec.T.push_back(computeT(mct));
                }
                rec.MCVertexX.push_back(nuw.Nu().EndX());
                rec.MCVertexY.push_back(nuw.Nu().EndY());
                rec.MCVertexZ.push_back(nuw.Nu().EndZ());
                rec.MCnuPx.push_back(nuw.Nu().Px());
                rec.MCnuPy.push_back(nuw.Nu().Py());
                rec.MCnuPz.push_back(nuw.Nu().Pz());
            }  // end MC info from MCTruth
        }
    }

    // save GTruth info
    for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl) {
        for ( auto const& gt : (*gthandlelist.at(igthl)) ) {
            rec.gint.push_back(gt.fGint);
            rec.tgtPDG.push_back(gt.ftgtPDG);
            rec.weight.push_back(gt.fweight);
            rec.gT.push_back(gt.fgT);
        }
    }

    // Save the particle list from the GENIE event record
    std::vector< art::Handle< std::vector<sdp::GenieParticle> > > gparthandlelist;
    e.getManyByType(gparthandlelist);

    for (size_t igphl = 0; igphl < gparthandlelist.size(); ++igphl) {
        unsigned int nGPart = 0;
        for ( auto const& gpart : (*gparthandlelist.at(igphl)) ) {
            rec.GPartIntIdx.push_back(gpart.InteractionIndex());
            rec.GPartIdx.push_back(gpart.Index());
            rec.GPartPdg.push_back(gpart.Pdg());
            rec.GPartStatus.push_back(gpart.Status());
            rec.GPartName.push_back(gpart.Name());
            rec.GPartFirstMom.push_back(gpart.FirstMother());
            rec.GPartLastMom.push_back(gpart.LastMother());
            rec.GPartFirstDaugh.push_back(gpart.FirstDaughter());
            rec.GPartLastDaugh.push_back(gpart.LastDaughter());
            rec.GPartPx.push_back(gpart.Px());
            rec.GPartPy.push_back(gpart.Py());
            rec.GPartPz.push_back(gpart.Pz());
            rec.GPartE.push_back(gpart.E());
            rec.GPartMass.push_back(gpart.Mass());
            nGPart++;
        }
        rec.nGPart.push_back(nGPart);
    }

    art::Handle< std::vector<simb::MCParticle> > MCPHandle;
    if (!e.getByLabel(fGeantLabel, MCPHandle)) {
        throw cet::exception("anatree") << " No simb::MCParticle branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // Save MCParticle info (post-GEANT; for pre-GEANT, maybe try MCTruth?)
    // Helps to have a map from the TrackId from generator to index number
    // into *MCPHandle, which is a vector<MCParticle>.  Not all TrackId value
    // appear in MCPHandle, although I think the missing ones are all from
    // the GEANT rather than the GENIE stage of GENIEGen.  (*MCPHandle).size()
    // is the same as eg fMCPStartX.size() by the time this entry is written.
    // TracIdToIndex isa class variable for accessibility by per-Track and
    // per-Cluster code
    Int_t index = 0;
    for ( auto const& mcp : (*MCPHandle) ) {
        int TrackId = mcp.TrackId();
        TrackIdToIndex[TrackId] = index++;
    }

    for ( auto const& mcp : (*MCPHandle) ) {
        rec.MCTrkID.push_back(mcp.TrackId());
        rec.MCPDG.push_back(mcp.PdgCode());

        // If mcp.Mother() == 0, particle is from initial vertex;
        // Set MCMotherIndex to -1 in that case.
        TrkId momTrkId = mcp.Mother();
        Int_t momIndex = -1;
        int momPDG = 0;
        if (momTrkId>0) {
            //Check if it exists!
            if(TrackIdToIndex.find(momTrkId) != TrackIdToIndex.end()){
                momIndex = TrackIdToIndex[momTrkId];
                momPDG   = (*MCPHandle).at(momIndex).PdgCode();
            } else {
                MF_LOG_DEBUG("Anatree_module")
                << " mcp trkid " << mcp.TrackId()
                << " pdg code " << mcp.PdgCode()
                << " could not find mother trk id " << momTrkId
                << " in the TrackIdToIndex map"
                << " creating process is [ " << mcp.Process() << " ]";
            }
        }

        rec.MCMotherIndex.push_back(momIndex);
        rec.MCMotherTrkID.push_back(mcp.Mother());//directly trackid not index
        rec.MCPDGMother.push_back(momPDG);

        const TLorentzVector& position = mcp.Position(0);
        const TLorentzVector& momentum = mcp.Momentum(0);
        rec.MCPStartX.push_back(position.X());
        rec.MCPStartY.push_back(position.Y());
        rec.MCPStartZ.push_back(position.Z());
        rec.MCPTime.push_back(mcp.T());
        rec.MCPStartPX.push_back(momentum.Px());
        rec.MCPStartPY.push_back(momentum.Py());
        rec.MCPStartPZ.push_back(momentum.Pz());

        const TLorentzVector& positionEnd = mcp.EndPosition();
        const TLorentzVector& momentumEnd = mcp.EndMomentum();
        rec.MCPEndX.push_back(positionEnd.X());
        rec.MCPEndY.push_back(positionEnd.Y());
        rec.MCPEndZ.push_back(positionEnd.Z());
        rec.MCPEndPX.push_back(momentumEnd.Px());
        rec.MCPEndPY.push_back(momentumEnd.Py());
        rec.MCPEndPZ.push_back(momentumEnd.Pz());
        rec.MCPProc.push_back(mcp.Process());
        rec.MCPEndProc.push_back(mcp.EndProcess());
    }

    // Get vertex for each MCParticle.  In principle, the mct.GetParticle()
    // method should produce all of them; filter out the ones with a non-stable
    // status code and the GENIE pseudo-particles (Pid>= 2000000000) and you
    // should have the contents of (*MCPHandle).  But you don't.  On the 1st
    // event I (Leo Bellantoni) tested, with 40 vertices and 351 primary
    // particles, particle 143, a proton, evidently decayed or something and is
    // not in (*MCPHandle).  So we use the following primitive earth technology.
    size_t nMCParticles = (*MCPHandle).size();
    size_t iMCParticle  = 0;
    rec.MCPVertIndex.resize(nMCParticles);
    for (; iMCParticle<nMCParticles; ++iMCParticle) {
        foundMCvert:
            // Assign noprimary to start with
            rec.MCPVertIndex[iMCParticle] = -1;
            // Do the primaries first
            if (rec.MCMotherIndex[iMCParticle]!=-1) break;
            Float_t trackX = rec.MCPStartX[iMCParticle];
            Float_t trackY = rec.MCPStartY[iMCParticle];
            Float_t trackZ = rec.MCPStartZ[iMCParticle];
            int vertexIndex = 0;
            for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
                for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
                    if (mct.NeutrinoSet()) {
                        simb::MCNeutrino nuw = mct.GetNeutrino();
                        Float_t vertX = nuw.Nu().EndX();
                        Float_t vertY = nuw.Nu().EndY();
                        Float_t vertZ = nuw.Nu().EndZ();
                        Float_t dist = std::hypot(trackX-vertX,trackY-vertY,trackZ-vertZ);
                        if ( dist <= fMatchMCPtoVertDist ) {
                            rec.MCPVertIndex[iMCParticle] = vertexIndex;
                            ++iMCParticle; goto foundMCvert;
                        }
                    }
                    ++vertexIndex;
                }
            }
    }

    // Now the secondaries.  As they are after the primaries, do not re-init iMCParticle
    for (; iMCParticle<nMCParticles; ++iMCParticle) {
        int momIndex = rec.MCMotherIndex[iMCParticle];
        int lastMCParticle = iMCParticle;
        while (momIndex != -1) {
            lastMCParticle = momIndex;
            momIndex       = rec.MCMotherIndex[momIndex];
        }
        rec.MCPVertIndex[iMCParticle] = rec.MCPVertIndex[lastMCParticle];
    }

    if (fWriteMCPTrajectory) {
        // It's in the MCParticle table
        Int_t mcpIndex = 0;
        for ( auto const& mcp : (*MCPHandle) ) {
            const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
            const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
            //No charge don't store the trajectory
            if (definition==nullptr || definition->Charge() == 0) continue;
            //TrackID of the mcp to keep track to which mcp this trajectory is
            int trackId = mcp.TrackId();
            for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++) {
                float xTraj = mcp.Trajectory().X(iTraj);
                float yTraj = mcp.Trajectory().Y(iTraj);
                float zTraj = mcp.Trajectory().Z(iTraj);
                float rTraj = std::hypot( yTraj - ItsInTulsa[1], zTraj - ItsInTulsa[2]);

                if (abs(xTraj - ItsInTulsa[0]) > xTPC) continue;
                if (rTraj > rTPC) continue;

                rec.TrajMCPX.push_back(xTraj);
                rec.TrajMCPY.push_back(yTraj);
                rec.TrajMCPZ.push_back(zTraj);
                rec.TrajMCPT.push_back(mcp.Trajectory().T(iTraj));
                if (fWriteMCPTrajMomenta) {
                    rec.TrajMCPPX.push_back(mcp.Trajectory().Px(iTraj));
                    rec.TrajMCPPY.push_back(mcp.Trajectory().Py(iTraj));
                    rec.TrajMCPPZ.push_back(mcp.Trajectory().Pz(iTraj));
                }
                rec.TrajMCPE.push_back(mcp.Trajectory().E(iTraj));
                rec.TrajMCPIndex.push_back(mcpIndex);
                rec.TrajMCPTrackID.push_back(trackId);
            }
            mcpIndex++;
        }
    }

    // Get handles for MCCaloInfo
    art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;//ecal
    art::Handle< std::vector<gar::sdp::CaloDeposit> > MuIDSimHitHandle;//ecal
    art::InputTag ecalgeanttag(fGeantLabel, fGeantInstanceCalo);
    art::InputTag muidgeanttag(fGeantLabel, fGeantInstanceMuID);
    if (fWriteMCCaloInfo) {
        if (!e.getByLabel(ecalgeanttag, SimHitHandle)) {
            throw cet::exception("anatree") << " No gar::sdp::CaloDeposit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fGeo->HasMuonDetector() && !e.getByLabel(muidgeanttag, MuIDSimHitHandle)) {
            throw cet::exception("anatree") << " No gar::sdp::CaloDeposit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // Save simulation ecal hit info
        for ( auto const& SimHit : (*SimHitHandle) ) {
            rec.SimnHits++;
            rec.SimHitX.push_back(SimHit.X());
            rec.SimHitY.push_back(SimHit.Y());
            rec.SimHitZ.push_back(SimHit.Z());
            rec.SimHitTime.push_back(SimHit.Time());
            rec.SimHitEnergy.push_back(SimHit.Energy());
            rec.SimHitTrackID.push_back(SimHit.TrackID());
            rec.SimHitCellID.push_back(SimHit.CellID());
            rec.SimEnergySum += SimHit.Energy();
        }

        // Save simulation muon system hit info
        if(fGeo->HasMuonDetector()) {
            for ( auto const& SimHit : (*MuIDSimHitHandle) ) {
                rec.SimnHits_MuID++;
                rec.SimHitX_MuID.push_back(SimHit.X());
                rec.SimHitY_MuID.push_back(SimHit.Y());
                rec.SimHitZ_MuID.push_back(SimHit.Z());
                rec.SimHitTime_MuID.push_back(SimHit.Time());
                rec.SimHitEnergy_MuID.push_back(SimHit.Energy());
                rec.SimHitTrackID_MuID.push_back(SimHit.TrackID());
                rec.SimHitCellID_MuID.push_back(SimHit.CellID());
                rec.SimEnergySum_MuID += SimHit.Energy();
            }
        }
    }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillRawInfo(art::Event const & e, caf::StandardRecord& rec) {

    // Get handle for CaloDigits
    art::InputTag ecalrawtag(fRawCaloHitLabel, fRawCaloHitInstanceCalo);
    art::InputTag muidrawtag(fRawCaloHitLabel, fRawCaloHitInstanceMuID);
    art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
    art::Handle< std::vector<gar::raw::CaloRawDigit> > MuIDRawHitHandle;

    if (!e.getByLabel(ecalrawtag, RawHitHandle)) {
        throw cet::exception("anatree") << " No :raw::CaloRawDigit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (fGeo->HasMuonDetector() && !e.getByLabel(muidrawtag, MuIDRawHitHandle)) {
        throw cet::exception("anatree") << " No :raw::CaloRawDigit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // save ecal raw digits info
    for ( auto const& DigiHit : (*RawHitHandle) ) {
        rec.DiginHits++;
        rec.DigiHitX.push_back(DigiHit.X());
        rec.DigiHitY.push_back(DigiHit.Y());
        rec.DigiHitZ.push_back(DigiHit.Z());
        rec.DigiHitTime.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
        rec.DigiHitADC.push_back(DigiHit.ADC().first);
        rec.DigiHitCellID.push_back(DigiHit.CellID());
    }

    // save muon system raw digits info
    if(fGeo->HasMuonDetector()) {
        for ( auto const& DigiHit : (*MuIDRawHitHandle) ) {
            rec.DiginHits_MuID++;
            rec.DigiHitX_MuID.push_back(DigiHit.X());
            rec.DigiHitY_MuID.push_back(DigiHit.Y());
            rec.DigiHitZ_MuID.push_back(DigiHit.Z());
            rec.DigiHitTime_MuID.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
            rec.DigiHitADC_MuID.push_back(DigiHit.ADC().first);
            rec.DigiHitCellID_MuID.push_back(DigiHit.CellID());
        }
    }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillRecoInfo(art::Event const & e, caf::StandardRecord& rec) {

    // Get handle for TPC hit data
    art::Handle< std::vector<rec::Hit> > HitHandle;
    if (fWriteHits) {
        if (!e.getByLabel(fHitLabel, HitHandle)) {
            throw cet::exception("anatree") << " No rec::Hit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // save hits in the TPC
        for ( auto const& Hit : (*HitHandle) ) {
            rec.HitX.push_back(Hit.Position()[0]);
            rec.HitY.push_back(Hit.Position()[1]);
            rec.HitZ.push_back(Hit.Position()[2]);
            rec.HitSig.push_back(Hit.Signal());
            rec.HitRMS.push_back(Hit.RMS());
            rec.HitChan.push_back(Hit.Channel());
        }
    }

    // Get handle for CaloHits
    art::InputTag ecalrecotag(fCaloHitLabel, fCaloHitInstanceCalo);
    art::InputTag muirecotag(fMuIDHitLabel, fCaloHitInstanceMuID);
    art::Handle< std::vector<rec::CaloHit> > RecoHitHandle;
    art::Handle< std::vector<rec::CaloHit> > MuIDRecoHitHandle;

    if (fWriteCaloHits) {

        if (!e.getByLabel(ecalrecotag, RecoHitHandle)) {
            throw cet::exception("anatree") << " No rec::CaloHit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fGeo->HasMuonDetector() && !e.getByLabel(muirecotag, MuIDRecoHitHandle)) {
            throw cet::exception("anatree") << " No rec::CaloHit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // save reco'd Calorimetry hits
        for ( auto const& Hit : (*RecoHitHandle) ) {
            rec.ReconHits++;
            rec.ReconHitIDNumber.push_back(Hit.getIDNumber());
            rec.RecoHitX.push_back(Hit.Position()[0]);
            rec.RecoHitY.push_back(Hit.Position()[1]);
            rec.RecoHitZ.push_back(Hit.Position()[2]);
            rec.RecoHitTime.push_back(Hit.Time().first);
            rec.RecoHitEnergy.push_back(Hit.Energy());
            rec.RecoHitCellID.push_back(Hit.CellID());
            rec.RecoEnergySum += Hit.Energy();
        }

        if(fGeo->HasMuonDetector()) {
            for ( auto const& Hit : (*MuIDRecoHitHandle) ) {
                rec.ReconHits_MuID++;
                rec.ReconHitIDNumber_MuID.push_back(Hit.getIDNumber());
                rec.RecoHitX_MuID.push_back(Hit.Position()[0]);
                rec.RecoHitY_MuID.push_back(Hit.Position()[1]);
                rec.RecoHitZ_MuID.push_back(Hit.Position()[2]);
                rec.RecoHitTime_MuID.push_back(Hit.Time().first);
                rec.RecoHitEnergy_MuID.push_back(Hit.Energy());
                rec.RecoHitCellID_MuID.push_back(Hit.CellID());
                rec.RecoEnergySum_MuID += Hit.Energy();
            }
        }
    }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillHighLevelRecoInfo(art::Event const & e, caf::StandardRecord& rec) {

    // Get handle for TPCClusters
    art::Handle< std::vector<rec::TPCCluster> > TPCClusterHandle;
    if (fWriteTPCClusters) {
        if (!e.getByLabel(fTPCClusterLabel, TPCClusterHandle)) {
            throw cet::exception("anatree") << " No rec::TPCCluster branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handles for Tracks and their ionizations; also Assn's to TPCClusters, TrackIoniz
    art::Handle< std::vector<rec::Track> > TrackHandle;
    art::Handle< std::vector<rec::TrackIoniz> > TrackIonHandle;
    art::Handle< std::vector<rec::TrackTrajectory> > TrackTrajHandle;
    art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
    art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
    if(fWriteTracks) {
        if (!e.getByLabel(fTrackLabel, TrackHandle)) {
            throw cet::exception("anatree") << " No rec::Track branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (!e.getByLabel(fTrackLabel, TrackIonHandle)) {
            throw cet::exception("anatree") << " No rec::TrackIoniz branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);
        findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,e,fTrackLabel);

        if(fWriteTrackTrajectories) {
             if (!e.getByLabel(fTrackTrajectoryLabel, TrackTrajHandle)) {
                throw cet::exception("anatree") << " No rec::TrackTrajectory branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }
    }

    // Get handle for Vertices; also Assn's to Tracks
    art::Handle< std::vector<rec::Vertex> > VertexHandle;
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
    if (fWriteVertices) {
        if (!e.getByLabel(fVertexLabel, VertexHandle)) {
            throw cet::exception("anatree") << " No rec::Vertex branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);
    }

    // Get handle for Vees; also Assn's to Tracks
    art::Handle< std::vector<rec::Vee> > VeeHandle;
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyVeeTrackEnd = NULL;
    if (fWriteVees) {
        if (!e.getByLabel(fVeeLabel, VeeHandle)) {
            throw cet::exception("anatree") << " No rec::Vee branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        findManyVeeTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VeeHandle,e,fVeeLabel);
    }

    // Get handle for CaloClusters; also Assn for matching tracks
    art::Handle< std::vector<rec::Cluster> > RecoClusterHandle;
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
    if (fWriteCaloClusters) {
        if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {
            throw cet::exception("anatree") << " No rec::Cluster branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fWriteTracks) {
            findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>
            (RecoClusterHandle,e,fECALAssnLabel);
        }
    }

    // save clusters in the TPC. For some reason, can't get FindOneP<rec::Track> or
    // FindManyP<rec::Track> to work; seems the underlying Assn isn't found.  Have
    // to FindManyP<TPCCluster> instead and  iterate if (fWriteTracks).  :(
    if (fWriteTPCClusters) {
        for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
            rec.TPCClusterX.push_back(TPCCluster.Position()[0]);
            rec.TPCClusterY.push_back(TPCCluster.Position()[1]);
            rec.TPCClusterZ.push_back(TPCCluster.Position()[2]);
            rec.TPCClusterSig.push_back(TPCCluster.Signal());
            rec.TPCClusterRMS.push_back(TPCCluster.RMS());
            const float* cov;
            cov = TPCCluster.CovMatPacked();
            rec.TPCClusterCovXX.push_back(cov[0]);
            rec.TPCClusterCovXY.push_back(cov[1]);
            rec.TPCClusterCovXZ.push_back(cov[2]);
            rec.TPCClusterCovYY.push_back(cov[3]);
            rec.TPCClusterCovYZ.push_back(cov[4]);
            rec.TPCClusterCovZZ.push_back(cov[5]);

            Int_t trackForThisTPCluster = -1;
            if (fWriteTracks) {
                size_t iTrack = 0;
                for ( auto const& track : (*TrackHandle) ) {
                    for (size_t iCluster=0; iCluster<track.NHits(); iCluster++) {
                        auto const& trackedCluster =
                        *(findManyTPCClusters->at(iTrack).at(iCluster));
                        if (TPCCluster==trackedCluster) {
                            trackForThisTPCluster = track.getIDNumber();
                            // No cluster is in 2 tracks (don't mess up dE/dx!)
                            goto pushit;   // break 2 loops
                        }
                    }
                    iTrack++;
                }
            }
            pushit:
                rec.TPCClusterTrkIDNumber.push_back(trackForThisTPCluster);
        }
    }

    // save per-track info
    if (fWriteTracks) {
        size_t iTrack = 0;
        for ( auto const& track : (*TrackHandle) ) {
            // track is a rec::Track, not a rec::TrackPar
            rec.TrackIDNumber.push_back(track.getIDNumber());

            rec.TrackStartX.push_back(track.Vertex()[0]);
            rec.TrackStartY.push_back(track.Vertex()[1]);
            rec.TrackStartZ.push_back(track.Vertex()[2]);
            rec.TrackStartPX.push_back(track.Momentum_beg()*track.VtxDir()[0]);
            rec.TrackStartPY.push_back(track.Momentum_beg()*track.VtxDir()[1]);
            rec.TrackStartPZ.push_back(track.Momentum_beg()*track.VtxDir()[2]);
            rec.TrackStartQ.push_back(track.ChargeBeg());

            rec.TrackEndX.push_back(track.End()[0]);
            rec.TrackEndY.push_back(track.End()[1]);
            rec.TrackEndZ.push_back(track.End()[2]);
            rec.TrackEndPX.push_back(track.Momentum_end()*track.EndDir()[0]);
            rec.TrackEndPY.push_back(track.Momentum_end()*track.EndDir()[1]);
            rec.TrackEndPZ.push_back(track.Momentum_end()*track.EndDir()[2]);
            rec.TrackEndQ.push_back(track.ChargeEnd());

            rec.TrackLenF.push_back(track.LengthForward());
            rec.TrackLenB.push_back(track.LengthBackward());
            rec.TrackChi2F.push_back(track.ChisqForward());
            rec.TrackChi2B.push_back(track.ChisqBackward());
            rec.NTPCClustersOnTrack.push_back(track.NHits());

            //Add the PID information based on Tom's parametrization
            TVector3 momF(track.Momentum_beg()*track.VtxDir()[0], track.Momentum_beg()*track.VtxDir()[1], track.Momentum_beg()*track.VtxDir()[2]);
            TVector3 momB(track.Momentum_end()*track.EndDir()[0], track.Momentum_end()*track.EndDir()[1], track.Momentum_end()*track.EndDir()[2]);

            //Reconstructed momentum forward and backward
            float pF = momF.Mag();
            float pB = momB.Mag();
            std::vector< std::pair<int, float> > pidF = processPIDInfo( pF );
            std::vector< std::pair<int, float> > pidB = processPIDInfo( pB );

            //Fill the pid and its probability
            for(size_t ipid = 0; ipid < pidF.size(); ipid++) {
                rec.TrackPIDF.push_back( pidF.at(ipid).first );
                rec.TrackPIDProbF.push_back( pidF.at(ipid).second );
            }
            for(size_t ipid = 0; ipid < pidB.size(); ipid++) {
                rec.TrackPIDB.push_back( pidB.at(ipid).first );
                rec.TrackPIDProbB.push_back( pidB.at(ipid).second );
            }

            // Matching MCParticle info
            std::vector<std::pair<simb::MCParticle*,float>> trakt;
            trakt = BackTrack->TrackToMCParticles( const_cast<rec::Track*>(&track) );
            int eileen = -1;
            if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
                eileen = TrackIdToIndex[trakt[0].first->TrackId()];
            }
            rec.TrackMCindex.push_back(eileen);
            if (eileen > -1) {
            	rec.TrackMCfrac.push_back(trakt[0].second);
            } else {
            	rec.TrackMCfrac.push_back(0.0);
            }

            if (findIonization->isValid()) {
                // No calibration for now.  Someday this should all be in reco
                rec::TrackIoniz ionization = *(findIonization->at(iTrack));
                float avgIonF, avgIonB;
                processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
                rec.TrackAvgIonF.push_back( avgIonF );
                rec.TrackAvgIonB.push_back( avgIonB );
            } else {
                // must push_back something so that rec.TrackAvgIonF,B are of correct size.
                rec.TrackAvgIonF.push_back( 0.0 );
                rec.TrackAvgIonB.push_back( 0.0 );
            }
            iTrack++;
        } // end loop over TrackHandle
    }

    //TrackTrajectories
    if(fWriteTrackTrajectories) {
        size_t iTrackTraj = 0;
        for ( auto const& tracktraj : (*TrackTrajHandle) ) {
            
            std::vector<TVector3> temp = tracktraj.getFWDTrajectory();
            for(size_t i = 0; i < temp.size(); i++) {
                rec.TrackTrajectoryFWDX.push_back(temp.at(i).X());
                rec.TrackTrajectoryFWDY.push_back(temp.at(i).Y());
                rec.TrackTrajectoryFWDZ.push_back(temp.at(i).Z());
                rec.TrackTrajectoryFWDID.push_back(iTrackTraj);
            }

            temp = tracktraj.getBAKTrajectory();
            for(size_t i = 0; i < temp.size(); i++) {
                rec.TrackTrajectoryBWDX.push_back(temp.at(i).X());
                rec.TrackTrajectoryBWDY.push_back(temp.at(i).Y());
                rec.TrackTrajectoryBWDZ.push_back(temp.at(i).Z());
                rec.TrackTrajectoryBWDID.push_back(iTrackTraj);
            }
            iTrackTraj++;
        }
    }

    // save Vertex and Track-Vertex association info
    if (fWriteVertices) {
        size_t iVertex = 0;
        for ( auto const& vertex : (*VertexHandle) ) {
            rec.VertexIDNumber.push_back(vertex.getIDNumber());
            rec.VertexX.push_back(vertex.Position()[0]);
            rec.VertexY.push_back(vertex.Position()[1]);
            rec.VertexZ.push_back(vertex.Position()[2]);
            rec.VertexT.push_back(vertex.Time());

            int nVertexedTracks = 0;
            if ( findManyTrackEnd->isValid() ) {
                nVertexedTracks = findManyTrackEnd->at(iVertex).size();
            }
            rec.VertexN.push_back(nVertexedTracks);

            int vertexCharge = 0;
            for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack) {
                rec.VTAssn_VertIDNumber.push_back(vertex.getIDNumber());

                // Get this vertexed track.
                rec::Track track = *(findManyTrackEnd->at(iVertex).at(iVertexedTrack));
                rec.VTAssn_TrackIDNumber.push_back(track.getIDNumber());

                // Get the end of the track in the vertex.  It isn't that odd for the end
                // of the track not used in the vertex to be closer to the vertex than the
                // one actually used; you might have a very short stub track in a 3 track
                // vertex with small opening angles and the other 2 tracks might pull the
                // vertex some distance towards the far end of the stub track
                rec::TrackEnd fee = *(findManyTrackEnd->data(iVertex).at(iVertexedTrack));
                // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
                rec.VTAssn_TrackEnd.push_back(fee);

                if (fee==rec::TrackEndBeg) {
                    vertexCharge += track.ChargeBeg();
                } else {
                    vertexCharge += track.ChargeEnd();
                }
            }
            rec.VertexQ.push_back(vertexCharge);
            ++iVertex;
        } // end loop over VertexHandle
    }

    // save Vee and Track-Vee association info
    if (fWriteVees) {
        size_t iVee = 0;
        for ( auto const& vee : (*VeeHandle) ) {
            rec.VeeIDNumber.push_back(vee.getIDNumber());
            rec.VeeX.push_back(vee.Position()[0]);
            rec.VeeY.push_back(vee.Position()[1]);
            rec.VeeZ.push_back(vee.Position()[2]);
            rec.VeeT.push_back(vee.Time());
            rec.VeePXKpipi.push_back(vee.FourMomentum(0).X());
            rec.VeePYKpipi.push_back(vee.FourMomentum(0).Y());
            rec.VeePZKpipi.push_back(vee.FourMomentum(0).Z());
            rec.VeeEKpipi.push_back(vee.FourMomentum(0).E());
            rec.VeeMKpipi.push_back(vee.FourMomentum(0).M());
            rec.VeePXLppi.push_back(vee.FourMomentum(1).X());
            rec.VeePYLppi.push_back(vee.FourMomentum(1).Y());
            rec.VeePZLppi.push_back(vee.FourMomentum(1).Z());
            rec.VeeELppi.push_back(vee.FourMomentum(1).E());
            rec.VeeMLppi.push_back(vee.FourMomentum(1).M());
            rec.VeePXLpip.push_back(vee.FourMomentum(2).X());
            rec.VeePYLpip.push_back(vee.FourMomentum(2).Y());
            rec.VeePZLpip.push_back(vee.FourMomentum(2).Z());
            rec.VeeELpip.push_back(vee.FourMomentum(2).E());
            rec.VeeMLpip.push_back(vee.FourMomentum(2).M());

            int nVeeTracks = 0;
            if ( findManyVeeTrackEnd->isValid() ) {
                nVeeTracks = findManyVeeTrackEnd->at(iVee).size();
            }

            for (int iVeeTrack=0; iVeeTrack<nVeeTracks; ++iVeeTrack) {
                rec.VeeTAssn_VeeIDNumber.push_back(vee.getIDNumber());

                // Get this vertexed track.
                rec::Track track = *(findManyVeeTrackEnd->at(iVee).at(iVeeTrack));
                rec.VeeTAssn_TrackIDNumber.push_back(track.getIDNumber());

                rec::TrackEnd fee = *(findManyVeeTrackEnd->data(iVee).at(iVeeTrack));
                // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
                rec.VeeTAssn_TrackEnd.push_back(fee);
            }
            ++iVee;
        } // end loop over VeeHandle
    }

    // save Cluster info
    if (fWriteCaloClusters) {
        for ( auto const& cluster : (*RecoClusterHandle) ) {
            rec.nCluster++;
            rec.ClusterIDNumber.push_back(cluster.getIDNumber());
            rec.ClusterNhits.push_back(cluster.CalorimeterHits().size());
            rec.ClusterEnergy.push_back(cluster.Energy());
            rec.ClusterTime.push_back(cluster.Time());
            rec.ClusterTimeDiffFirstLast.push_back(cluster.TimeDiffFirstLast());
            rec.ClusterX.push_back(cluster.Position()[0]);
            rec.ClusterY.push_back(cluster.Position()[1]);
            rec.ClusterZ.push_back(cluster.Position()[2]);
            rec.ClusterTheta.push_back(cluster.ITheta());
            rec.ClusterPhi.push_back(cluster.IPhi());
            rec.ClusterPID.push_back(cluster.ParticleID());
            // rec.ClusterShape.push_back(cluster.Shape());
            rec.ClusterMainAxisX.push_back(cluster.EigenVectors()[0]);
            rec.ClusterMainAxisY.push_back(cluster.EigenVectors()[1]);
            rec.ClusterMainAxisZ.push_back(cluster.EigenVectors()[2]);

            // Matching MCParticle info
            std::vector<std::pair<simb::MCParticle*,float>> trakt;
            trakt = BackTrack->ClusterToMCParticles( const_cast<rec::Cluster*>(&cluster) );
            int eileen = -1;
            if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
                eileen = TrackIdToIndex[trakt[0].first->TrackId()];
            }
            rec.ClusterMCindex.push_back(eileen);
            if (eileen > -1) {
            	rec.ClusterMCfrac.push_back(trakt[0].second);
            } else {
            	rec.ClusterMCfrac.push_back(0.0);
            }
        }
    }

    // Write info for ECAL-matched tracks
    if (fWriteMatchedTracks) {
        size_t iCluster = 0;
        for ( auto const& cluster : (*RecoClusterHandle) ) {
            int nCALedTracks(0);
            if ( findManyCALTrackEnd->isValid() ) {
                nCALedTracks = findManyCALTrackEnd->at(iCluster).size();
            }
            for (int iCALedTrack=0; iCALedTrack<nCALedTracks; ++iCALedTrack) {
                rec.CALAssn_ClusIDNumber.push_back(cluster.getIDNumber());
                rec::Track track  = *(findManyCALTrackEnd->at(iCluster).at(iCALedTrack));
                rec.CALAssn_TrackIDNumber.push_back( track.getIDNumber() );

                rec::TrackEnd fee = *(findManyCALTrackEnd->data(iCluster).at(iCALedTrack));
                rec.CALAssn_TrackEnd.push_back(fee);    // The rec::TrackEnd (see Track.h) that extrapolated to cluster
            }
            iCluster++;
        }
    } // end branch on fWriteCaloInfo

    return;
} // end :anatree::FillVectors



//==============================================================================
//==============================================================================
//==============================================================================
// Process ionization.  Eventually this moves into the reco code.
void gar::anatree::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
float& forwardIonVal, float& backwardIonVal) {

    // NO CALIBRATION SERVICE FOR NOW

    std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
    forwardIonVal = processOneDirection(SigData, ionizeTruncate);

    SigData = ion.getBAK_dSigdXs();
    backwardIonVal = processOneDirection(SigData, ionizeTruncate);

    return;
}



float gar::anatree::processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

    std::vector<std::pair<float,float>> dEvsX;    // Ionization vs distance along track

    // The first hit on the track never had its ionization info stored.  Not a problem
    // really.  Each pair is a hit and the step along the track that ends at the hit
    // For the last hit, just take the step from the n-1 hit; don't guess some distance
    // to (nonexistant!) n+1 hit.  Using pointer arithmetic because you are a real K&R
    // C nerd!  Except that C++ doesn't know you are such a nerd and if
    //  SigData.size()==0, then SigData.end()-1 is 0xFFFFFFFFFFFFFFF8.
    if (SigData.size()==0) return 0.0;
    float distAlongTrack = 0;
    std::vector<std::pair<float,float>>::iterator littlebit = SigData.begin();
    for (; littlebit<(SigData.end()-1); ++littlebit) {
        float dE =   std::get<0>(*littlebit);
        // tpctrackfit2_module.cc fills the TrackIoniz data product so that
        // this quantity is really dL > 0 not dX, a coordinate on the drift axis
        float dX  = std::get<1>(*littlebit);
        distAlongTrack += dX;    // But count full step to get hit position on track
        // Take dX to be 1/2 the previous + last segment
        dX += std::get<1>(*(littlebit+1));
        float dEdX = dE/(0.5*dX);

        std::pair pushme = std::make_pair(dEdX,distAlongTrack);
        dEvsX.push_back( pushme );
    }

    // Get the truncated mean; first sort then take mean
    std::sort(dEvsX.begin(),dEvsX.end(), lessThan_byE);

    // Get the dEdX vs length data, truncated.
    int goUpTo = ionizeTruncate * dEvsX.size() +0.5;
    if (goUpTo > (int)dEvsX.size()) goUpTo = dEvsX.size();
    int i = 1;        float returnvalue = 0;
    littlebit = dEvsX.begin();
    for (; littlebit<dEvsX.end(); ++littlebit) {
        returnvalue += std::get<0>(*littlebit);
        ++i;
        if (i>goUpTo) break;
    }
    returnvalue /= goUpTo;
    return returnvalue;
}



//==============================================================================
//==============================================================================
//==============================================================================
// Coherent pion analysis specific code
float gar::anatree::computeT( simb::MCTruth theMCTruth ) {
    // Warning.  You probably want the absolute value of t, not t.
    int nPart = theMCTruth.NParticles();
    enum { nu, mu, pi};
    float E[3], Px[3], Py[3], Pz[3];
    E[nu] = E[mu] = E[pi] = -1e42;

    for (int i=0; i<3;++i) {
        Px[i] = 0;
        Py[i] = 0;
        Pz[i] = 0;
        E[i]  = 0;
    }
    // Find t from the MCParticles via the
    for (int iPart=0; iPart<nPart; iPart++) {
        simb::MCParticle Part = theMCTruth.GetParticle(iPart);
        int code = Part.PdgCode();
        int mom  = Part.Mother();

        // get the neutrino
        if ( abs(code) == 12 || abs(code) == 14 || abs(code) == 16 ) {
            if (mom == -1) {
                E[nu] = Part.E();   Px[nu] = Part.Px();   Py[nu] = Part.Py();   Pz[nu] = Part.Pz();
            }
        }

        // get the lepton
        if ( abs(code) == 11 || abs(code) == 13 || abs(code) == 15 ) {
            if (mom == 0) {
                E[mu] = Part.E();   Px[mu] = Part.Px();   Py[mu] = Part.Py();   Pz[mu] = Part.Pz();
            }
        }

        // get the pion
        if ( code==111 || abs(code)==211 ) {
            if (mom == 1) {
                E[pi] = Part.E();   Px[pi] = Part.Px();   Py[pi] = Part.Py();   Pz[pi] = Part.Pz();
            }
        }

        // get outa here
        if ( E[nu]!=0 && E[mu]!=0 && E[pi]!=0) break;

    }

    // Compute t; reuse nu 4-vector to get first q, then t.
    E[nu] -= E[mu];   Px[nu] -= Px[mu];   Py[nu] -= Py[mu];   Pz[nu] -= Pz[mu];
    E[nu] -= E[pi];   Px[nu] -= Px[pi];   Py[nu] -= Py[pi];   Pz[nu] -= Pz[pi];
    float t = E[nu]*E[nu] -Px[nu]*Px[nu] -Py[nu]*Py[nu] -Pz[nu]*Pz[nu];
    return t;
}



//==============================================================================
//==============================================================================
//==============================================================================
std::vector< std::pair<int, float> > gar::anatree::processPIDInfo( float p ) {

    std::vector<std::string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};
    std::vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
    std::vector< std::pair<int, float> > pid;
    pid.resize(6);

    int qclosest = 0;
    float dist = 100000000.;
    CLHEP::RandFlat FlatRand(fEngine);

    for (int q = 0; q < 501; ++q) {
        //Check the title and the reco momentum take only the one that fits
        std::string fulltitle = m_pidinterp[q]->GetTitle();
        unsigned first = fulltitle.find("=");
        unsigned last = fulltitle.find("GeV");
        std::string substr = fulltitle.substr(first+1, last - first-1);
        float pidinterp_mom = std::atof(substr.c_str());
        //calculate the distance between the bin and mom, store the q the closest
        float disttemp = std::abs(pidinterp_mom - p);

        if( disttemp < dist ) {
            dist = disttemp;
            qclosest = q;
        }
    } // closes the "pidmatrix" loop

    //Compute all the probabities for each type of true to reco
    //loop over the columns (true pid)
    for (int pidm = 0; pidm < 6; ++pidm) {

        //loop over the columns (true pid)
        std::vector< std::pair<float, std::string> > v_prob;

        //loop over the rows (reco pid)
        for (int pidr = 0; pidr < 6; ++pidr) {
            std::string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
            float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1, pidr+1);
            //Need to check random number value and prob value then associate the recopdg to the reco prob
            v_prob.push_back( std::make_pair(prob, recoparticlename) );
        }

        //Compute the pid from it
        if (v_prob.size() > 1) {
            //Order the vector of prob
            std::sort(v_prob.begin(), v_prob.end());
            //Throw a random number between 0 and 1
            float random_number = FlatRand.fire();
            //Make cumulative sum to get the range
            std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});

            for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++) {
                if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first ) {
                    pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ), v_prob.at(ivec+1).first) );
                }
            }
        } else {
            pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ), v_prob.at(0).first) );
        }
    }

    //return a vector of pid and prob
    return pid;
}



DEFINE_ART_MODULE(gar::anatree)

