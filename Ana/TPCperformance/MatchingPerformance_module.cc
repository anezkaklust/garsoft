////////////////////////////////////////////////////////////////////////////////
// Class:        MatchingPerformance
// Plugin Type: analyzer (art v2_11_02)
// File:        MatchingPerformance_module.cc
//
// Generated 25 Sep 2020 by Leo Bellantoni
// 
// For each reco track:
//     Use backtracker to get matching MCparticles.  Take charged & stable
//         matching MCParticle with largest energy fraction for the track
//     Only use MCParticles which go eventually into the ECAL?
//     Use backtracker info to determine "starting" end of the track.
//     Record chi-2 and NTPC clusters in the track (and other stuff)
//     List all the clusters that match back to the MCParticle, and all
//         that match to the track.  Tabulate the number and energy sum
//             For clusters on track & not MCParticle
//             For clusters on MCParticle & not on track
//             For clusters on both
//
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "Math/ProbFuncMathCore.h"

#include <string>
#include <vector>
#include <unordered_map>





namespace gar {

    class MatchingPerformance : public art::EDAnalyzer {
    public:
        explicit MatchingPerformance(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        MatchingPerformance(MatchingPerformance const &) = delete;
        MatchingPerformance(MatchingPerformance &&) = delete;
        MatchingPerformance & operator = (MatchingPerformance const &) = delete;
        MatchingPerformance & operator = (MatchingPerformance &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;



    private:
        void ClearVectors();
        void FillVectors(art::Event const & e);



        // Position of TPC from geometry service; 1 S Boston Ave.
        float ItsInTulsa[3];

        // Input data labels
        std::string fGeantLabel;
        std::string fTrackLabel;
        std::string fClusterLabel;
        std::string fECALAssnLabel;

        // the analysis tree
        TTree *fTree;

        // Geometry service
        const gar::geo::GeometryCore* fGeo;
        // Backtracker service
        cheat::BackTrackerCore* fBack;
        // PDG database from ROOT
        TDatabasePDG* pdgInstance;



        // The analysis tree. Use Rtypes.h here, as these data get used by root
        // global event info
        Int_t    fRun;
        Int_t    fSubRun;
        Int_t    fEvent;

        // This analysis starts from reco'd tracks
        std::vector<ULong64_t>          fTrackIDNumber;
        std::vector<Float_t>            fTrackX;
        std::vector<Float_t>            fTrackY;
        std::vector<Float_t>            fTrackZ;
        std::vector<Float_t>            fTrackPX;
        std::vector<Float_t>            fTrackPY;
        std::vector<Float_t>            fTrackPZ;
        std::vector<Float_t>            fTrackEend;
        std::vector<Int_t>              fTrackQ;
        std::vector<Float_t>            fTrackLen;
        std::vector<Float_t>            fTrackChi2;
        std::vector<Int_t>              fNTPCClustersOnTrack;
        std::vector<Float_t>            fTrackPval;
        std::vector<Double_t>           fTrackTime;

        // Which match to qualified MCParticles
        std::vector<Int_t>              fMCPDG;
        std::vector<Int_t>              fMCPDGMother;
        std::vector<Float_t>            fMCPStartX;
        std::vector<Float_t>            fMCPStartY;
        std::vector<Float_t>            fMCPStartZ;
        std::vector<Float_t>            fMCPStartPX;
        std::vector<Float_t>            fMCPStartPY;
        std::vector<Float_t>            fMCPStartPZ;
        std::vector<Float_t>            fMCPStartE;
        std::vector<std::string>        fMCPProc;
        std::vector<std::string>        fMCPEndProc;
        std::vector<Float_t>            fMCPTime;

        // And here's the track - ECAL matching info
        std::vector<Float_t>            fEassocNoInTruth;
        std::vector<Int_t>              fNassocNoInTruth;
        std::vector<Float_t>            fEunassocInTruth;
        std::vector<Int_t>              fNunassocInTruth;
        std::vector<Float_t>            fEisassocInTruth;
        std::vector<Int_t>              fNisassocInTruth;

    };
}



//==============================================================================
//==============================================================================
//==============================================================================
gar::MatchingPerformance::MatchingPerformance(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    fGeantLabel    = p.get<std::string>("GEANTLabel",   "geant");
    fTrackLabel    = p.get<std::string>("TrackLabel",   "track");
    fClusterLabel  = p.get<std::string>("ClusterLabel", "calocluster");
    fECALAssnLabel = p.get<std::string>("ECALAssnLabel","trkecalassn");

    pdgInstance = TDatabasePDG::Instance();

    consumesMany<std::vector<simb::MCTruth> >();
    consumes<std::vector<simb::MCParticle> >(fGeantLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);
    consumes<std::vector<rec::Cluster> >(fClusterLabel);
    consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

    return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::beginJob() {

    fGeo = gar::providerFrom<geo::GeometryGAr>();
    ItsInTulsa[0] = fGeo->TPCXCent();        // 1 S Boston Ave
    ItsInTulsa[1] = fGeo->TPCYCent();
    ItsInTulsa[2] = fGeo->TPCZCent();

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");



    fTree->Branch("Run",                    &fRun,              "Run/I");
    fTree->Branch("SubRun",                 &fSubRun,           "SubRun/I");
    fTree->Branch("Event",                  &fEvent,            "Event/I");

    fTree->Branch("TrackIDNumber",          &fTrackIDNumber);
    fTree->Branch("TrackX",                 &fTrackX);
    fTree->Branch("TrackY",                 &fTrackY);
    fTree->Branch("TrackZ",                 &fTrackZ);
    fTree->Branch("TrackPX",                &fTrackPX);
    fTree->Branch("TrackPY",                &fTrackPY);
    fTree->Branch("TrackPZ",                &fTrackPZ);
    fTree->Branch("TrackEend",              &fTrackEend);
    fTree->Branch("TrackQ",                 &fTrackQ);
    fTree->Branch("TrackLen",               &fTrackLen);
    fTree->Branch("TrackChi2",              &fTrackChi2);
    fTree->Branch("NTPCClustersOnTrack",    &fNTPCClustersOnTrack);
    fTree->Branch("TrackPval",              &fTrackPval);
    fTree->Branch("TrackTime",              &fTrackTime);

    fTree->Branch("PDG",                    &fMCPDG);
    fTree->Branch("PDGMother",              &fMCPDGMother);
    fTree->Branch("MCPStartX",              &fMCPStartX);
    fTree->Branch("MCPStartY",              &fMCPStartY);
    fTree->Branch("MCPStartZ",              &fMCPStartZ);
    fTree->Branch("MCPStartPX",             &fMCPStartPX);
    fTree->Branch("MCPStartPY",             &fMCPStartPY);
    fTree->Branch("MCPStartPZ",             &fMCPStartPZ);
    fTree->Branch("MCPStartE",              &fMCPStartE);
    fTree->Branch("MCPProc",                &fMCPProc);
    fTree->Branch("MCPEndProc",             &fMCPEndProc);
    fTree->Branch("MCPTime",                &fMCPTime);

    fTree->Branch("EassocNoInTruth",        &fEassocNoInTruth);
    fTree->Branch("NassocNoInTruth",        &fNassocNoInTruth);
    fTree->Branch("EunassocInTruth",        &fEunassocInTruth);
    fTree->Branch("NunassocInTruth",        &fNunassocInTruth);
    fTree->Branch("EisassocInTruth",        &fEisassocInTruth);
    fTree->Branch("NisassocInTruth",        &fNisassocInTruth);

    return;
}  // End of :MatchingPerformance::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::analyze(art::Event const & event) {

    // Is the backtracker good to go?  Need to fill fBack on per-event basis
    // not in beginJob()
    cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
    fBack = const_cast<cheat::BackTrackerCore*>(const_bt);
    if ( !fBack->HasTracks() || !fBack->HasClusters() ) {
        throw cet::exception("MatchingPerformance") << " Not enough backtracker info"
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    ClearVectors();
    FillVectors(event);
    fTree->Fill();

    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::ClearVectors() {
    // clear out all our vectors

    fTrackIDNumber.clear();
    fTrackX.clear();
    fTrackY.clear();
    fTrackZ.clear();
    fTrackPX.clear();
    fTrackPY.clear();
    fTrackPZ.clear();
    fTrackEend.clear();
    fTrackQ.clear();
    fTrackLen.clear();
    fTrackChi2.clear();
    fNTPCClustersOnTrack.clear();
    fTrackPval.clear();
    fTrackTime.clear();

    fMCPDG.clear();
    fMCPDGMother.clear();
    fMCPStartX.clear();
    fMCPStartY.clear();
    fMCPStartZ.clear();
    fMCPStartPX.clear();
    fMCPStartPY.clear();
    fMCPStartPZ.clear();
    fMCPStartE.clear();
    fMCPProc.clear();
    fMCPEndProc.clear();
    fMCPTime.clear();

    fEassocNoInTruth.clear();
    fNassocNoInTruth.clear();
    fEunassocInTruth.clear();
    fNunassocInTruth.clear();
    fEisassocInTruth.clear();
    fNisassocInTruth.clear();

    return;
} // end :MatchingPerformance::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::FillVectors(art::Event const& event) {

    // =============  Get art handles ==========================================
    // Get handles for Tracks, clusters, associations between them
    art::Handle< std::vector<rec::Track> > TrackHandle;
    if (!event.getByLabel(fTrackLabel, TrackHandle)) {
        throw cet::exception("MatchingPerformance") << " No rec::Track"
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    art::Handle< std::vector<rec::Cluster> > ClusterHandle;
    if (!event.getByLabel(fClusterLabel, ClusterHandle)) {
        throw cet::exception("MatchingPerformance") << " No rec::Cluster branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    art::FindMany<rec::Cluster, rec::TrackEnd>* findManyTrackEndCAL = NULL;
    findManyTrackEndCAL = new art::FindMany<rec::Cluster, rec::TrackEnd>
            (TrackHandle,event,fECALAssnLabel);
    if ( !findManyTrackEndCAL->isValid() ) {
        throw cet::exception("MatchingPerformance") << " Bad TrackEnd-ECAL Assn."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
 
    // Get handles for MCinfo, also good for MCPTrajectory.  Want to get all
    // MCTruths, regardless of generator label.
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mctruthHandles;
    art::Handle< std::vector<simb::MCParticle> > MCPHandle;
    event.getManyByType(mctruthHandles);
    if (mctruthHandles.size()!=1) {
        throw cet::exception("MatchingPerformance") << " Need just 1 simb::MCTruth"
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (!event.getByLabel(fGeantLabel, MCPHandle)) {
        throw cet::exception("MatchingPerformance") << " No simb::MCParticle"
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Need something to look at!
    if ( TrackHandle->size()==0 || ClusterHandle->size()==0 ) return;

    // Make the usual map for the MC info
    typedef int TrkId;
    std::unordered_map<TrkId, Int_t> TrackIdToIndex;
    Int_t index = 0;
    for ( auto const& mcp : (*MCPHandle) ) {
        int TrackId = mcp.TrackId();
        TrackIdToIndex[TrackId] = index++;
    }

    // Will need a vector of art::Ptrs to all the ECAL clusters for the
    // backtracker calls.
    std::vector<art::Ptr<rec::Cluster>> ClusterGrabber;
    art::PtrMaker<rec::Cluster> makeClusterPtr(event,ClusterHandle.id());
    for (size_t iCluster=0; iCluster<ClusterHandle->size(); ++iCluster) {
        art::Ptr<rec::Cluster> aPtr = makeClusterPtr(iCluster);
        ClusterGrabber.push_back(aPtr);
    }
    // Will need a vector<art::Ptr<rec::Track>> of tracks also
    std::vector<art::Ptr<rec::Track>> TrackGrabber;
    art::PtrMaker<rec::Track> makeTrackPtr(event,TrackHandle.id());
    for (size_t iTrack=0; iTrack<TrackHandle->size(); ++iTrack) {
        art::Ptr<rec::Track> aPtr = makeTrackPtr(iTrack);
        TrackGrabber.push_back(aPtr);
    }





    // =============  Pull art's handles =======================================
    fRun    = event.run();
    fSubRun = event.subRun();
    fEvent  = event.id().event();



    // Try to kill stub or other undesirable tracks first
    std::vector<std::pair<rec::Track,simb::MCParticle>> niceTracks;
    for ( auto const& track : (*TrackHandle) ) {

        /// Want to be able to know what you've matched back to at MCTruth
        std::vector<std::pair<simb::MCParticle*,float>> MCsfromTrack;
        gar::rec::Track* nonconst_track = const_cast<gar::rec::Track*>(&track);
        MCsfromTrack = fBack->TrackToMCParticles(nonconst_track);
        int nMCsfromTrack = MCsfromTrack.size();
        if (nMCsfromTrack==0) continue;

        simb::MCParticle theMCPart;
        bool broke = false;
        for (int iMCfromTrack =0; iMCfromTrack<nMCsfromTrack; ++iMCfromTrack) {
            // Plausible MCParticle to make this track?
            theMCPart = *(MCsfromTrack[iMCfromTrack].first);
            TParticlePDG* particle = pdgInstance->GetParticle(theMCPart.PdgCode());
            if (particle==NULL         ) continue;   // What causes this err? If anything.
            if (particle->Charge()==0.0) continue;
            if (particle->Stable()==0  ) continue;
            broke = true;
            break;            // TrackToMCParticles returns MCParts sorted by how much
                              // the MCParts contribute to the track; take highest fraction.
        }
        if (!broke) continue;

        // Does theMCParticle go out to the ECAL?
        int nTraj = theMCPart.NumberTrajectoryPoints();
        TVector3 doink;            bool whackThatECAL = false;
        for (int iTraj=0; iTraj<nTraj; ++iTraj) {
            doink.SetXYZ(theMCPart.Vx(iTraj),theMCPart.Vy(iTraj),theMCPart.Vz(iTraj));
            if ( fGeo->PointInECALBarrel(doink) ) {
                whackThatECAL = true;
                break;
            }
        }
        if (!whackThatECAL) continue;

        // This track is nice so far; have to check for stubbiness
        niceTracks.push_back(std::make_pair(track,theMCPart));
    }
 
    for (size_t iNiceTrk=0; iNiceTrk<niceTracks.size(); ++iNiceTrk) {
        simb::MCParticle theMCPart = niceTracks[iNiceTrk].second;

        // Look at all the other tracks that have contributions from this MCParticle
        // & make sure this is the "biggest" one in terms of deposited charge; take the
        // smaller-contributing (stub) tracks out of consideration
        std::vector<art::Ptr<rec::Track>> tracksFromSameMCP;
        tracksFromSameMCP = fBack->MCParticleToTracks(&theMCPart,TrackGrabber);

        std::vector<std::pair<ULong64_t,float>> trkID_Epairs;
        for (size_t iTrkFromSame=0; iTrkFromSame<tracksFromSameMCP.size(); ++iTrkFromSame) {
            // index & dereference art::Ptr.
            rec::Track tmpTrk = *(tracksFromSameMCP[iTrkFromSame]);
            std::vector<art::Ptr<rec::Hit>> hitList = fBack->TrackToHits( &tmpTrk );
            float thisTracksE = 0;
            for ( auto& ptrHit : hitList ) thisTracksE += ptrHit->Signal();
            ULong64_t thisTracksID = tmpTrk.getIDNumber();
            trkID_Epairs.push_back( std::make_pair(thisTracksID,thisTracksE) );
        }
        // Next, sort trkID_Epairs by value.  1st, declare the type of the sorting predicate
        typedef std::function<bool(std::pair<ULong64_t,float>, std::pair<ULong64_t,float>)>
            Comparator;
        // Create a set to contain the pairs sorted using lambda func
        std::set<std::pair<ULong64_t,float>, Comparator> sortedTrkID_Epairs(
            trkID_Epairs.begin(), trkID_Epairs.end(),
            [](std::pair<ULong64_t,float> a ,std::pair<ULong64_t,float> b) {
                return a.second > b.second;
            }
        );
        // Keep the first one; the rest are stubs to remove.
        if (sortedTrkID_Epairs.size()>0) sortedTrkID_Epairs.erase(sortedTrkID_Epairs.begin());
        for ( auto& iTrkID_E : sortedTrkID_Epairs ) {
            std::vector<std::pair<rec::Track,simb::MCParticle>>::iterator notNice;
            notNice = niceTracks.begin();
            for (; notNice<niceTracks.end(); ++notNice) {
                if ( notNice->first.getIDNumber() == iTrkID_E.first ) {
                    // notNice points to a (Track,MCParticle) pair where the Track
                    // is a stub
                    notNice = niceTracks.erase(notNice);
                }
            }
        }
    }



    for (size_t iNiceTrk=0; iNiceTrk<niceTracks.size(); ++iNiceTrk) {
        rec::Track       track     = niceTracks[iNiceTrk].first;
        simb::MCParticle theMCPart = niceTracks[iNiceTrk].second;

        // Which end of this track corresponds to the start of the reco track?
        // This calc in the geometry's coordinates
        const TLorentzVector& positionMCP = theMCPart.Position(0);
        float distStart = std::hypot(track.Vertex()[1] -positionMCP[1],
                                     track.Vertex()[2] -positionMCP[2]);
        float distEnd   = std::hypot(track.End()[1]    -positionMCP[1],
                                     track.End()[2]    -positionMCP[2]);
        rec::TrackEnd endTowardECAL = (distStart>distEnd) ? rec::TrackEndBeg
                                                          : rec::TrackEndEnd;



        // Record some info about the track at the end near the ECAL
        fTrackIDNumber.push_back(track.getIDNumber());
        Float_t saveChi2;		
		TParticlePDG* particle = pdgInstance->GetParticle(theMCPart.PdgCode());
        float MCmass = particle->Mass();
        if ( endTowardECAL==rec::TrackEndBeg ) {
            fTrackX.push_back   ( track.Vertex()[0] -ItsInTulsa[0]);
            fTrackY.push_back   ( track.Vertex()[1] -ItsInTulsa[1]);
            fTrackZ.push_back   ( track.Vertex()[2] -ItsInTulsa[2]);
            fTrackPX.push_back  (-track.Momentum_beg()*track.VtxDir()[0]);
            fTrackPY.push_back  (-track.Momentum_beg()*track.VtxDir()[1]);
            fTrackPZ.push_back  (-track.Momentum_beg()*track.VtxDir()[2]);
            fTrackEend.push_back( std::hypot(track.Momentum_beg(),MCmass));
            fTrackQ.push_back   ( track.ChargeEnd());
            fTrackLen.push_back ( track.LengthBackward());
            saveChi2 = track.ChisqBackward();
            fTrackChi2.push_back( saveChi2);
        } else {
            fTrackX.push_back   ( track.End()[0] -ItsInTulsa[0]);
            fTrackY.push_back   ( track.End()[1] -ItsInTulsa[1]);
            fTrackZ.push_back   ( track.End()[2] -ItsInTulsa[2]);
            fTrackPX.push_back  (-track.Momentum_end()*track.EndDir()[0]);
            fTrackPY.push_back  (-track.Momentum_end()*track.EndDir()[1]);
            fTrackPZ.push_back  (-track.Momentum_end()*track.EndDir()[2]);
            fTrackEend.push_back( std::hypot(track.Momentum_end(),MCmass));
            fTrackQ.push_back   ( track.ChargeBeg());
            fTrackLen.push_back ( track.LengthForward());        
            saveChi2 = track.ChisqForward();
            fTrackChi2.push_back( saveChi2);
        }
        Int_t nHits = track.NHits();
        fNTPCClustersOnTrack.push_back(nHits);
        Float_t pVal = ROOT::Math::chisquared_cdf_c(saveChi2,nHits);
        fTrackPval.push_back(pVal);
        fTrackTime.push_back(track.Time());

        // Record some info about the matching MC
        fMCPDG.push_back(theMCPart.PdgCode());
        TrkId momTrkId = theMCPart.Mother();
        int momPDG = 0;
        if (momTrkId>0) {
            int momIndex = TrackIdToIndex[momTrkId];
            momPDG   = (*MCPHandle).at(momIndex).PdgCode();
        }
        fMCPDGMother.push_back(momPDG);

        fMCPStartX.push_back(positionMCP.X());
        fMCPStartY.push_back(positionMCP.Y());
        fMCPStartZ.push_back(positionMCP.Z());
        const TLorentzVector& momentumMCP = theMCPart.Momentum(0);
        fMCPStartPX.push_back(momentumMCP.Px());
        fMCPStartPY.push_back(momentumMCP.Py());
        fMCPStartPZ.push_back(momentumMCP.Pz());
        fMCPStartE.push_back(theMCPart.E());

        fMCPProc.push_back(theMCPart.Process());
        fMCPEndProc.push_back(theMCPart.EndProcess());
        fMCPTime.push_back(theMCPart.T());    



        // Finally, look for the ECAL-track matching info
        std::vector<art::Ptr<rec::Cluster>> clustersOnMCP;
        clustersOnMCP = fBack->MCParticleToClusters(&theMCPart, ClusterGrabber);

        float eAssocNoInTruth = 0.0;            int nAssocNoInTruth = 0;
        float eUnassocInTruth = 0.0;            int nUnassocInTruth = 0;
        float eIsassocInTruth = 0.0;            int nIsassocInTruth = 0;
        size_t   nCALedClusts =   0;
        std::vector<rec::Cluster  const*> clustersFromAssn;
        std::vector<rec::TrackEnd const*> trackEndsFromAssn;

        // The art::FindMany object findManyTrackEndCAL is indexed over 
        // all the tracks in *TrackHandle (which is bigger than niceTracks).
        int iTrack = -1;
        for (size_t iAllTrk=0; iAllTrk<TrackHandle->size(); ++iAllTrk) {
            iTrack++;
            if ( (*TrackHandle)[iAllTrk] == track ) break;
        }
        findManyTrackEndCAL->get(iTrack, clustersFromAssn,trackEndsFromAssn);
        nCALedClusts = clustersFromAssn.size();

        for ( auto const& cluster : (*ClusterHandle) ) {

            bool clusterOnTrack = false;
            for (size_t iCALedClust=0; iCALedClust<nCALedClusts; ++iCALedClust) {
                if (*trackEndsFromAssn[iCALedClust] != endTowardECAL) continue;
                if (*clustersFromAssn[iCALedClust] == cluster) {
                    clusterOnTrack = true;
                    break;
                }
            }

            bool clusterOnMCP = false;
            for (size_t iClusOnMCP = 0; iClusOnMCP<clustersOnMCP.size(); ++iClusOnMCP) {
                if ( *(clustersOnMCP[iClusOnMCP]) == cluster ) {
                    clusterOnMCP = true;
                    break;
                }
            }

            if ( clusterOnTrack &&!clusterOnMCP ) {
                eAssocNoInTruth += cluster.Energy();
                nAssocNoInTruth += cluster.CalorimeterHits().size();
            }
            if (!clusterOnTrack && clusterOnMCP ) {
                eUnassocInTruth += cluster.Energy();
                nUnassocInTruth += cluster.CalorimeterHits().size();
            }
            if ( clusterOnTrack && clusterOnMCP ) {
                eIsassocInTruth += cluster.Energy();
                nIsassocInTruth += cluster.CalorimeterHits().size();
            }
        }

        fEassocNoInTruth.push_back(eAssocNoInTruth);
        fNassocNoInTruth.push_back(nAssocNoInTruth);
        fEunassocInTruth.push_back(eUnassocInTruth);
        fNunassocInTruth.push_back(nUnassocInTruth);
        fEisassocInTruth.push_back(eIsassocInTruth);
        fNisassocInTruth.push_back(nIsassocInTruth);
    
    } // end loop over TrackHandle





    return;
} // end MatchingPerformance::FillVectors





DEFINE_ART_MODULE(gar::MatchingPerformance)
