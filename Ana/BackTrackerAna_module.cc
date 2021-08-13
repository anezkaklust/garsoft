////////////////////////////////////////////////////////////////////////
// Class:       BackTrackerAna
// Plugin Type: analyzer (art v3_06_03)
// File:        BackTrackerAna_module.cc
//
// Generated at Tue Jul 20 09:48:51 2021 by Christopher Hilgenberg using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

// framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// simb includes
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// garsoft includes
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/Vee.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "CoreUtils/ServiceUtil.h"
#include "MCCheater/BackTracker.h"
#include "Geometry/GeometryGAr.h"

// root includes
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

// c++ includes
#include <vector>
#include <string>
#include <map>
#include <iostream>

using namespace gar;
using namespace std;
using art::Assns;

class BackTrackerAna;


class BackTrackerAna : public art::EDAnalyzer {
public:
  explicit BackTrackerAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BackTrackerAna(BackTrackerAna const&) = delete;
  BackTrackerAna(BackTrackerAna&&) = delete;
  BackTrackerAna& operator=(BackTrackerAna const&) = delete;
  BackTrackerAna& operator=(BackTrackerAna&&) = delete;

  virtual void beginJob() override;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // input product labels
  string fGeantLabel; ///< module label for geant4 simulated hits
  string fGeantInstanceCalo; ///< Instance name ECAL for sdp::CaloDeposit
  string fHitLabel; ///< module label for reco TPC hits rec::Hit
  string fTPCClusterLabel; ///< module label for TPC Clusters rec::TPCCluster
  string fTrackLabel; ///< module label for TPC Tracks rec:Track
  string fVertexLabel; ///< module label for vertexes rec:Vertex
  string fVeeLabel; ///< module label for conversion/decay vertexes rec:Vee
  string fClusterLabel; ///< module label for calo clusters rec::Cluster
  string fCaloHitLabel; ///< module label for reco calo hits rec::CaloHit
  string fCaloHitInstanceCalo; ///< Instance name ECAL for rec::CaloHit

  // ana objects
  //TTree *fTree;
  TH1F  *fHistGenDiff;
  TH1F  *fHistG4Diff;
  TH1F  *fHistHitComplete;
  TH1F  *fHistHitEnergy;
  TH1F  *fHistTrackHitDiff;
  TH1F  *fHistTrackHitMult;
  TH1F  *fHistClustHitMult;
  TH1F  *fHistClustMult;
  TH1F  *fHistNHitAssocTrackNHitTotDiff;
  TH1F  *fHistEdepRelKE;


  TH2F  *fHistG4Diff2d;

  // Geometry
  const geo::GeometryCore* fGeo; ///< pointer to the geometry service

  // backtracker for determining inter product associations
  cheat::BackTrackerCore* fBt;

}; // class header


BackTrackerAna::BackTrackerAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
    fGeantLabel        = p.get<string>("GEANTLabel","geant");
    fGeantInstanceCalo = p.get<string>("GEANTInstanceCalo","ECAL");
    fHitLabel         = p.get<string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<string>("TrackLabel","track");
    fVertexLabel      = p.get<string>("VertexLabel","vertex");
    fVeeLabel         = p.get<string>("VeeLabel","veefinder1");
    fCaloHitLabel        = p.get<string>("CaloHitLabel","sipmhit");
    fCaloHitInstanceCalo = p.get<string>("CaloHitInstanceCalo","ECAL");
    fClusterLabel     = p.get<string>("ClusterLabel","calocluster");

    // gen
    consumesMany<vector<simb::MCTruth> >();
    consumesMany<vector<simb::GTruth> >();

    // g4
    consumes<Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
    consumes<vector<simb::MCParticle> >(fGeantLabel);
    consumes<vector<sdp::EnergyDeposit> >(fGeantLabel);
    art::InputTag ecalgeanttag(fGeantLabel, fGeantInstanceCalo);
    consumes<vector<sdp::CaloDeposit> >(ecalgeanttag);

    // reco
    consumes<vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<vector<rec::Hit> >(fHitLabel);
    consumes<vector<rec::Track> >(fTrackLabel);
    consumes<vector<rec::Vertex> >(fVertexLabel);
    consumes<Assns<rec::Track, rec::Vertex> >(fVertexLabel);
    consumes<vector<rec::Vee> >(fVeeLabel);
    consumes<Assns<rec::Track, rec::Vee> >(fVeeLabel);
    consumes<vector<rec::TrackIoniz>>(fTrackLabel);
    consumes<Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);
    consumes<vector<rec::Cluster> >(fClusterLabel);
    //consumes<Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

    fGeo = providerFrom<geo::GeometryGAr>();

}// end constructor()

void BackTrackerAna::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;
    //fTree = tfs->make<TTree>("backtrackTree", "analyze garsoft::BackTracker performance");
    fHistGenDiff = tfs->make<TH1F>("hGenDiff","gen stage;E_{final} - E_{initial} [MeV]",101,-500,500); //genie conserves energy?
    fHistG4Diff = tfs->make<TH1F>("hG4Diff","G4 stage;E_{final} - E_{initial} [GeV]",301,-1000,2000);//geant conserves energy?
    fHistHitComplete = tfs->make<TH1F>("hHitComplete","reco stage: TPC hit;remaining MCParticle energy [GeV]",101,-100,100);
    fHistHitEnergy = tfs->make<TH1F>("hHitEnergy", "TPC Hit -> deposited energy; energy per hit [MeV]",1000,0,10);
    fHistTrackHitDiff = tfs->make<TH1F>("hTrackHitDiff","reco stage: TPC track vs. hit;track - sum hit energy [GeV]",51,-2,2);
    //fHistTrackTotHitSum = tfs->make<TH1F>("hTrackTotHitSum","reco stage: TPC track vs hit; 
    fHistTrackHitMult = tfs->make<TH1F>("hTrackHitMult","reco stage: TPC Track #rightarrow hit;hit duplication",30,-0.5,29.5);
    fHistClustHitMult = tfs->make<TH1F>("hClustHitMult","reco stage: TPC cluster #rightarrow hit;hit duplication",30,-0.5,29.5);
    fHistClustMult = tfs->make<TH1F>("hClustMult","reco stage: TPC track #rightarrow cluster;cluster duplication",30,-0.5,29.5);
    fHistNHitAssocTrackNHitTotDiff = tfs->make<TH1F>("hNHitAssocTrackNHitTotDiff","TPC Track #rightarrow hit;total N^{0} hits - total N^{0} hits assoc'd to tracks",11001,-1000,10000);

    fHistEdepRelKE = tfs->make<TH1F>("hEdepRelKE","Track -> MCParticles;E_{dep}/KE_{initial}",30,0,30);
    fHistG4Diff2d = tfs->make<TH2F>("hG4Diff2d","G4 stage;E_{initial} [GeV]; E_{final} [GeV]",100,0,1000,100,0,1000);
}// end beginJob()

void BackTrackerAna::analyze(art::Event const& e)
{
    cheat::BackTrackerCore const* const_bt = providerFrom<cheat::BackTracker>();
    fBt = const_cast<cheat::BackTrackerCore*>(const_bt);

    // compare neutrino energy to FS particles (not expected to be conserved due to nuclear effects)

    bool warnClusterHitNotFound = true, warnTrackHitNotFound = true;

    //vector< art::Handle< vector<simb::MCTruth> > > mcthandlelist;
    vector< art::Handle< vector<simb::GTruth> > > gthandlelist;
    //e.getManyByType(mcthandlelist);
    e.getManyByType(gthandlelist);

    //const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
    //vector<double> e0, ef;

    //for(auto const& mcth : mcthandlelist) {
    for(auto const& gth : gthandlelist) {

        // loop over single neutrino interactions
        for(auto const& gt : *gth) {
        /*for(auto const& mct : *mcth) {
            const simb::MCNeutrino& mcnu = mct.GetNeutrino();
            const simb::MCParticle& nu = mcnu.Nu();
            double efinal = 0.;

            // loop over final state particles
            for(int ifsp=0; ifsp < mct.NParticles(); ifsp++) {
                const simb::MCParticle& mcp = mct.GetParticle(ifsp);
                cout << "PDG code: " << mcp.PdgCode() << " w/total E = " << mcp.E() << endl;
                if(mcp.PdgCode()!=2212 && mcp.PdgCode()!=2112 && mcp.PdgCode() < 3000) {
                    efinal += mcp.E();
                }
                else if(mcp.PdgCode()==2212 || mcp.PdgCode()==2112){
                    const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
                    efinal += mcp.E() - definition->Mass();
                }
                else if(mcp.PdgCode()>1000000000) {
                    efinal += 
                }
            }*/
            fHistGenDiff->Fill((gt.fFShadSystP4.E()+gt.fFSleptonP4.E()-gt.fProbeP4.E()-gt.fTgtP4.E()));
            //fHistGenDiff->Fill(1000*(efinal - nu.E()));
        }
    }


    // compare FS particle energies from gen stage to initial primary energies in g4 stage

    art::InputTag mcptag(fGeantLabel);
    auto MCPHandle = e.getValidHandle<vector<simb::MCParticle>>(mcptag);
    double eprim = 0., efinal = 0.;

    map<int,double> mcpIdToKE;
    map<int,double> trkIdToERemain;

    //loop over MCParticles
    for ( auto mcp : *MCPHandle ) {

        mcpIdToKE[mcp.TrackId()] = mcp.E()-mcp.Momentum().M(); //initial KE

        if(mcp.PdgCode()>1e9) continue; //ignore nuclei
    
        // is it primary?
        if(mcp.Process()=="primary"){
            eprim += mcp.E();
        }

        if(fGeo->PointInGArTPC(mcp.Position().Vect())) {
            trkIdToERemain[mcp.TrackId()] = mcp.E();
        }

        efinal += mcp.EndE();
        //cout << "MCP (PDG " << mcp.PdgCode() << ") end energy: " << mcp.EndE() << " GeV" << endl;

        if(mcp.E() < mcp.EndE())
            std::cout << "MCParticle with final energy > intial energy ??? Ef-Ei = " << mcp.EndE()-mcp.E() << std::endl;

    }// for MCParticles

    fHistG4Diff->Fill(efinal-eprim);
    fHistG4Diff2d->Fill(eprim,efinal);

    // loop over ALL hits, fill map to count instances in associations to TPCClusters and Tracks
    art::InputTag tpchittag(fHitLabel);
    auto TPCHitHandle = e.getValidHandle< vector<rec::Hit> >(tpchittag);
    size_t nhit = TPCHitHandle->size();
    map<rec::Hit,int> track_hitcounts, cluster_hitcounts;
    for ( rec::Hit const hit : *TPCHitHandle ) {
        //std::cout << hit << std::endl;
        for(cheat::HitIDE const& ide : fBt->HitToHitIDEs(hit)) { //hitIdes) {
            trkIdToERemain[ide.trackID] -= ide.energyTot;
            fHistHitEnergy->Fill(ide.energyTot*1000.);
        }
        //if(track_hitcounts.find(hit.getIDNumber())!=track_hitcounts.end())
        if(track_hitcounts.find(hit)!=track_hitcounts.end()){
            cout << "duplicate hits found in hit list" << endl;
        }
        //track_hitcounts[hit.getIDNumber()] = 0;
        //cluster_hitcounts[hit.getIDNumber()] = 0;
        track_hitcounts[hit] = 0;
        cluster_hitcounts[hit] = 0;
    }//for all hits


    for(auto const& pair : trkIdToERemain) {
        fHistHitComplete->Fill(pair.second);
    }
    //-------------------TPCCluster -> hits ----------------------------------------
    art::InputTag tpcclustag(fTPCClusterLabel);
    auto TPCClusterHandle = e.getValidHandle< vector<rec::TPCCluster> >(tpcclustag);
    map<const rec::IDNumber,int> tpcclustcounts; //map of all TPCClusters to num instances in associations
    int nclusterhit =0, nclusterhit_notfound = 0;
    for(auto clust : *TPCClusterHandle) {
        
        if(tpcclustcounts.find(clust.getIDNumber())!=tpcclustcounts.end())
            cout << "duplicate TPCClusters found!" << endl;
        tpcclustcounts[clust.getIDNumber()] = 0;

        // for hits associated to this cluster
        for(art::Ptr<rec::Hit> const& hit : fBt->TPCClusterToHits(&clust)) {
            nclusterhit++;
            if(cluster_hitcounts.find(*hit)==cluster_hitcounts.end()) {
                nclusterhit_notfound++;
                if(warnClusterHitNotFound) {
                    cout << "TPC hit associated to TPC cluster not found in hit list" << endl; //*********
                    warnClusterHitNotFound = false;
                }
                continue;
            }
            cluster_hitcounts[*hit] += 1;
        }// for matched hits
    }// for all TPCClusters


    cout << "TPCCluster->hits: " << nclusterhit_notfound << " hits out of " << nclusterhit << " not found" << endl;

    //----------- TPC track->hits --------------------------------------
    art::InputTag tracktag(fTrackLabel);
    auto TrackHandle = e.getValidHandle< vector<rec::Track> >(tracktag);
    size_t ntrackhit = 0;
    for(auto track : *TrackHandle ) {

        // ------- TPC track -> clusters -------------------------
        for(auto const& clust : fBt->TrackToClusters(&track) ) {
            if(tpcclustcounts.find(clust->getIDNumber())==tpcclustcounts.end()) {
                //if(warnTrackHitNotFound){
                cout << "TPCCluster associated to track not found in cluster list" << endl;
                    //warnTrackHitNotFound = false;
                //}
                continue;
            }
            tpcclustcounts[clust->getIDNumber()] += 1;
        }

        std::vector<std::pair<simb::MCParticle*,float>> trkmcps = fBt->TrackToMCParticles(&track);
        double etot = fBt->TrackToTotalEnergy(&track);
        double ehits = 0.;
        vector<art::Ptr<rec::Hit>> const hits = fBt->TrackToHits(&track);
        ntrackhit += hits.size();

        for(std::pair<simb::MCParticle*,float> const& mcp : trkmcps) {
            if(mcpIdToKE.find(mcp.first->TrackId())==mcpIdToKE.end()){
                std::cout << "MCParticle matched to track not found in MCParticle list" << std::endl;
                continue;
            }
            double eratio = mcp.second*etot/mcpIdToKE[mcp.first->TrackId()];
            fHistEdepRelKE->Fill(eratio);
            if(eratio>1) std::cout << "bad eratio for trackID " << mcp.first->TrackId() << " with associated energy/initial KE, "
                                   << eratio << std::endl;
        }

        for(auto const& hit : hits) {

            if(track_hitcounts.find(*hit)==track_hitcounts.end()) {
                if(warnTrackHitNotFound){
                    cout << "TPC hit associated to track not found in hit list" << endl;//*********
                    //std::cout << *hit << std::endl;
                    warnTrackHitNotFound = false;
                }
                continue;
            }
            //track_hitcounts[hit.getIDNumber()] += 1;
            track_hitcounts[*hit] += 1;

            for(auto const& ide : fBt->HitToHitIDEs(hit)) {
                ehits += ide.energyTot;
            }
        }
        fHistTrackHitDiff->Fill(etot-ehits);
    }

    for(auto const& cts : track_hitcounts) {
        fHistTrackHitMult->Fill(cts.second);
    }

    for(auto const& cts : cluster_hitcounts) {
        fHistClustHitMult->Fill(cts.second);
    }
    for(auto const& cts : tpcclustcounts) {
        fHistClustMult->Fill(cts.second);
    }
    fHistNHitAssocTrackNHitTotDiff->Fill((int)nhit - (int)ntrackhit);
    //cout << "total # hits - # hits associated with tracks = " << (int)nhit - (int)ntrackhit << endl;
}// end analyze()

DEFINE_ART_MODULE(BackTrackerAna)
