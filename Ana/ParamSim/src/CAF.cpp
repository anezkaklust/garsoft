#include "CAF.h"

#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <functional>
#include <exception>

CAF::CAF()
: cafFile(nullptr), _intfile(nullptr), _inttree(nullptr), _util(new Utils()), _inputfile(""), _outputFile(""), _correct4origin(0)
{

}

CAF::CAF( std::string infile, std::string filename, int correct4origin, double *originTPC)
: cafFile(nullptr), _intfile(nullptr), _inttree(nullptr), _util(new Utils()), _inputfile(infile), _outputFile(filename), _correct4origin(correct4origin)
{
    _util->SetOrigin(originTPC);
}

CAF::~CAF()
{
    delete _util;
    delete cafMVA;
    delete cafFile;
    delete _inttree;
    delete _intfile;
}

bool CAF::BookTFile()
{
    _intfile = new TFile( _inputfile.c_str() );
    if(nullptr == _intfile){
        std::cout << "Cannot open file " << _inputfile.c_str() << std::endl;
        return false;
    }

    _inttree = (TTree*) _intfile->Get( "anatree/GArAnaTree" );
    if(nullptr == _inttree){
        std::cout << "Cannot find tree anatree/GArAnaTree" << std::endl;
        return false;
    }

    cafFile = new TFile(_outputFile.c_str(), "RECREATE");
    if(nullptr == cafFile){
        return false;
    }
    else
    {
        cafMVA = new TTree( "caf", "caf tree" );
        cafMVA->SetDirectory(0);

        //Event number
        cafMVA->Branch("Event", &_Event);
        //Run number
        cafMVA->Branch("Run", &_Run);
        //Sub Run number
        cafMVA->Branch("SubRun", &_SubRun);

        //MC Truth info (Generator)
        cafMVA->Branch("mode", &mode);
        cafMVA->Branch("q2", &q2);
        cafMVA->Branch("w", &w);
        cafMVA->Branch("y", &y);
        cafMVA->Branch("x", &x);
        cafMVA->Branch("theta", &theta);
        cafMVA->Branch("t", &t);
        cafMVA->Branch("ntype", &ntype);
        cafMVA->Branch("ccnc", &ccnc);
        cafMVA->Branch("gint", &gint);
        cafMVA->Branch("tgtpdg", &tgtpdg);
        cafMVA->Branch("weight", &weight);
        cafMVA->Branch("gt_t", &gt_t);
        cafMVA->Branch("intert", &intert);
        cafMVA->Branch("mcnupx", &mcnupx);
        cafMVA->Branch("mcnupy", &mcnupy);
        cafMVA->Branch("mcnupz", &mcnupz);
        cafMVA->Branch("vertx", &vertx);
        cafMVA->Branch("verty", &verty);
        cafMVA->Branch("vertz", &vertz);

        //List of particles from GENIE from the interaction and their status flag
        cafMVA->Branch("nGPart", &nGPart);
        cafMVA->Branch("GPartName", &GPartName);
        cafMVA->Branch("GPartPdg", &GPartPdg);
        cafMVA->Branch("GPartStatus", &GPartStatus);
        cafMVA->Branch("GPartFirstMom", &GPartFirstMom);
        cafMVA->Branch("GPartLastMom", &GPartLastMom);
        cafMVA->Branch("GPartFirstDaugh", &GPartFirstDaugh);
        cafMVA->Branch("GPartLastDaugh", &GPartLastDaugh);
        cafMVA->Branch("GPartPx", &GPartPx);
        cafMVA->Branch("GPartPy", &GPartPy);
        cafMVA->Branch("GPartPz", &GPartPz);
        cafMVA->Branch("GPartE", &GPartE);
        cafMVA->Branch("GPartMass", &GPartMass);

        //Number of final state particle (primaries)
        cafMVA->Branch("nFSP", &_nFSP);
        //MC Particle info
        cafMVA->Branch("detected", &detected);
        cafMVA->Branch("pdgmother", &pdgmother);
        cafMVA->Branch("MCPTime", &mctime);
        cafMVA->Branch("MCPStartX", &_MCPStartX);
        cafMVA->Branch("MCPStartY", &_MCPStartY);
        cafMVA->Branch("MCPStartZ", &_MCPStartZ);
        cafMVA->Branch("motherid", &motherid);
        cafMVA->Branch("mctrkid", &mctrkid);
        cafMVA->Branch("truepx", &truepx);
        cafMVA->Branch("truepy", &truepy);
        cafMVA->Branch("truepz", &truepz);
        cafMVA->Branch("MCPEndX", &_MCPEndX);
        cafMVA->Branch("MCPEndY", &_MCPEndY);
        cafMVA->Branch("MCPEndZ", &_MCPEndZ);
        cafMVA->Branch("MCProc", &_MCProc);
        cafMVA->Branch("MCEndProc", &_MCEndProc);
        cafMVA->Branch("angle", &_angle);
        cafMVA->Branch("truep", &truep);
        cafMVA->Branch("truepdg", &truepdg);
        //Reco info
        cafMVA->Branch("recopid", &recopid);
        cafMVA->Branch("recopidecal", &recopidecal);
        cafMVA->Branch("trkLen", &trkLen);
        cafMVA->Branch("trkLenPerp", &trkLenPerp);
        cafMVA->Branch("preco", &_preco);
        cafMVA->Branch("anglereco", &anglereco);
        cafMVA->Branch("erecon", &erecon);
        cafMVA->Branch("etime", &etime);
        cafMVA->Branch("prob_arr", &prob_arr);
        //Geometry
        cafMVA->Branch("isFidStart", &isFidStart);
        cafMVA->Branch("isTPCStart", &isTPCStart);
        cafMVA->Branch("isCaloStart", &isCaloStart);
        cafMVA->Branch("isInBetweenStart", &isInBetweenStart);
        cafMVA->Branch("isThroughCaloStart", &isThroughCaloStart);
        cafMVA->Branch("isBarrelStart", &isBarrelStart);
        cafMVA->Branch("isEndcapStart", &isEndcapStart);

        cafMVA->Branch("isFidEnd", &isFidEnd);
        cafMVA->Branch("isTPCEnd", &isTPCEnd);
        cafMVA->Branch("isCaloEnd", &isCaloEnd);
        cafMVA->Branch("isThroughCaloEnd", &isThroughCaloEnd);
        cafMVA->Branch("isInBetweenEnd", &isInBetweenEnd);
        cafMVA->Branch("isBarrelEnd", &isBarrelEnd);
        cafMVA->Branch("isEndcapEnd", &isEndcapEnd);

        return true;
    }
}

void CAF::CloseTFile()
{
    if(cafFile != nullptr) {
        cafFile->Close();
        return;
    }
    else{ return; }
}

void CAF::WriteTTree()
{
    if(cafMVA != nullptr) {
        cafFile->cd();
        cafMVA->Write();
        cafFile->Close();
    }

    return;
}

void CAF::FillTTree()
{
    if(cafMVA != nullptr) {
        cafFile->cd();
        cafMVA->Fill();
    }

    return;
}

void CAF::ClearVectors()
{
    //Generator values
    mode.clear();
    ccnc.clear();
    ntype.clear();
    gint.clear();
    weight.clear();
    tgtpdg.clear();
    gt_t.clear();
    intert.clear();
    q2.clear();
    w.clear();
    y.clear();
    x.clear();
    theta.clear();
    t.clear();
    mcnupx.clear();
    mcnupy.clear();
    mcnupz.clear();
    vertx.clear();
    verty.clear();
    vertz.clear();

    nGPart.clear();
    GPartPdg.clear();
    GPartStatus.clear();
    GPartName.clear();
    GPartFirstMom.clear();
    GPartLastMom.clear();
    GPartFirstDaugh.clear();
    GPartLastDaugh.clear();
    GPartPx.clear();
    GPartPy.clear();
    GPartPz.clear();
    GPartE.clear();
    GPartMass.clear();

    //MC Particle values
    _nFSP.clear();
    pdgmother.clear();
    truepdg.clear();
    mctime.clear();
    mctrkid.clear();
    motherid.clear();
    _MCPStartX.clear();
    _MCPStartY.clear();
    _MCPStartZ.clear();
    _MCPEndX.clear();
    _MCPEndY.clear();
    _MCPEndZ.clear();
    _MCProc.clear();
    _MCEndProc.clear();
    trkLen.clear();
    trkLenPerp.clear();
    truep.clear();
    truepx.clear();
    truepy.clear();
    truepz.clear();
    _angle.clear();

    //Reco values
    detected.clear();
    recopid.clear();
    recopidecal.clear();
    prob_arr.clear();
    anglereco.clear();
    _preco.clear();
    erecon.clear();
    etime.clear();

    //Geometry
    isFidStart.clear();
    isTPCStart.clear();
    isCaloStart.clear();
    isThroughCaloStart.clear();
    isInBetweenStart.clear();
    isBarrelStart.clear();
    isEndcapStart.clear();

    isFidEnd.clear();
    isTPCEnd.clear();
    isCaloEnd.clear();
    isThroughCaloEnd.clear();
    isInBetweenEnd.clear();
    isBarrelEnd.clear();
    isEndcapEnd.clear();
}

bool CAF::CheckVectorSize()
{
    bool isOK = true;

    if(_nFSP.at(0) != pdgmother.size()) isOK = false;
    if(_nFSP.at(0) != truepdg.size()) isOK = false;
    if(_nFSP.at(0) != mctime.size()) isOK = false;
    if(_nFSP.at(0) != mctrkid.size()) isOK = false;
    if(_nFSP.at(0) != motherid.size()) isOK = false;
    if(_nFSP.at(0) != _MCPStartX.size()) isOK = false;
    if(_nFSP.at(0) != _MCPStartY.size()) isOK = false;
    if(_nFSP.at(0) != _MCPStartZ.size()) isOK = false;
    if(_nFSP.at(0) != _MCPEndX.size()) isOK = false;
    if(_nFSP.at(0) != _MCPEndY.size()) isOK = false;
    if(_nFSP.at(0) != _MCPEndZ.size()) isOK = false;
    if(_nFSP.at(0) != _MCProc.size()) isOK = false;
    if(_nFSP.at(0) != _MCEndProc.size()) isOK = false;
    if(_nFSP.at(0) != trkLen.size()) isOK = false;
    if(_nFSP.at(0) != trkLenPerp.size()) isOK = false;
    if(_nFSP.at(0) != truep.size()) isOK = false;
    if(_nFSP.at(0) != truepx.size()) isOK = false;
    if(_nFSP.at(0) != truepy.size()) isOK = false;
    if(_nFSP.at(0) != truepz.size()) isOK = false;
    if(_nFSP.at(0) != _angle.size()) isOK = false;

    //Reco values
    if(_nFSP.at(0) != recopid.size()) isOK = false;
    if(_nFSP.at(0) != detected.size()) isOK = false;
    if(_nFSP.at(0) != recopidecal.size()) isOK = false;
    // if(_nFSP.at(0) != prob_arr.size()) isOK = false;
    if(_nFSP.at(0) != anglereco.size()) isOK = false;
    if(_nFSP.at(0) != _preco.size()) isOK = false;
    if(_nFSP.at(0) != erecon.size()) isOK = false;
    if(_nFSP.at(0) != etime.size()) isOK = false;

    //Geometry
    if(_nFSP.at(0) != isFidStart.size()) isOK = false;
    if(_nFSP.at(0) != isTPCStart.size()) isOK = false;
    if(_nFSP.at(0) != isCaloStart.size()) isOK = false;
    if(_nFSP.at(0) != isThroughCaloStart.size()) isOK = false;
    if(_nFSP.at(0) != isInBetweenStart.size()) isOK = false;
    if(_nFSP.at(0) != isBarrelStart.size()) isOK = false;
    if(_nFSP.at(0) != isEndcapStart.size()) isOK = false;

    if(_nFSP.at(0) != isFidEnd.size()) isOK = false;
    if(_nFSP.at(0) != isTPCEnd.size()) isOK = false;
    if(_nFSP.at(0) != isCaloEnd.size()) isOK = false;
    if(_nFSP.at(0) != isThroughCaloEnd.size()) isOK = false;
    if(_nFSP.at(0) != isInBetweenEnd.size()) isOK = false;
    if(_nFSP.at(0) != isBarrelEnd.size()) isOK = false;
    if(_nFSP.at(0) != isEndcapEnd.size()) isOK = false;

    return isOK;
}

float CAF::calcGluck(double sigmaX, double B, double X0, float nHits, double mom, double length, double& ratio)
{
  double sig_meas = sqrt(720./(nHits+4))*( (0.01*sigmaX*mom)/(0.0001*0.3*B*length*length) );
  double sig_mcs = (0.052/B)*sqrt( 1.43/(0.0001*X0*length) );
  double sigma = sqrt(sig_meas*sig_meas + sig_mcs*sig_mcs);
  ratio = sig_meas/sig_mcs;
 
  return sigma;
}

// main loop function
void CAF::loop()
{
    //double gastpc_len = 5.; // track length cut in cm
    const float gastpc_len = 2.; // new track length cut in cm based on Thomas' study of low energy protons
    // dont care about electrons -- check momentum and see if hit ECAL
    const float gastpc_B = 0.5; // B field strength in Tesla
    const float gastpc_padPitch = 1.0; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
    const float gastpc_X0 = 1300.; // cm = 13m radiation length
    //Resolution for short tracks //TODO check this numbers!
    const float sigmaP_short = 0.1; //in GeV
    // point resolution
    const float sigma_x = 0.1;

    std::vector<float> v;
    for (float pit = 0.040; pit < 20.0; pit += 0.001)
    {
        v.push_back(pit);
    }

    //as function of KE
    //0 -> 50 MeV ~ 20%
    //> 50 MeV ~ 40%
    const float NeutronECAL_detEff[2] = {0.2, 0.4};
    const float sigmaNeutronECAL_first = 0.11;
    //TODO fraction of rescatters
    // float sigmaNeutronECAL_rescatter = 0.26;
    //ECAL energy resolution sigmaE/E
    const float ECAL_stock = 0.06; //in %
    const float ECAL_const = 0.02;
    TF1 *fRes = new TF1("fRes", "TMath::Sqrt ( [0]*[0]/x + [1]*[1] )");
    fRes->FixParameter(0, ECAL_stock);
    fRes->FixParameter(1, ECAL_const);
    //ECAL sampling fraction
    // double sampling_frac = 4.32;
    //ECAL nlayers
    const int nLayers = 60;
    //ECAL MIP resolution (based on AHCAL numbers)
    const double ECAL_MIP_Res = 0.23;
    //MIP2GeV conversion factor
    const double MIP2GeV_factor = 0.814 / 1000;
    //float ECAL_pi0_resolution = 0.13; //sigmaE/E in between at rest (17%) and high energy (~few %)
    const float ECAL_time_resolution = 1.; // 1 ns time resolution
    TParticlePDG *neutron = TDatabasePDG::Instance()->GetParticle(2112);
    const float neutron_mass = neutron->Mass(); //in GeV
    //TParticlePDG *pi0 = TDatabasePDG::Instance()->GetParticle(111);
    //float pi0_mass = pi0->Mass(); //in GeV

    TString filename="${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
    TFile infile(filename,"READ");

    m_pidinterp.clear();
    char str[11];
    for (int q = 0; q < 501; ++q)
    {
        sprintf(str, "%d", q);
        std::string s = "pidmatrix";
        s.append(str);
        // read the 500 histograms one by one; each histogram is a
        // 6 by 6 matrix of probabilities for a given momentum value
        m_pidinterp.insert( std::make_pair(q, (TH2F*) infile.Get(s.c_str())->Clone("pidinterp")) );
    }

    // infile.Close();

    //------------------------------------------------------------------------

    int           	 Event = 0;
    int           	 SubRun = 0;
    int           	 Run = 0;

    std::vector<float> *MC_Q2 = 0;
    TBranch            *b_MC_Q2 = 0;
    std::vector<float> *MC_W = 0;
    TBranch            *b_MC_W = 0;
    std::vector<float> *MC_Y = 0;
    TBranch            *b_MC_Y = 0;
    std::vector<float> *MC_X = 0;
    TBranch            *b_MC_X = 0;
    std::vector<float> *MC_Theta = 0;
    TBranch            *b_MC_Theta = 0;
    std::vector<float> *MCVertX = 0;
    TBranch            *b_MCVertX = 0;
    std::vector<float> *MCVertY = 0;
    TBranch            *b_MCVertY = 0;
    std::vector<float> *MCVertZ = 0;
    TBranch            *b_MCVertZ = 0;
    std::vector<float> *MCNuPx = 0;
    TBranch            *b_MCNuPx = 0;
    std::vector<float> *MCNuPy = 0;
    TBranch            *b_MCNuPy = 0;
    std::vector<float> *MCNuPz = 0;
    TBranch            *b_MCNuPz = 0;
    std::vector<int>   *NType = 0;
    TBranch            *b_NType = 0;
    std::vector<int>   *CCNC = 0;
    TBranch            *b_CCNC = 0;
    std::vector<int>   *Mode = 0;
    TBranch            *b_Mode = 0;
    std::vector<int>   *Gint=0;
    TBranch            *b_Gint = 0;
    std::vector<int>   *TgtPDG = 0;
    TBranch            *b_TgtPDG = 0;
    std::vector<int>   *GT_T = 0;
    TBranch            *b_GT_T = 0;
    std::vector<int>   *InterT=0;
    TBranch            *b_InterT = 0;
    std::vector<float> *Weight = 0;
    TBranch            *b_Weight = 0;

    std::vector<int> *_nGPart = 0;
    TBranch            *b_nGPart = 0;
    std::vector<int> *_GPartPdg = 0;
    TBranch            *b_GPartPdg = 0;
    std::vector<int> *_GPartStatus = 0;
    TBranch            *b_GPartStatus = 0;
    std::vector<std::string> *_GPartName = 0;
    TBranch            *b_GPartName = 0;
    std::vector<int> *_GPartFirstMom = 0;
    TBranch            *b_GPartFirstMom = 0;
    std::vector<int> *_GPartLastMom = 0;
    TBranch            *b_GPartLastMom = 0;
    std::vector<int> *_GPartFirstDaugh = 0;
    TBranch            *b_GPartFirstDaugh = 0;
    std::vector<int> *_GPartLastDaugh = 0;
    TBranch            *b_GPartLastDaugh = 0;
    std::vector<float> *_GPartPx = 0;
    TBranch            *b_GPartPx = 0;
    std::vector<float> *_GPartPy = 0;
    TBranch            *b_GPartPy = 0;
    std::vector<float> *_GPartPz = 0;
    TBranch            *b_GPartPz = 0;
    std::vector<float> *_GPartE = 0;
    TBranch            *b_GPartE = 0;
    std::vector<float> *_GPartMass = 0;
    TBranch            *b_GPartMass = 0;

    std::vector<int>     *PDG = 0;
    TBranch              *b_PDG = 0;
    std::vector<int>     *MCPTrkID = 0;
    TBranch              *b_MCPTrkID = 0;
    std::vector<int>     *MCMotherTrkID = 0;
    TBranch              *b_MCMotherTrkID = 0;
    std::vector<int>     *PDGMother = 0;
    TBranch              *b_PDGMother = 0;
    std::vector<float>   *MCPTime = 0;
    TBranch              *b_MCPTime = 0;
    std::vector<float>   *MCPStartX = 0;
    TBranch              *b_MCPStartX = 0;
    std::vector<float>   *MCPStartY = 0;
    TBranch              *b_MCPStartY = 0;
    std::vector<float>   *MCPStartZ = 0;
    TBranch              *b_MCPStartZ = 0;
    std::vector<float>   *MCPStartPX = 0;
    TBranch              *b_MCPStartPX = 0;
    std::vector<float>   *MCPStartPY = 0;
    TBranch              *b_MCPStartPY = 0;
    std::vector<float>   *MCPStartPZ = 0;
    TBranch              *b_MCPStartPZ = 0;
    std::vector<std::string>   *MCPProc = 0;
    TBranch                    *b_MCPProc = 0;
    std::vector<std::string>   *MCPEndProc = 0;
    TBranch                    *b_MCPEndProc = 0;
    std::vector<float>   *MCPEndX = 0;
    TBranch              *b_MCPEndX = 0;
    std::vector<float>   *MCPEndY = 0;
    TBranch              *b_MCPEndY = 0;
    std::vector<float>   *MCPEndZ = 0;
    TBranch              *b_MCPEndZ = 0;
    std::vector<float>   *TrajMCPX = 0;
    TBranch              *b_TrajMCPX = 0;
    std::vector<float>   *TrajMCPY = 0;
    TBranch              *b_TrajMCPY = 0;
    std::vector<float>   *TrajMCPZ = 0;
    TBranch              *b_TrajMCPZ = 0;
    std::vector<int>     *TrajMCPTrajIndex = 0;
    TBranch              *b_TrajMCPTrajIndex = 0;

    _inttree->SetBranchAddress("Event", &Event);
    _inttree->SetBranchAddress("SubRun", &SubRun);
    _inttree->SetBranchAddress("Run", &Run);

    //Generator info
    _inttree->SetBranchAddress("NType", &NType, &b_NType);
    _inttree->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
    _inttree->SetBranchAddress("MC_Q2", &MC_Q2, &b_MC_Q2);
    _inttree->SetBranchAddress("MC_W", &MC_W, &b_MC_W);
    _inttree->SetBranchAddress("MC_Y", &MC_Y, &b_MC_Y);
    _inttree->SetBranchAddress("MC_X", &MC_X, &b_MC_X);
    _inttree->SetBranchAddress("MC_Theta", &MC_Theta, &b_MC_Theta);
    _inttree->SetBranchAddress("Mode", &Mode, &b_Mode);
    _inttree->SetBranchAddress("Gint", &Gint, &b_Gint);
    _inttree->SetBranchAddress("TgtPDG", &TgtPDG, &b_TgtPDG);
    _inttree->SetBranchAddress("GT_T", &GT_T, &b_GT_T);
    _inttree->SetBranchAddress("MCVertX", &MCVertX, &b_MCVertX);
    _inttree->SetBranchAddress("MCVertY", &MCVertY, &b_MCVertY);
    _inttree->SetBranchAddress("MCVertZ", &MCVertZ, &b_MCVertZ);
    _inttree->SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
    _inttree->SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
    _inttree->SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
    _inttree->SetBranchAddress("InterT", &InterT, &b_InterT);
    _inttree->SetBranchAddress("Weight", &Weight, &b_Weight);

    _inttree->SetBranchAddress("nGPart", &_nGPart, &b_nGPart);
    _inttree->SetBranchAddress("GPartPdg", &_GPartPdg, &b_GPartPdg);
    _inttree->SetBranchAddress("GPartStatus", &_GPartStatus, &b_GPartStatus);
    _inttree->SetBranchAddress("GPartName", &_GPartName, &b_GPartName);
    _inttree->SetBranchAddress("GPartFirstMom", &_GPartFirstMom, &b_GPartFirstMom);
    _inttree->SetBranchAddress("GPartLastMom", &_GPartLastMom, &b_GPartLastMom);
    _inttree->SetBranchAddress("GPartFirstDaugh", &_GPartFirstDaugh, &b_GPartFirstDaugh);
    _inttree->SetBranchAddress("GPartLastDaugh", &_GPartLastDaugh, &b_GPartLastDaugh);
    _inttree->SetBranchAddress("GPartPx", &_GPartPx, &b_GPartPx);
    _inttree->SetBranchAddress("GPartPy", &_GPartPy, &b_GPartPy);
    _inttree->SetBranchAddress("GPartPz", &_GPartPz, &b_GPartPz);
    _inttree->SetBranchAddress("GPartE", &_GPartE, &b_GPartE);
    _inttree->SetBranchAddress("GPartMass", &_GPartMass, &b_GPartMass);

    //MC info
    _inttree->SetBranchAddress("PDG", &PDG, &b_PDG);
    _inttree->SetBranchAddress("MCTrkID", &MCPTrkID, &b_MCPTrkID);
    _inttree->SetBranchAddress("MotherTrkID", &MCMotherTrkID, &b_MCMotherTrkID);
    _inttree->SetBranchAddress("MCPTime", &MCPTime, &b_MCPTime);
    _inttree->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    _inttree->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    _inttree->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
    _inttree->SetBranchAddress("MCPEndX", &MCPEndX, &b_MCPEndX);
    _inttree->SetBranchAddress("MCPEndY", &MCPEndY, &b_MCPEndY);
    _inttree->SetBranchAddress("MCPEndZ", &MCPEndZ, &b_MCPEndZ);
    _inttree->SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    _inttree->SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    _inttree->SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    _inttree->SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    _inttree->SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);
    _inttree->SetBranchAddress("MCPEndProc", &MCPEndProc, &b_MCPEndProc);
    _inttree->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    _inttree->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    _inttree->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    // _inttree->SetBranchAddress("TrajMCPTrajIndex", &TrajMCPTrajIndex);
    _inttree->SetBranchAddress("TrajMCPTrackID", &TrajMCPTrajIndex, &b_TrajMCPTrajIndex);

    //gamma, neutron, pi0, k0L, k0S, k0, delta0
    std::vector<int> neutrinos = {12, 14, 16};
    std::vector<int> pdg_neutral = {22, 2112, 111, 130, 310, 311, 2114};
    //pion, muon, proton, kaon, deuteron, electron
    std::vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
    //-------------------------------------------------------------------

    // Main event loop
    for( int entry = 0; entry < _inttree->GetEntries(); entry++ )
    {
        _inttree->GetEntry(entry);
        this->ClearVectors();

        //Filling MCTruth values
        // std::cout << "Event " << Event << " Run " << Run << std::endl;
        // if(Event < 259) continue;

        _Event = Event;
        _Run = Run;
        _SubRun = SubRun;

        for(size_t i = 0; i < NType->size(); i++)
        {
            ccnc.push_back(CCNC->at(i));
            ntype.push_back(NType->at(i));
            q2.push_back(MC_Q2->at(i));
            w.push_back(MC_W->at(i));
            y.push_back(MC_Y->at(i));
            x.push_back(MC_X->at(i));
            theta.push_back(MC_Theta->at(i));
            mode.push_back(Mode->at(i));
            intert.push_back(InterT->at(i));
            if(_correct4origin){
                vertx.push_back(MCVertX->at(i) - _util->GetOrigin()[0]);
                verty.push_back(MCVertY->at(i) - _util->GetOrigin()[1]);
                vertz.push_back(MCVertZ->at(i) - _util->GetOrigin()[2]);
            } else {
                vertx.push_back(MCVertX->at(i));
                verty.push_back(MCVertY->at(i));
                vertz.push_back(MCVertZ->at(i));
            }
            mcnupx.push_back(MCNuPx->at(i));
            mcnupy.push_back(MCNuPy->at(i));
            mcnupz.push_back(MCNuPz->at(i));
        }

        for(size_t i = 0; i < Gint->size(); i++)
        {
            gint.push_back(Gint->at(i));
            tgtpdg.push_back(TgtPDG->at(i));//target pdg
            gt_t.push_back(GT_T->at(i));
            weight.push_back(Weight->at(i));
        }

        for(size_t i = 0; i < _nGPart->size(); i++)
        {
            nGPart.push_back(_nGPart->at(i));
            for(int j = 0; j < _nGPart->at(i); j++) {
                GPartPdg.push_back(_GPartPdg->at(j));
                GPartStatus.push_back(_GPartStatus->at(j));
                GPartName.push_back(_GPartName->at(j));
                GPartFirstMom.push_back(_GPartFirstMom->at(j));
                GPartLastMom.push_back(_GPartLastMom->at(j));
                GPartFirstDaugh.push_back(_GPartFirstDaugh->at(j));
                GPartLastDaugh.push_back(_GPartLastDaugh->at(j));
                GPartPx.push_back(_GPartPx->at(j));
                GPartPy.push_back(_GPartPy->at(j));
                GPartPz.push_back(_GPartPz->at(j));
                GPartE.push_back(_GPartE->at(j));
                GPartMass.push_back(_GPartMass->at(j));
            }
        }

        //--------------------------------------------------------------------------
        // Start of Parameterized Reconstruction
        //--------------------------------------------------------------------------
        unsigned int nFSP = 0;

        //TODO we should skip particles that we don't see or cannot reconstruct! change filling caf with i to index counting the particle
        //have to be careful with indexes and continue
        //---------------------------------------------------------------
        // all Gluckstern calculations happen in the following loop
        for(size_t i = 0; i < MCPStartPX->size(); i++ )
        {
            //Get the creating process
            std::string mcp_process = MCPProc->at(i);
            //Get ending process
            std::string mcp_endprocess = MCPEndProc->at(i);
            int mctrackid = MCPTrkID->at(i);
            int pdg = PDG->at(i);
            const TVector3 spoint(MCPStartX->at(i) - _util->GetOrigin()[0], MCPStartY->at(i) - _util->GetOrigin()[1], MCPStartZ->at(i) - _util->GetOrigin()[2]);
            const TVector3 epoint(MCPEndX->at(i) - _util->GetOrigin()[0], MCPEndY->at(i) - _util->GetOrigin()[1], MCPEndZ->at(i) - _util->GetOrigin()[2]);

            TVector3 mcp(MCPStartPX->at(i), MCPStartPY->at(i), MCPStartPZ->at(i));
            float ptrue = (mcp).Mag();
            float pypz = std::sqrt( MCPStartPY->at(i)*MCPStartPY->at(i) + MCPStartPZ->at(i)*MCPStartPZ->at(i));

            //need to ignore neutrals for this - put the value to 0
            auto result = std::find(pdg_neutral.begin(), pdg_neutral.end(), abs(pdg));
            bool isNeutral = (result != pdg_neutral.end()) ? true : false;

            //start track length
            //***************************************************************************************************************/

            if( isNeutral )
            {
                trkLen.push_back(-1);
                trkLenPerp.push_back(-1);
            }
            else
            {
                // calculate the total and the transverse track lengths and restrict the
                // tracklength to be above the gas TPC track length threshold
                double tracklen = 0.;
                double tracklen_perp = 0.;

                //CAREFUL No offset for the trajectory points (origin for them is the TPC?)??????
                //TODO check if the mcp point is within the TPC volume! Skip for mcp in the ECAL (showers)
                //TODO Link showers to original mcp?
                for(size_t itraj = 1; itraj < TrajMCPX->size(); itraj++)
                {
                    //check that it is the correct mcp
                    if(TrajMCPTrajIndex->at(itraj) == mctrackid)
                    {
                        //Traj point+1
                        TVector3 point(TrajMCPX->at(itraj) - _util->GetOrigin()[0], TrajMCPY->at(itraj) - _util->GetOrigin()[1], TrajMCPZ->at(itraj) - _util->GetOrigin()[2]);

                        //point is not in the TPC anymore - stop traj loop
                        if(not _util->PointInTPC(point))
                        {
                            // // std::cout << "Point not within the TPC: " << point.X() << " r " << std::sqrt(point.Y()*point.Y() + point.Z()*point.Z()) << std::endl;
                            continue;
                        }

                        // find the length of the track by getting the distance between each hit
                        TVector3 diff(TrajMCPX->at(itraj) - TrajMCPX->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1), TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1));
                        // perp length
                        TVector2 tracklen_perp_vec(TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1));
                        // Summing up
                        tracklen += diff.Mag();
                        tracklen_perp += tracklen_perp_vec.Mod();
                    }
                }

                trkLen.push_back(tracklen);
                trkLenPerp.push_back(tracklen_perp);
            }

            //end track length
            //***************************************************************************************************************/

            TVector3 xhat(1, 0, 0);
            // float pz = mcp.Z();
            // float pt = (mcp.Cross(xhat)).Mag();
            // float px = mcp.X();
            // float py = mcp.Y();
            //float mctrackid = MCPTrkID->at(i);
            // angle with respect to the incoming neutrino
            float angle  = atan(mcp.X() / mcp.Z());
            float ecaltime = _util->GaussianSmearing(MCPTime->at(i), ECAL_time_resolution);
            float time = MCPTime->at(i);

            //Check where start point is for the mcp
            isFidStart.push_back(_util->PointInFiducial(spoint));
            isTPCStart.push_back(_util->PointInTPC(spoint));
            isCaloStart.push_back(_util->PointInCalo(spoint));
            isThroughCaloStart.push_back(_util->isThroughCalo(spoint));
            isInBetweenStart.push_back(_util->PointStopBetween(spoint));
            isBarrelStart.push_back(_util->isBarrel(spoint));
            isEndcapStart.push_back(_util->isEndcap(spoint));
            //Check where endpoint of mcp is
            isFidEnd.push_back(_util->PointInFiducial(epoint));
            isTPCEnd.push_back(_util->PointInTPC(epoint));
            isCaloEnd.push_back(_util->PointInCalo(epoint));
            isThroughCaloEnd.push_back(_util->isThroughCalo(epoint));
            isInBetweenEnd.push_back(_util->PointStopBetween(epoint));
            isBarrelEnd.push_back(_util->isBarrel(epoint));
            isEndcapEnd.push_back(_util->isEndcap(epoint));

            //start tpc
            //***************************************************************************************************************/

            //Visible in the TPC
            if( trkLen.at(nFSP) > gastpc_len )
            {
                for (int pidm = 0; pidm < 6; ++pidm)
                {
                    if ( abs(pdg) == pdg_charged.at(pidm) )
                    {
                        // // std::cout << "Entered reco TPC" << std::endl;

                        //Use range instead of Gluckstern for stopping tracks
                        //TODO is that correct? What if it is a scatter in the TPC? Need to check if daughter is same particle
                        float preco = 0;

                        // save the true PDG, parametrized PID comes later
                        truepdg.push_back(pdg);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        if(_correct4origin){
                            _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                            _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                            _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                            _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                            _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                            _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                        } else {
                            _MCPStartX.push_back(MCPStartX->at(i));
                            _MCPStartY.push_back(MCPStartY->at(i));
                            _MCPStartZ.push_back(MCPStartZ->at(i));
                            _MCPEndX.push_back(MCPEndX->at(i));
                            _MCPEndY.push_back(MCPEndY->at(i));
                            _MCPEndZ.push_back(MCPEndZ->at(i));
                        }
                        pdgmother.push_back(PDGMother->at(i));
                        // save the true momentum
                        truep.push_back(ptrue);
                        // save the true angle
                        _angle.push_back(angle);
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                        mctime.push_back(time);
                        mctrkid.push_back(MCPTrkID->at(i));
                        motherid.push_back(MCMotherTrkID->at(i));

                        //Case for range, the end point of the mcp is in the TPC, does not reach the ecal
                        if( _util->PointInTPC(epoint) )
                        {
                            // calculate number of trackpoints
                            float nHits = round (trkLen.at(nFSP) / gastpc_padPitch);
                            // angular resolution first term
                            float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / trkLen.at(nFSP)*trkLen.at(nFSP)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
                            // scattering term in Gluckstern formula
                            float sigma_angle_2 = (0.015*0.015 / (3. * ptrue * ptrue)) * (trkLen.at(nFSP)/gastpc_X0);
                            // angular resolution from the two terms above
                            float sigma_angle_short = sqrt(sigma_angle_1 + sigma_angle_2);
                            //reconstructed angle
                            float angle_reco = _util->GaussianSmearing(angle, sigma_angle_short);
                            //reconstructed momentum
                            preco = _util->GaussianSmearing( ptrue, sigmaP_short );

                            if(preco > 0)
                            _preco.push_back(preco);
                            else
                            _preco.push_back(-1);

                            anglereco.push_back(angle_reco);

                            erecon.push_back(-1);
                            recopidecal.push_back(-1);
                            etime.push_back(-1);
                            detected.push_back(-1);
                        }
                        else
                        {
                            //Case where the endpoint is not in the TPC, should be able to use the Gluckstern formula
                            // calculate number of trackpoints
                            // float nHits = round (trkLen.at(nFSP) / gastpc_padPitch);
                            // // measurement term in Gluckstern formula
                            // float fracSig_meas = sqrt(720./(nHits+4)) * ((0.01*gastpc_padPitch*ptrue) / (0.3 * gastpc_B * 0.0001 *trkLenPerp.at(nFSP)*trkLenPerp.at(nFSP)));
                            // // multiple Coulomb scattering term in Gluckstern formula
                            // float fracSig_MCS = (0.052*sqrt(1.43)) / (gastpc_B * sqrt(gastpc_X0*trkLenPerp.at(nFSP)*0.0001));
                            // // momentum resoltion from the two terms above
                            // float sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
                            // // now Gaussian smear the true momentum using the momentum resolution
                            // preco = _util->GaussianSmearing( ptrue, sigmaP );

                            //Vivek Gluckstern
                            int nHits = round (trkLen.at(nFSP) / gastpc_padPitch);
                            double ratio = 0.;
                            float sigmaPt  = pypz*calcGluck(sigma_x, gastpc_B, gastpc_X0, nHits, pypz, trkLenPerp.at(nFSP), ratio);
                            float tan_lambda = MCPStartPX->at(i)/pypz;
                            float lambda = atan2(MCPStartPX->at(i), pypz); // between -pi and +pi
                            float del_lambda = 0.0062;
                            float temp = pypz*tan_lambda*del_lambda;
                            float sigmaP = (1./cos(lambda))*sqrt(sigmaPt*sigmaPt + temp*temp);
                            // now Gaussian smear the true momentum using the momentum resolution
                            preco = _util->GaussianSmearing( ptrue, sigmaP );

                            // measurement term in the Gluckstern formula for calculating the
                            // angular resolution
                            float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / trkLen.at(nFSP)*trkLen.at(nFSP)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
                            // scattering term in Gluckstern formula
                            float sigma_angle_2 = (0.015*0.015 / (3. * ptrue * ptrue)) * (trkLen.at(nFSP)/gastpc_X0);
                            // angular resolution from the two terms above
                            float sigma_angle = sqrt(sigma_angle_1 + sigma_angle_2);
                            // now Gaussian smear the true angle using the angular resolution
                            float angle_reco = _util->GaussianSmearing(angle, sigma_angle);

                            // save reconstructed momentum and angle to cafanatree
                            if(preco > 0)
                            _preco.push_back(preco);
                            else
                            _preco.push_back(-1);

                            anglereco.push_back(angle_reco);

                            //Reaches the ECAL and stops there
                            if( _util->PointInCalo(epoint) )
                            {
                                //Need energy measurement in ecal
                                TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(pdg));
                                if(nullptr == part)
                                {
                                    // std::cout << "Could not find particle in root pdg table, pdg " << pdg << std::endl;
                                    //deuteron
                                    if( pdg == 1000010020 ) {
                                        float mass = 1.8756;//in GeV mass deuteron
                                        float etrue = std::sqrt(ptrue*ptrue + mass*mass) - mass;//needs to be KE here
                                        float ECAL_resolution = fRes->Eval(etrue)*etrue;
                                        float ereco = _util->GaussianSmearing(etrue, ECAL_resolution);
                                        erecon.push_back(ereco);
                                        recopidecal.push_back(-1);
                                        detected.push_back(1);
                                        etime.push_back(ecaltime);
                                    }
                                    else
                                    {
                                        erecon.push_back(-1);
                                        recopidecal.push_back(-1);
                                        detected.push_back(0);
                                        etime.push_back(-1);
                                    }
                                }
                                else
                                {
                                    //by default should be tagged as an electron as it has a track,
                                    //otherwise tag as gamma if not track -> need mis association rate, and use dE/dX in Scintillator?
                                    //separation between e and mu/pi should be around 100%
                                    //separation mu/pi -> based on Chris study with only the ECAL (no Muon ID detector)
                                    //separation with p and mu/pi/e ?? high energy -> confusion with mu/pi, low energy confusion with e
                                    //using E/p to ID?
                                    float mass = part->Mass();//in GeV
                                    float etrue = std::sqrt(ptrue*ptrue + mass*mass);//Just energy here is fine
                                    float ECAL_resolution = fRes->Eval(etrue)*etrue;
                                    float ereco = _util->GaussianSmearing(etrue, ECAL_resolution);
                                    erecon.push_back((ereco > 0) ? ereco : 0.);
                                    detected.push_back(1);
                                    etime.push_back(ecaltime);
                                    // // std::cout << "E/p " << ereco/preco << " true pdg " << pdg << std::endl;

                                    //Electron
                                    if( abs(pdg) == 11 ){
                                        recopidecal.push_back(11);
                                    }
                                    else if( abs(pdg) == 13 || abs(pdg) == 211 )
                                    {
                                        //Muons and Pions
                                        //ptrue < 480 MeV/c 100% separation
                                        //80% from 480 to 750
                                        //90% up to 750 to 900
                                        //95% over 900
                                        float random_number = _util->GetRamdomNumber();

                                        if(ptrue < 0.48) {
                                            recopidecal.push_back(abs(pdg));//100% efficiency by range
                                        }
                                        else if(ptrue >= 0.48 && ptrue < 0.75)
                                        {
                                            //case muon
                                            if(abs(pdg) == 13)
                                            {
                                                if(random_number > (1 - 0.8)) {
                                                    recopidecal.push_back(13);
                                                }
                                                else{
                                                    recopidecal.push_back(211);
                                                }
                                            }

                                            //case pion
                                            if(abs(pdg) == 211)
                                            {
                                                if(random_number > (1 - 0.8)) {
                                                    recopidecal.push_back(211);
                                                }
                                                else{
                                                    recopidecal.push_back(13);
                                                }
                                            }
                                        }
                                        else if(ptrue >= 0.75 && ptrue < 0.9)
                                        {
                                            //case muon
                                            if(abs(pdg) == 13){
                                                if(random_number > (1 - 0.9)) {
                                                    recopidecal.push_back(13);
                                                }
                                                else{
                                                    recopidecal.push_back(211);
                                                }
                                            }
                                            //case pion
                                            if(abs(pdg) == 211) {
                                                if(random_number > (1 - 0.9)) {
                                                    recopidecal.push_back(211);
                                                }
                                                else{
                                                    recopidecal.push_back(13);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            //case muon
                                            if(abs(pdg) == 13){
                                                if(random_number > (1 - 0.95)) {
                                                    recopidecal.push_back(13);
                                                }
                                                else{
                                                    recopidecal.push_back(211);
                                                }
                                            }
                                            //case pion
                                            if(abs(pdg) == 211){
                                                if(random_number > (1 - 0.95)) {
                                                    recopidecal.push_back(211);
                                                }
                                                else{
                                                    recopidecal.push_back(13);
                                                }
                                            }
                                        }
                                    }
                                    else if( abs(pdg) == 2212 )
                                    {
                                        recopidecal.push_back(2212);//TODO for p/pi separation
                                    }
                                    else {
                                        recopidecal.push_back(-1);
                                    }
                                }
                            }
                            else if( _util->isThroughCalo(epoint) )
                            {
                                //Case the endpoint is outside the CALO -> it went through the ECAL (mu/pi/p possible)
                                //the ECAL will see 60 MIPs on average
                                double Evis = (double)nLayers; //in MIP
                                //Smearing to account for Gaussian detector noise (Landau negligible)
                                Evis = _util->GaussianSmearing(Evis, ECAL_MIP_Res);
                                //1 MIP = 0.814 MeV
                                double Erec = Evis * MIP2GeV_factor;
                                erecon.push_back((Erec > 0) ? Erec : 0.);
                                etime.push_back(ecaltime);
                                detected.push_back(1);

                                //Muon/Pions/Protons are reco as Muons (without MuID detector)
                                if( abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 ) {
                                    recopidecal.push_back(13);
                                }
                                else{
                                    recopidecal.push_back(-1);
                                }
                            }
                            else
                            {
                                //Does not reach the ECAL???
                                erecon.push_back(-1);
                                recopidecal.push_back(-1);
                                etime.push_back(-1);
                                detected.push_back(0);
                            }
                        } //end endpoint is not in TPC

                        //end tpc
                        //***************************************************************************************************************/

                        //start pid
                        //***************************************************************************************************************/

                        //--------------------------------------------------------------------------
                        // Start of PID Parametrization
                        //--------------------------------------------------------------------------

                        float p = preco;
                        // read the PID parametrization ntuple from T. Junk
                        // TString filename="pid.root";

                        std::vector<double> vec;
                        std::vector<std::string> pnamelist     = {"#pi", "#mu", "p", "K", "d", "e"};
                        std::vector<std::string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};

                        int qclosest = 0;
                        float dist = 100000000.;

                        // // std::cout << "preco " << p << " qclosest " << qclosest << std::endl;

                        for (int q = 0; q < 501; ++q)
                        {
                            //Check the title and the reco momentum take only the one that fits
                            std::string fulltitle = m_pidinterp[q]->GetTitle();
                            unsigned first = fulltitle.find("=");
                            unsigned last = fulltitle.find("GeV");
                            std::string substr = fulltitle.substr(first+1, last - first-1);
                            float pidinterp_mom = std::atof(substr.c_str());
                            //calculate the distance between the bin and mom, store the q the closest
                            float disttemp = std::abs(pidinterp_mom - p);
                            //// std::cout << disttemp << " " << dist << std::endl;

                            // // std::cout << "preco " << p << " ptitle " << pidinterp_mom << " dist " << disttemp << " q " << q << std::endl;

                            if( disttemp < dist ){
                                dist = disttemp;
                                qclosest = q;
                                // // std::cout << "pid mom " << pidinterp_mom << " reco mom " << p << " dist " << dist << " qclosest " << qclosest << std::endl;
                            }
                        } // closes the "pidmatrix" loop

                        // // std::cout << "Started pid" << std::endl;
                        //loop over the columns (true pid)
                        std::vector< std::pair<float, std::string> > v_prob;

                        // // std::cout << "pidm " << pidm << std::endl;
                        //get true particle name
                        std::string trueparticlename = m_pidinterp[qclosest]->GetXaxis()->GetBinLabel(pidm+1);
                        // // std::cout << trueparticlename << std::endl;
                        if ( trueparticlename == pnamelist[pidm] )
                        {
                            //loop over the rows (reco pid)
                            for (int pidr = 0; pidr < 6; ++pidr)
                            {
                                // // std::cout << "pidr " << pidr << std::endl;
                                std::string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
                                if (recoparticlename == recopnamelist[pidr])
                                {
                                    float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1,pidr+1);
                                    prob_arr.push_back(prob);

                                    // // std::cout << "true part " << trueparticlename << " true pid " << pdg_charged.at(pidm) << " reco name " << recoparticlename << " reco part list "
                                    // << recopnamelist[pidr] <<  " true mom " << ptrue << " reco mom " <<  p << " prob " << pidinterp->GetBinContent(pidm+1,pidr+1) << '\n';
                                    //Need to check random number value and prob value then associate the recopdg to the reco prob
                                    v_prob.push_back( std::make_pair(prob, recoparticlename) );
                                }
                            }

                            int pid = -1;
                            if(v_prob.size() > 1)
                            {
                                //Order the vector of prob
                                std::sort(v_prob.begin(), v_prob.end());
                                //Throw a random number between 0 and 1
                                float random_number = _util->GetRamdomNumber();
                                //Make cumulative sum to get the range
                                std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});

                                for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++)
                                {
                                    if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first )
                                    {
                                        pid = pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) );
                                    }
                                }
                            }
                            else
                            {
                                pid = pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) );
                            }

                            recopid.push_back( pid );
                        } // closes the if statement

                        //end pid
                        //***************************************************************************************************************/

                    } // closes the conditional statement of trueparticlename == MC true pdg
                    else
                    {
                        //not in the pdglist of particles but visible in TPC?
                        auto found = std::find(pdg_charged.begin(), pdg_charged.end(), abs(pdg));
                        if(found == pdg_charged.end())
                        {
                            // // std::cout << "Maybe visible but not {#pi, #mu, p, K, d, e};" << std::endl;
                            // // std::cout << "pdg " << pdg << std::endl;

                            truepdg.push_back(pdg);
                            detected.push_back(0);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            // save the true momentum
                            truep.push_back(ptrue);
                            // save the true angle
                            _angle.push_back(angle);
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(-1);
                            erecon.push_back(-1);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);
                            recopid.push_back(-1);
                            recopidecal.push_back(-1);
                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                            break; //break out of the loop over the vertical bining!
                        }
                    }
                } // closes the vertical bining loop of the pid matrix
            }//close if track_length > tpc_min_length
            else
            {
                //Not visible in the TPC
                //for neutrons
                if(std::abs(pdg) == 2112)
                {
                    //start neutrons
                    //***************************************************************************************************************/

                    if(_util->PointInCalo(epoint)) //needs to stop in the ECAL
                    {
                        //check if it can be detected by the ECAL
                        //Assumes 40% efficiency to detect
                        float random_number = _util->GetRamdomNumber();
                        float true_KE = std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass) - neutron_mass;//KE here
                        // float true_KE = ptrue*ptrue / (2*neutron_mass); // in GeV
                        int index = (true_KE >= 0.05) ? 1 : 0;

                        // // std::cout << "KE " << true_KE << " index " << index << " 1 - eff " << 1-NeutronECAL_detEff[index] << " rdnm " << random_number << std::endl;

                        if(random_number > (1 - NeutronECAL_detEff[index]) && true_KE > 0.003)//Threshold of 3 MeV
                        {
                            //TODO random is first interaction or rescatter and smear accordingly to Chris's study
                            //Detected in the ECAL
                            // recopid.push_back(2112);
                            recopid.push_back(-1); //reco pid set to 0?
                            detected.push_back(1);
                            float eres = sigmaNeutronECAL_first * true_KE;
                            float ereco_KE = _util->GaussianSmearing( true_KE, eres );
                            float ereco = ereco_KE + neutron_mass;
                            erecon.push_back(ereco > 0 ? ereco : 0.);
                            // // std::cout << "true part n true energy " << std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass) << " ereco " << erecon[i] << std::endl;
                            truepdg.push_back(pdg);
                            truep.push_back(ptrue);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(ecaltime);
                            _angle.push_back(angle);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);
                            recopidecal.push_back(2112);
                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                        else
                        {
                            //neutron not detected
                            detected.push_back(0);
                            truep.push_back(ptrue);
                            recopid.push_back(-1);
                            erecon.push_back(-1);
                            truepdg.push_back(pdg);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(-1);
                            _angle.push_back(angle);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);
                            recopidecal.push_back(-1);
                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                    } //endpoint is in ECAL
                    else
                    {
                        //Endpoint is not in calo (TPC/isInBetween or outside Calo)
                        detected.push_back(0);
                        truep.push_back(ptrue);
                        recopid.push_back(-1);
                        erecon.push_back(-1);
                        truepdg.push_back(pdg);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        if(_correct4origin){
                            _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                            _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                            _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                            _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                            _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                            _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                        } else {
                            _MCPStartX.push_back(MCPStartX->at(i));
                            _MCPStartY.push_back(MCPStartY->at(i));
                            _MCPStartZ.push_back(MCPStartZ->at(i));
                            _MCPEndX.push_back(MCPEndX->at(i));
                            _MCPEndY.push_back(MCPEndY->at(i));
                            _MCPEndZ.push_back(MCPEndZ->at(i));
                        }
                        pdgmother.push_back(PDGMother->at(i));
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                        mctime.push_back(time);
                        etime.push_back(-1);
                        _angle.push_back(angle);
                        _preco.push_back(-1);
                        anglereco.push_back(-1);
                        recopidecal.push_back(-1);
                        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                        mctrkid.push_back(MCPTrkID->at(i));
                        motherid.push_back(MCMotherTrkID->at(i));
                    }
                    //End neutrons
                    //***************************************************************************************************************/
                }
                else if(std::abs(pdg) == 111)    //for pi0s
                {
                    //start pi0s
                    //***************************************************************************************************************/
                    erecon.push_back(-1);
                    recopid.push_back(-1);
                    detected.push_back(0);
                    recopidecal.push_back(-1);

                    truep.push_back(ptrue);
                    truepdg.push_back(pdg);
                    truepx.push_back(MCPStartPX->at(i));
                    truepy.push_back(MCPStartPY->at(i));
                    truepz.push_back(MCPStartPZ->at(i));
                    if(_correct4origin){
                        _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                        _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                        _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                        _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                        _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                        _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                    } else {
                        _MCPStartX.push_back(MCPStartX->at(i));
                        _MCPStartY.push_back(MCPStartY->at(i));
                        _MCPStartZ.push_back(MCPStartZ->at(i));
                        _MCPEndX.push_back(MCPEndX->at(i));
                        _MCPEndY.push_back(MCPEndY->at(i));
                        _MCPEndZ.push_back(MCPEndZ->at(i));
                    }
                    pdgmother.push_back(PDGMother->at(i));
                    _angle.push_back(angle);
                    //Save MC process
                    _MCProc.push_back(mcp_process);
                    _MCEndProc.push_back(mcp_endprocess);
                    mctime.push_back(time);
                    etime.push_back(-1);
                    _preco.push_back(-1);
                    anglereco.push_back(-1);

                    for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);

                    mctrkid.push_back(MCPTrkID->at(i));
                    motherid.push_back(MCMotherTrkID->at(i));

                    //end pi0s
                    //***************************************************************************************************************/
                }
                else if(std::abs(pdg) == 22)//for gammas
                {
                    //start gammas
                    //***************************************************************************************************************/

                    if( PDGMother->at(i) != 111 )
                    {
                        //Endpoint is in the ECAL
                        if(_util->PointInCalo(epoint))
                        {
                            //if they hit the ECAL and smear their energy
                            float ECAL_resolution = fRes->Eval(ptrue)*ptrue;
                            float ereco = _util->GaussianSmearing(ptrue, ECAL_resolution);
                            erecon.push_back( (ereco > 0) ? ereco : 0. );
                            recopid.push_back(-1);
                            detected.push_back(1);

                            truep.push_back(ptrue);
                            truepdg.push_back(pdg);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            _angle.push_back(angle);
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(ecaltime);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);

                            //reach the ECAL, should be tagged as gamma
                            recopidecal.push_back(22);

                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                        else if(_util->PointInTPC(epoint) || _util->PointStopBetween(epoint) || _util->isThroughCalo(epoint))
                        {
                            //case endpoint is in the TPC (Converted!) or in between the TPC/ECAL
                            erecon.push_back(-1);
                            recopid.push_back(-1);
                            detected.push_back(0);

                            truep.push_back(ptrue);
                            truepdg.push_back(pdg);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            _angle.push_back(angle);
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(-1);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);
                            //converted so not seen in ECAL
                            recopidecal.push_back(-1);
                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                        else{
                        }
                    }
                    else
                    {
                        //case they are from pi0
                        //Endpoint is not in the tracker, reaches the ecal
                        if(_util->PointInCalo(epoint))
                        {
                            //if they hit the ECAL and smear their energy
                            float ECAL_resolution = fRes->Eval(ptrue)*ptrue;
                            float ereco = _util->GaussianSmearing(ptrue, ECAL_resolution);
                            erecon.push_back((ereco > 0) ? ereco : 0.);
                            recopid.push_back(-1);
                            detected.push_back(1);

                            truep.push_back(ptrue);
                            truepdg.push_back(pdg);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            _angle.push_back(angle);
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(ecaltime);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);

                            //reaches the ecal
                            recopidecal.push_back(22);

                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                        else if(_util->PointInTPC(epoint) || _util->PointStopBetween(epoint) || _util->isThroughCalo(epoint))
                        {
                            //from pi0 and converted in TPC or stopped between TPC/ECAL
                            erecon.push_back(-1);
                            recopid.push_back(-1);
                            detected.push_back(0);

                            truep.push_back(ptrue);
                            truepdg.push_back(pdg);
                            truepx.push_back(MCPStartPX->at(i));
                            truepy.push_back(MCPStartPY->at(i));
                            truepz.push_back(MCPStartPZ->at(i));
                            if(_correct4origin){
                                _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                                _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                                _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                                _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                                _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                                _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                            } else {
                                _MCPStartX.push_back(MCPStartX->at(i));
                                _MCPStartY.push_back(MCPStartY->at(i));
                                _MCPStartZ.push_back(MCPStartZ->at(i));
                                _MCPEndX.push_back(MCPEndX->at(i));
                                _MCPEndY.push_back(MCPEndY->at(i));
                                _MCPEndZ.push_back(MCPEndZ->at(i));
                            }
                            pdgmother.push_back(PDGMother->at(i));
                            _angle.push_back(angle);
                            //Save MC process
                            _MCProc.push_back(mcp_process);
                            _MCEndProc.push_back(mcp_endprocess);
                            mctime.push_back(time);
                            etime.push_back(-1);
                            _preco.push_back(-1);
                            anglereco.push_back(-1);
                            //converted not seen by ecal
                            recopidecal.push_back(-1);
                            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                            mctrkid.push_back(MCPTrkID->at(i));
                            motherid.push_back(MCMotherTrkID->at(i));
                        }
                        else{
                        }
                    }
                    //end gammas
                    //***************************************************************************************************************/
                }
                else if(std::find(neutrinos.begin(), neutrinos.end(), abs(pdg)) != neutrinos.end())
                {
                    //Case for neutrinos from NC interactions
                    truepdg.push_back(pdg);
                    detected.push_back(0);
                    truepx.push_back(MCPStartPX->at(i));
                    truepy.push_back(MCPStartPY->at(i));
                    truepz.push_back(MCPStartPZ->at(i));
                    if(_correct4origin){
                        _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                        _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                        _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                        _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                        _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                        _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                    } else {
                        _MCPStartX.push_back(MCPStartX->at(i));
                        _MCPStartY.push_back(MCPStartY->at(i));
                        _MCPStartZ.push_back(MCPStartZ->at(i));
                        _MCPEndX.push_back(MCPEndX->at(i));
                        _MCPEndY.push_back(MCPEndY->at(i));
                        _MCPEndZ.push_back(MCPEndZ->at(i));
                    }
                    pdgmother.push_back(PDGMother->at(i));
                    // save the true momentum
                    truep.push_back(ptrue);
                    // save the true angle
                    _angle.push_back(angle);
                    //Save MC process
                    _MCProc.push_back(mcp_process);
                    _MCEndProc.push_back(mcp_endprocess);
                    mctime.push_back(time);

                    etime.push_back(-1);
                    erecon.push_back(0.);
                    recopidecal.push_back(-1);
                    _preco.push_back(-1);
                    anglereco.push_back(-1);
                    recopid.push_back(-1);
                    for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);

                    mctrkid.push_back(MCPTrkID->at(i));
                    motherid.push_back(MCMotherTrkID->at(i));
                }
                else
                {
                    //Case for particles that stop or go through ECAL (problematic particles with no track length????)
                    //Not visible in the TPC and not neutron or gamma or pi0 (otherwise it has been already done above)

                    if(_util->PointInCalo(epoint))
                    {
                        truepdg.push_back(pdg);
                        detected.push_back(1);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        if(_correct4origin){
                            _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                            _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                            _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                            _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                            _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                            _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                        } else {
                            _MCPStartX.push_back(MCPStartX->at(i));
                            _MCPStartY.push_back(MCPStartY->at(i));
                            _MCPStartZ.push_back(MCPStartZ->at(i));
                            _MCPEndX.push_back(MCPEndX->at(i));
                            _MCPEndY.push_back(MCPEndY->at(i));
                            _MCPEndZ.push_back(MCPEndZ->at(i));
                        }
                        pdgmother.push_back(PDGMother->at(i));
                        // save the true momentum
                        truep.push_back(ptrue);
                        // save the true angle
                        _angle.push_back(angle);
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                        mctime.push_back(time);

                        etime.push_back(ecaltime);

                        TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(pdg));
                        float mass = 0.;
                        if(nullptr != part)
                        mass = part->Mass();//in GeV

                        float etrue = std::sqrt(ptrue*ptrue + mass*mass);
                        float ECAL_resolution = fRes->Eval(etrue)*etrue;
                        float ereco = _util->GaussianSmearing(etrue, ECAL_resolution);
                        erecon.push_back((ereco > 0) ? ereco : 0.);

                        //Electron
                        if( abs(pdg) == 11 ){
                            recopidecal.push_back(11);
                        }
                        else if( abs(pdg) == 13 || abs(pdg) == 211 )
                        {
                            //Muons and Pions
                            //ptrue < 480 MeV/c 100% separation
                            //80% from 480 to 750
                            //90% up to 750 to 900
                            //95% over 900
                            float random_number = _util->GetRamdomNumber();

                            if(ptrue < 0.48) {
                                recopidecal.push_back(abs(pdg));//100% efficiency by range
                            }
                            else if(ptrue >= 0.48 && ptrue < 0.75)
                            {
                                //case muon
                                if(abs(pdg) == 13)
                                {
                                    if(random_number > (1 - 0.8)) {
                                        recopidecal.push_back(13);
                                    }
                                    else{
                                        recopidecal.push_back(211);
                                    }
                                }

                                //case pion
                                if(abs(pdg) == 211)
                                {
                                    if(random_number > (1 - 0.8)) {
                                        recopidecal.push_back(211);
                                    }
                                    else{
                                        recopidecal.push_back(13);
                                    }
                                }
                            }
                            else if(ptrue >= 0.75 && ptrue < 0.9)
                            {
                                //case muon
                                if(abs(pdg) == 13){
                                    if(random_number > (1 - 0.9)) {
                                        recopidecal.push_back(13);
                                    }
                                    else{
                                        recopidecal.push_back(211);
                                    }
                                }
                                //case pion
                                if(abs(pdg) == 211) {
                                    if(random_number > (1 - 0.9)) {
                                        recopidecal.push_back(211);
                                    }
                                    else{
                                        recopidecal.push_back(13);
                                    }
                                }
                            }
                            else
                            {
                                //case muon
                                if(abs(pdg) == 13){
                                    if(random_number > (1 - 0.95)) {
                                        recopidecal.push_back(13);
                                    }
                                    else{
                                        recopidecal.push_back(211);
                                    }
                                }
                                //case pion
                                if(abs(pdg) == 211){
                                    if(random_number > (1 - 0.95)) {
                                        recopidecal.push_back(211);
                                    }
                                    else{
                                        recopidecal.push_back(13);
                                    }
                                }
                            }
                        }
                        else if( abs(pdg) == 2212 )
                        {
                            recopidecal.push_back(2212);//TODO for p/pi separation
                        }
                        else {
                            recopidecal.push_back(-1);
                        }

                        _preco.push_back(-1);
                        anglereco.push_back(-1);
                        recopid.push_back(-1);
                        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);

                        mctrkid.push_back(MCPTrkID->at(i));
                        motherid.push_back(MCMotherTrkID->at(i));
                    }
                    else if (_util->isThroughCalo(epoint))
                    {
                        truepdg.push_back(pdg);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        if(_correct4origin){
                            _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                            _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                            _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                            _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                            _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                            _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                        } else {
                            _MCPStartX.push_back(MCPStartX->at(i));
                            _MCPStartY.push_back(MCPStartY->at(i));
                            _MCPStartZ.push_back(MCPStartZ->at(i));
                            _MCPEndX.push_back(MCPEndX->at(i));
                            _MCPEndY.push_back(MCPEndY->at(i));
                            _MCPEndZ.push_back(MCPEndZ->at(i));
                        }
                        pdgmother.push_back(PDGMother->at(i));
                        // save the true momentum
                        truep.push_back(ptrue);
                        // save the true angle
                        _angle.push_back(angle);
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                        mctime.push_back(time);
                        //Case the endpoint is outside the CALO -> it went through the ECAL (mu/pi/p possible)
                        //the ECAL will see 60 MIPs on average
                        double Evis = (double)nLayers; //in MIP
                        //Smearing to account for Gaussian detector noise (Landau negligible)
                        Evis = _util->GaussianSmearing(Evis, ECAL_MIP_Res);
                        //1 MIP = 0.814 MeV
                        double Erec = Evis * MIP2GeV_factor;
                        erecon.push_back((Erec > 0) ? Erec : 0.);
                        etime.push_back(ecaltime);
                        detected.push_back(1);

                        //Muon/Pions/Protons are reco as Muons (without MuID detector)
                        if( abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 ) {
                            recopidecal.push_back(13);
                        }
                        else {
                            recopidecal.push_back(-1);
                        }

                        _preco.push_back(-1);
                        anglereco.push_back(-1);
                        recopid.push_back(-1);
                        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);

                        mctrkid.push_back(MCPTrkID->at(i));
                        motherid.push_back(MCMotherTrkID->at(i));
                    }
                    else if(_util->PointInTPC(epoint) || _util->PointStopBetween(epoint))
                    {
                        truepdg.push_back(pdg);
                        detected.push_back(0);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        if(_correct4origin){
                            _MCPStartX.push_back(MCPStartX->at(i) - _util->GetOrigin()[0]);
                            _MCPStartY.push_back(MCPStartY->at(i) - _util->GetOrigin()[1]);
                            _MCPStartZ.push_back(MCPStartZ->at(i) - _util->GetOrigin()[2]);
                            _MCPEndX.push_back(MCPEndX->at(i) - _util->GetOrigin()[0]);
                            _MCPEndY.push_back(MCPEndY->at(i) - _util->GetOrigin()[1]);
                            _MCPEndZ.push_back(MCPEndZ->at(i) - _util->GetOrigin()[2]);
                        } else {
                            _MCPStartX.push_back(MCPStartX->at(i));
                            _MCPStartY.push_back(MCPStartY->at(i));
                            _MCPStartZ.push_back(MCPStartZ->at(i));
                            _MCPEndX.push_back(MCPEndX->at(i));
                            _MCPEndY.push_back(MCPEndY->at(i));
                            _MCPEndZ.push_back(MCPEndZ->at(i));
                        }
                        pdgmother.push_back(PDGMother->at(i));
                        // save the true momentum
                        truep.push_back(ptrue);
                        // save the true angle
                        _angle.push_back(angle);
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                        mctime.push_back(time);
                        etime.push_back(-1);
                        erecon.push_back(-1);
                        _preco.push_back(-1);
                        anglereco.push_back(-1);
                        recopid.push_back(-1);
                        recopidecal.push_back(-1);
                        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
                        mctrkid.push_back(MCPTrkID->at(i));
                        motherid.push_back(MCMotherTrkID->at(i));
                    }
                    else {

                    }
                }// end is not neutron, pi0 or gamma
            }// end not visible in TPC
            nFSP++;
        } // closes the MC truth loop

        _nFSP.push_back(nFSP);

        //Check if vectors have good size
        if(this->CheckVectorSize()) {
            this->FillTTree();
        } else {
            std::cerr << "Event " << _Event << std::endl;
            std::cerr << "Number of FSP " << nFSP << std::endl;

            std::cerr << "Size of pdgmother " << pdgmother.size() << std::endl;
            std::cerr << "Size of truepdg " << truepdg.size() << std::endl;
            std::cerr << "Size of mctime " << mctime.size() << std::endl;
            std::cerr << "Size of mctrkid " << mctrkid.size() << std::endl;
            std::cerr << "Size of motherid " << motherid.size() << std::endl;
            std::cerr << "Size of _MCPStartX " << _MCPStartX.size() << std::endl;
            std::cerr << "Size of _MCPStartY " << _MCPStartY.size() << std::endl;
            std::cerr << "Size of _MCPStartZ " << _MCPStartZ.size() << std::endl;
            std::cerr << "Size of _MCPEndX " << _MCPEndX.size() << std::endl;
            std::cerr << "Size of _MCPEndY " << _MCPEndY.size() << std::endl;
            std::cerr << "Size of _MCPEndZ " << _MCPEndZ.size() << std::endl;
            std::cerr << "Size of _MCProc " << _MCProc.size() << std::endl;
            std::cerr << "Size of _MCEndProc " << _MCEndProc.size() << std::endl;
            std::cerr << "Size of trkLen " << trkLen.size() << std::endl;
            std::cerr << "Size of trkLenPerp " << trkLenPerp.size() << std::endl;
            std::cerr << "Size of truep " << truep.size() << std::endl;
            std::cerr << "Size of truepx " << truepx.size() << std::endl;
            std::cerr << "Size of truepy " << truepy.size() << std::endl;
            std::cerr << "Size of truepz " << truepz.size() << std::endl;
            std::cerr << "Size of _angle " << _angle.size() << std::endl;

            //Reco values
            std::cerr << "Size of detected " << detected.size() << std::endl;
            std::cerr << "Size of recopid " << recopid.size() << std::endl;
            std::cerr << "Size of recopidecal " << recopidecal.size() << std::endl;
            // std::cerr << "Size of prob_arr " << prob_arr.size() << std::endl;
            std::cerr << "Size of anglereco " << anglereco.size() << std::endl;
            std::cerr << "Size of _preco " << _preco.size() << std::endl;
            std::cerr << "Size of erecon " << erecon.size() << std::endl;
            std::cerr << "Size of etime " << etime.size() << std::endl;

            //Geometry
            std::cerr << "Size of isFidStart " << isFidStart.size() << std::endl;
            std::cerr << "Size of isTPCStart " << isTPCStart.size() << std::endl;
            std::cerr << "Size of isCaloStart " << isCaloStart.size() << std::endl;
            std::cerr << "Size of isThroughCaloStart " << isThroughCaloStart.size() << std::endl;
            std::cerr << "Size of isInBetweenStart " << isInBetweenStart.size() << std::endl;
            std::cerr << "Size of isBarrelStart " << isBarrelStart.size() << std::endl;
            std::cerr << "Size of isEndcapStart " << isEndcapStart.size() << std::endl;

            std::cerr << "Size of isFidEnd " << isFidEnd.size() << std::endl;
            std::cerr << "Size of isTPCEnd " << isTPCEnd.size() << std::endl;
            std::cerr << "Size of isCaloEnd " << isCaloEnd.size() << std::endl;
            std::cerr << "Size of isThroughCaloEnd " << isThroughCaloEnd.size() << std::endl;
            std::cerr << "Size of isInBetweenEnd " << isInBetweenEnd.size() << std::endl;
            std::cerr << "Size of isBarrelEnd " << isBarrelEnd.size() << std::endl;
            std::cerr << "Size of isEndcapEnd " << isEndcapEnd.size() << std::endl;

            std::cerr << "Event with wrong vector sizes... skipped" << std::endl;
        }
    } // closes the event loop
} // closes the main loop function
