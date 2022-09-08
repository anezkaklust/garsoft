////////////////////////////////////////////////////////////////////////////////
// Class:       anatest
// Plugin Type: analyzer (art v2_11_02)
// File:        anatest_module.cc
//
// Generated at Tue Mar 3 11:11:11 2020 by Leo Bellantoni
//
// a plugin for testing whatever needs be tested today.
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



#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/BitFieldCoder.h"



#include "TTree.h"



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

        // Position of TPC from geometry service; 1 S Boston Ave.
        double ItsInTulsa[3];

        //Geometry
        const geo::GeometryCore* fGeo; ///< pointer to the geometry
/*        std::string fECALEncoding;
        std::string fMuIDEncoding;
         geo::BitFieldCoder *fFieldDecoder_ECAL;
         geo::BitFieldCoder *fFieldDecoder_MuID;*/



        // Output tree
        TTree *fTree;

        // global event info
        Int_t fEvent;
        Int_t fRun;
        Int_t fSubRun;



        // Decipher the calorimeter geometry
        // Positions in this node
        Float_t fX;
        Float_t fY;
        Float_t fZ;

        unsigned int fECALorMuID, fSystem, fModule, fStave, fLayer, fSlice;

        unsigned int getECALorMuIDNumber(std::string volname) const;
        unsigned int getSystemNumber    (std::string volname) const;
        unsigned int getModuleNumber    (std::string volname) const;
        unsigned int getStaveNumber     (std::string volname) const;
        unsigned int getLayerNumber     (std::string volname) const;
        unsigned int getSliceNumber     (std::string volname) const;

    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatest::anatest(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    fGeo     = providerFrom<geo::GeometryGAr>();


/*    fECALEncoding = fGeo->GetECALCellIDEncoding();
      fFieldDecoder_ECAL = new gar::geo::BitFieldCoder( fECALEncoding );

    if (fGeo->HasMuonDetector()) {
        fMuIDEncoding = fGeo->GetMuIDCellIDEncoding();
        fFieldDecoder_MuID = new gar::geo::BitFieldCoder( fMuIDEncoding );
    }*/




    fECALorMuID = p.get<unsigned int>("PickECALorMuID");
    fSystem     = p.get<unsigned int>("PickSystem");
    fModule     = p.get<unsigned int>("PickModule");
    fStave      = p.get<unsigned int>("PickStave");
    fLayer      = p.get<unsigned int>("PickLayer");
    fSlice      = p.get<unsigned int>("PickSlice");



    fTree   = nullptr;

    return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::beginJob() {

    ItsInTulsa[0] = fGeo->TPCXCent();
    ItsInTulsa[1] = fGeo->TPCYCent();
    ItsInTulsa[2] = fGeo->TPCZCent();



    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("MPDtestTree","MPDtestTree");



    fTree->Branch("Run",      &fRun,      "Run/I");
    fTree->Branch("SubRun",   &fSubRun,   "SubRun/I");
    fTree->Branch("Event",    &fEvent,    "Event/I");
    
    fTree->Branch("X",        &fX,        "X/F");
    fTree->Branch("Y",        &fY,        "Y/F");
    fTree->Branch("Z",        &fZ,        "Z/F");

    return;
}  // End of :anatest::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::analyze(art::Event const & e) {



    fRun    = e.run();            // Hardly matters right now
    fSubRun = e.subRun();
    fEvent  = e.id().event();


    float Xrange = floor(fGeo->GetMPDHalfWidth());
    float Yrange = floor(fGeo->GetMPDHalfHeight());
    float Zrange = floor(fGeo->GetMPDLength()/2.0);

    for (float Xo = -Xrange; Xo <= +Xrange; Xo+=1.0){
        for (float Yo  = -Yrange; Yo <= +Yrange; Yo+=1.0){
            for (float Zo  = -Zrange; Zo <= +Zrange; Zo+=1.0){

                float X = Xo +ItsInTulsa[0];
                float Y = Yo +ItsInTulsa[1];
                float Z = Zo +ItsInTulsa[2];

                TGeoNode* node = fGeo->FindNode(X, Y, Z);
                if (node==NULL) continue;
                std::string VolumeName  = node->GetVolume()->GetName();

                unsigned int lECALorMuID = getECALorMuIDNumber(VolumeName);
                unsigned int lSystem     = getSystemNumber(VolumeName);
                unsigned int lModule     = getModuleNumber(VolumeName);
                unsigned int lStave      = getStaveNumber(VolumeName);
                // Don't select on Layer or Slice, for now.
                // unsigned int lLayer      = getLayerNumber(VolumeName);
                // unsigned int lSlice      = getSliceNumber(VolumeName);

               if (lECALorMuID==fECALorMuID && lSystem==fSystem && 
                    lModule==fModule         && lStave==fStave) {
                        fX = Xo;    fY = Yo;    fZ = Zo;
                        fTree->Fill();
                }
            }
        }
    }
    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
unsigned int gar::anatest::getECALorMuIDNumber(std::string volname) const {
    unsigned int retval = 0;
    if( volname.find("ECAL") !=  std::string::npos ) retval = 1;
    if( volname.find("ECal") !=  std::string::npos ) retval = 1;
    if( volname.find("Yoke") !=  std::string::npos ) retval = 2;
    return retval;
}

//------------------------------------------------------------------------------
unsigned int gar::anatest::getSystemNumber(std::string volname) const {
    unsigned int det_id = 0;
    if( volname.find("Barrel") !=  std::string::npos ) det_id = 1;
    if( volname.find("Endcap") !=  std::string::npos ) det_id = 2;
    return det_id;
}

//------------------------------------------------------------------------------
unsigned int gar::anatest::getModuleNumber(std::string volname) const {
    return std::atoi( (volname.substr( volname.find("_module") + 7, 2)).c_str() );
}

//------------------------------------------------------------------------------
unsigned int gar::anatest::getStaveNumber(std::string volname) const {
    return std::atoi( (volname.substr( volname.find("_stave") + 6, 2)).c_str() );
}

//------------------------------------------------------------------------------
unsigned int gar::anatest::getLayerNumber(std::string volname) const {
    return std::atoi( (volname.substr( volname.find("layer") + 5, 2)).c_str() );
}

//------------------------------------------------------------------------------
unsigned int gar::anatest::getSliceNumber(std::string volname) const {
    return std::atoi( (volname.substr( volname.find("slice") + 5, 1)).c_str() );
}


DEFINE_ART_MODULE(gar::anatest)
