#ifndef CONVEDEP2ART_H
#define CONVEDEP2ART_H

// C++ Includes
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <limits>
#include <exception>

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
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "art/Persistency/Common/PtrMaker.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "nutools/EventGeneratorBase/GENIE/EVGBAssociationUtil.h"

//GENIE
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "PDG/PDGLibrary.h"

// GArSoft Includes
// #include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/Geometry.h"
#include "Geometry/LocalTransformation.h"
#include "Geometry/BitFieldCoder.h"
#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorProperties.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/ECALProperties.h"
#include "DetectorInfo/ECALPropertiesService.h"

// ROOT Includes
#include "TChain.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"

//Edep-Sim includes
#include "Utilities/Converter/edep-io/TG4Event.h"
#include "Utilities/Converter/edep-io/TG4HitSegment.h"
#include "Utilities/Converter/edep-io/TG4PrimaryVertex.h"
#include "Utilities/Converter/edep-io/TG4Trajectory.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

namespace util {

    class ConvertEdep2Art : public ::art::EDProducer{
    public:

        /// Standard constructor and destructor for an FMWK module.
        explicit ConvertEdep2Art(fhicl::ParameterSet const& pset);
        virtual ~ConvertEdep2Art();

        /// The main routine of this module: Fetch the primary particles
        /// from the event, simulate their evolution in the detctor, and
        /// produce the detector response.
        void produce (::art::Event& evt);
        void beginJob();
        void beginRun(::art::Run& run);

    private:

        unsigned int GetDetNumber(std::string volname) const;
        unsigned int GetStaveNumber(std::string volname) const;
        unsigned int GetModuleNumber(std::string volname) const;
        unsigned int GetLayerNumber(std::string volname) const;
        unsigned int GetSliceNumber(std::string volname) const;

        double VisibleEnergyDeposition(const TG4HitSegment *hit, bool applyBirks) const;
        bool CheckProcess( std::string process_name ) const;
        unsigned int GetParentage( unsigned int trkid ) const;
        void AddHits(const std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits);
        void AddHitsMinerva(const std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits);

        std::map<int, size_t> TrackIDToMCTruthIndexMap() const { return fTrackIDToMCTruthIndex; }

        bool isMCPMatch(const simb::MCParticle& p1, const simb::MCParticle& p2) const;
        size_t FindMCTruthIndex(std::vector<simb::MCTruth> *mctruthcol, const simb::MCParticle& part) const;

        std::string fEDepSimfile;
        std::string fGhepfile;
        bool fOverlay;
        std::string fECALMaterial;
        std::string fTPCMaterial;
        double fEnergyCut;
        bool fApplyBirks;

        TChain* fTreeChain;
        unsigned int nEntries;

        TChain* fGTreeChain;
        unsigned int nEntriesGhep;
        genie::NtpMCEventRecord*                 fMCRec;

        TG4Event* fEvent;
        double fBirksCoeff;

        bool fkeepEMShowers;
        std::string fEMShowerDaughterMatRegex;   ///< keep EM shower daughters only in these materials

        genie::PDGLibrary *pdglib;
        const gar::geo::GeometryCore*            fGeo;               ///< geometry information
        const gar::detinfo::ECALProperties*      fEcalProp;

        sim::ParticleList* fParticleList;
        bool fHasGHEP;

        std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_ECALDeposits;
        std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_TrackerDeposits;
        std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_MuIDDeposits;

        std::vector<gar::sdp::CaloDeposit> fECALDeposits;
        std::vector<gar::sdp::EnergyDeposit> fGArDeposits;
        std::vector<gar::sdp::CaloDeposit> fTrackerDeposits;
        std::vector<gar::sdp::CaloDeposit> fMuIDDeposits;

        int fSpillCount;
        std::vector<int> fStartSpill;
        std::vector<int> fStopSpill;

        std::map<unsigned int, unsigned int> fTrkIDParent;
        std::map<int, size_t> fTrackIDToMCTruthIndex;
    };

} // namespace util

namespace util {

    //----------------------------------------------------------------------
    // Constructor
    ConvertEdep2Art::ConvertEdep2Art(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset},
    fEDepSimfile( pset.get< std::string >("EDepSimFile", "") ),
    fGhepfile( pset.get< std::string >("GhepFile", "") ),
    fOverlay( pset.get< bool >("OverlayFile", true) ),
    fECALMaterial( pset.get< std::string >("ECALMaterial", "Scintillator") ),
    fTPCMaterial( pset.get< std::string >("TPCMaterial", "HP_ArCH4") ),
    fEnergyCut( pset.get< double >("EnergyCut", 0.000001) ),
    fApplyBirks( pset.get< bool >("applyBirks", true) ),
    fTreeChain(new TChain("EDepSimEvents")),
    fGTreeChain(new TChain("gtree")),
    fMCRec(nullptr),
    fEvent(nullptr),
    fkeepEMShowers( pset.get< bool >("keepEMShowers", true) ),
    fEMShowerDaughterMatRegex( pset.get< std::string >("EMShowerDaughterMatRegex", ".*") ),
    fParticleList(new sim::ParticleList()),
    fHasGHEP(false)
    {
        pdglib = genie::PDGLibrary::Instance();

        if(fEDepSimfile.empty())
        {
            throw cet::exception("ConvertEdep2Art")
            << "Empty edep-sim file";
        }

        //If ghep file is provided
        if(not fGhepfile.empty()) {
            fGTreeChain->Add(fGhepfile.c_str());
            nEntriesGhep = fTreeChain->GetEntries();
            fGTreeChain->SetBranchAddress("gmcrec", &fMCRec);
            fHasGHEP = true;
        }

        fTreeChain->Add(fEDepSimfile.c_str());
        nEntries = fTreeChain->GetEntries();
        fTreeChain->SetBranchAddress("Event", &fEvent);

        fGeo = gar::providerFrom<gar::geo::Geometry>();
        fEcalProp = gar::providerFrom<gar::detinfo::ECALPropertiesService>();
        fBirksCoeff = fEcalProp->ScintBirksConstant(); //CLHEP::mm/CLHEP::MeV;

        fSpillCount = 0;

        if(!fkeepEMShowers){
            LOG_DEBUG("ConvertEdep2Art")
            << " Will not keep EM shower daughters!";
        }

        produces< std::vector<simb::MCTruth> >();
        if(fHasGHEP) {
            produces< std::vector<simb::GTruth>  >();
            produces< art::Assns<simb::MCTruth, simb::GTruth> >();
        }
        produces< art::Assns<simb::MCTruth, simb::MCParticle> >();
        produces< std::vector<simb::MCParticle> >();

        produces< std::vector<gar::sdp::EnergyDeposit> >();
        produces< std::vector<gar::sdp::CaloDeposit> >("ECAL");
        produces< std::vector<gar::sdp::CaloDeposit> >("TrackerSc");
        produces< std::vector<gar::sdp::CaloDeposit> >("MuID");
        // produces< std::vector<gar::sdp::LArDeposit> >();

        produces< art::Assns<gar::sdp::EnergyDeposit, simb::MCParticle> >();
        produces< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> >("ECAL");
        produces< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> >("TrackerSc");
        produces< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> >("MuID");
        // produces< art::Assns<gar::sdp::LArDeposit, simb::MCParticle> >();
    }

    //----------------------------------------------------------------------
    // Destructor
    ConvertEdep2Art::~ConvertEdep2Art() {
        if ( fTreeChain ) delete fTreeChain;
        if ( fGTreeChain ) delete fGTreeChain;
    }

    //----------------------------------------------------------------------
    void ConvertEdep2Art::beginJob() {
        if(fOverlay && fHasGHEP){
            //need to get where is the end of the spill
            for(int ientry = 0; ientry < fGTreeChain->GetEntries(); ientry++)
            {
                fGTreeChain->GetEntry(ientry);

                genie::NtpMCRecHeader rec_header = fMCRec->hdr;
                genie::EventRecord *event = fMCRec->event;

                LOG_DEBUG("ConvertEdep2Art")
                << rec_header
                << *event;

                //Need to get the Rootino showing the end of spill
                genie::GHepParticle *neutrino = event->Probe();
                if(nullptr == neutrino){
                    if(fSpillCount == 0)
                    fStartSpill.insert(fStartSpill.begin()+fSpillCount, 0);
                    else{
                        fStartSpill.insert(fStartSpill.begin()+fSpillCount, fStopSpill.at(fSpillCount-1)+1);
                    }

                    fStopSpill.insert(fStopSpill.begin()+fSpillCount, ientry);
                    fSpillCount++;
                }
            }

            LOG_INFO("ConvertEdep2Art::beginJob()")
            << "Number of spills in the ghep file "
            << fSpillCount;
        }

        return;
    }

    //--------------------------------------------------------------------------
    void ConvertEdep2Art::beginRun(::art::Run& run) {
        return;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetLayerNumber(std::string volname) const {
        return std::atoi( (volname.substr( volname.find("layer_") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetSliceNumber(std::string volname) const {
        return std::atoi( (volname.substr( volname.find("slice") + 5, 1)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetDetNumber(std::string volname) const {
        unsigned int det_id = 0;
        if( volname.find("Barrel") !=  std::string::npos )
        det_id = 1;
        if( volname.find("Endcap") !=  std::string::npos )
        det_id = 2;
        return det_id;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetStaveNumber(std::string volname) const {
        return std::atoi( (volname.substr( volname.find("_stave") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetModuleNumber(std::string volname) const {
        return std::atoi( (volname.substr( volname.find("_module") + 7, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    double ConvertEdep2Art::VisibleEnergyDeposition(const TG4HitSegment *hit, bool applyBirks) const {
        double edep = hit->GetEnergyDeposit();//MeV
        double niel = hit->GetSecondaryDeposit();//MeV
        double length = hit->GetTrackLength();//mm
        double evis = edep;

        if(applyBirks)
        {
            if(fBirksCoeff > 0)
            {
                double nloss = niel;
                if(nloss < 0.0) nloss = 0.0;
                double eloss = edep - nloss;

                if(eloss < 0.0 || length <= 0.0) {
                    nloss = edep;
                    eloss = 0.0;
                }
                if(eloss > 0.0) { eloss /= (1.0 + fBirksCoeff*eloss/length); }

                if(nloss > 0){
                    nloss /= (1.0 + fBirksCoeff*nloss/length);
                }

                evis = eloss + nloss;
            }
        }

        return evis;//MeV
    }

    //--------------------------------------------------------------------------
    bool ConvertEdep2Art::CheckProcess( std::string process_name ) const {
        bool isEMShowerProcess = false;

        //avoid false positive
        if(process_name.find("photonNuclear") != std::string::npos)
        return isEMShowerProcess;

        if( process_name.find("conv")        != std::string::npos  ||
        process_name.find("LowEnConversion") != std::string::npos  ||
        process_name.find("Pair")            != std::string::npos  ||
        process_name.find("compt")           != std::string::npos  ||
        process_name.find("Compt")           != std::string::npos  ||
        process_name.find("Brem")            != std::string::npos  ||
        process_name.find("phot")            != std::string::npos  ||
        process_name.find("Photo")           != std::string::npos  ||
        process_name.find("hIoni")           != std::string::npos  ||
        process_name.find("eIoni")           != std::string::npos  ||
        process_name.find("ionIoni")         != std::string::npos  ||
        (process_name.find("Ion")            != std::string::npos  && process_name.find("mu") != std::string::npos) ||
        process_name.find("annihil")         != std::string::npos ) {
            isEMShowerProcess = true;
        }


        return isEMShowerProcess;
    }

    //--------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetParentage(unsigned int trkid) const {
        unsigned int parentid = -1;
        std::map<unsigned int, unsigned int>::const_iterator it = fTrkIDParent.find(trkid);
        while ( it != fTrkIDParent.end() )
        {
            LOG_DEBUG("ConvertEdep2Art::GetParentage")
            << "parentage for track id " << trkid
            << " is " << (*it).second;

            parentid = (*it).second;
            it = fTrkIDParent.find(parentid);
        }
        return parentid;
    }

    //------------------------------------------------------------------------------
    void ConvertEdep2Art::AddHits(const std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits)
    {
        //Loop over the hits in the map and add them together

        for(auto const &it : m_Deposits) {

            gar::raw::CellID_t cellID = it.first;
            std::vector<gar::sdp::CaloDeposit> vechit = it.second;
            std::sort(vechit.begin(), vechit.end()); //sort per time

            float esum = 0.;
            float time = vechit.at(0).Time();
            int trackID = vechit.at(0).TrackID();
            double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

            for(auto const &hit : vechit) {
                esum += hit.Energy();
            }

            fDeposits.emplace_back( trackID, time, esum, pos, cellID );
        }
    }

    //------------------------------------------------------------------------------
    void ConvertEdep2Art::AddHitsMinerva(const std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits)
    {
        std::string fEncoding = fGeo->GetMinervaCellIDEncoding();
        gar::geo::BitFieldCoder *fFieldDecoder = new gar::geo::BitFieldCoder( fEncoding );

        //Loop over the hits in the map and add them together
        for(auto const &it : m_Deposits) {
            gar::raw::CellID_t cellID = it.first;

            //Check this cellID
            //if it is a triangle = 0 or 3, add all
            //if it is triangle = 1 or 2, need to look for the complementary cellID, add them to this hit and remove from the deposits (to avoid double hit creation)

            int triangleNb = fFieldDecoder->get(cellID, "triangle");
            if(triangleNb == 0 || triangleNb == 3) {
                std::vector<gar::sdp::CaloDeposit> vechit = it.second;
                std::sort(vechit.begin(), vechit.end()); //sort per time

                float esum = 0.;
                float time = vechit.at(0).Time();
                int trackID = vechit.at(0).TrackID();
                double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

                for(auto const &hit : vechit) {
                    esum += hit.Energy();
                }

                fDeposits.emplace_back( trackID, time, esum, pos, cellID );
            }

            // if(triangleNb == 1 || triangleNb == 2)
            // {
            //     if(triangleNb == 1) {
            //         std::vector<gar::sdp::CaloDeposit> vechit = it.second;
            //         std::sort(vechit.begin(), vechit.end()); //sort per time
            //
            //         float esum = 0.;
            //         float time = vechit.at(0).Time();
            //         int trackID = vechit.at(0).TrackID();
            //         double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };
            //
            //         for(auto const &hit : vechit) {
            //             esum += hit.Energy();
            //         }
            //
            //         gar::raw::CellID_t complementary_cellID = 0.;
            //         auto find = m_Deposits.find( complementary_cellID );
            //         if(find != m_Deposits.end()) {
            //             std::vector<gar::sdp::CaloDeposit> vechit_comp = m_Deposits[complementary_cellID];
            //             std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time
            //
            //             for(auto const &hit : vechit_comp) {
            //                 esum += hit.Energy();
            //             }
            //         }
            //         fDeposits.emplace_back( trackID, time, esum, pos, cellID );
            //     }
            //
            //     if(triangleNb == 2) {
            //         std::vector<gar::sdp::CaloDeposit> vechit = it.second;
            //         std::sort(vechit.begin(), vechit.end()); //sort per time
            //
            //         float esum = 0.;
            //         float time = vechit.at(0).Time();
            //         int trackID = vechit.at(0).TrackID();
            //         double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };
            //
            //         for(auto const &hit : vechit) {
            //             esum += hit.Energy();
            //         }
            //
            //         gar::raw::CellID_t complementary_cellID = 0.;
            //         auto find = m_Deposits.find( complementary_cellID );
            //         if(find != m_Deposits.end()) {
            //             std::vector<gar::sdp::CaloDeposit> vechit_comp = m_Deposits[complementary_cellID];
            //             std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time
            //
            //             for(auto const &hit : vechit_comp) {
            //                 esum += hit.Energy();
            //             }
            //         }
            //         fDeposits.emplace_back( trackID, time, esum, pos, cellID );
            //     }
            // }
        }

        delete fFieldDecoder;
    }

    //------------------------------------------------------------------------------
    bool ConvertEdep2Art::isMCPMatch(const simb::MCParticle& p1, const simb::MCParticle& p2) const
    {
        int ulp = 2;
        if( std::fabs(p1.E() - p2.E()) <= std::numeric_limits<float>::epsilon() * std::fabs(p1.E()+p2.E()) * ulp && p1.PdgCode() == p2.PdgCode() )
        return true;

        return false;
    }

    //------------------------------------------------------------------------------
    size_t ConvertEdep2Art::FindMCTruthIndex(std::vector<simb::MCTruth>* mctruthcol, const simb::MCParticle& part) const {

        size_t mctruthindex = 0;
        bool found = false;
        for(size_t index = 0; index < mctruthcol->size(); index++) {
            simb::MCTruth mctrh = mctruthcol->at(index);
            for( size_t ipart = 0; ipart < (size_t)mctrh.NParticles(); ipart++ ) {
                if( isMCPMatch(mctrh.GetParticle(ipart), part) ) {

                    LOG_DEBUG("ConvertEdep2Art") << "FindMCTruthIndex() \n"
                    << mctrh.GetParticle(ipart) << "\n"
                    << part << "\n";

                    mctruthindex = index;
                    found = true;
                    break;
                }
            }
            if( found ) break;
        }

        return mctruthindex;
    }

    //--------------------------------------------------------------------------
    void ConvertEdep2Art::produce(::art::Event& evt) {
        LOG_DEBUG("ConvertEdep2Art") << "produce()";
        art::EventNumber_t eventnumber = evt.id().event();

        if( eventnumber > nEntries ){
            //stop...
            return;
        }

        std::unique_ptr< std::vector<simb::MCTruth> > mctruthcol( new std::vector<simb::MCTruth> );
        std::unique_ptr< std::vector<simb::GTruth> > gtruthcol( new std::vector<simb::GTruth> );
        std::unique_ptr< art::Assns<simb::MCTruth, simb::GTruth> > tgassn( new art::Assns<simb::MCTruth, simb::GTruth> );

        art::PtrMaker<simb::MCTruth> makeMCTruthPtr(evt);
        std::vector< ::art::Ptr<simb::MCTruth> > mctPtrs;

        //--------------------------------------------------------------------------
        //Get the event
        //Starts at 0, evt starts at 1
        fTreeChain->GetEntry(eventnumber-1);

        //--------------------------------------------------------------------------
        if(fOverlay && fHasGHEP){
            size_t index = 0;
            for(int ientry = fStartSpill[eventnumber-1]; ientry < fStopSpill[eventnumber-1]; ientry++)
            {
                fGTreeChain->GetEntry(ientry);

                genie::NtpMCRecHeader rec_header = fMCRec->hdr;
                genie::EventRecord *event = fMCRec->event;

                genie::GHepParticle *neutrino = event->Probe();
                //avoid rootino events
                if(nullptr == neutrino) continue;

                simb::MCTruth mctruth;
                simb::GTruth  gtruth;

                evgb::FillMCTruth(event, 0., mctruth);
                evgb::FillGTruth(event, gtruth);

                mctruthcol->push_back(mctruth);
                gtruthcol->push_back(gtruth);

                //Make a vector of mctruth art ptr
                art::Ptr<simb::MCTruth> MCTruthPtr = makeMCTruthPtr(mctruthcol->size() - 1);
                mctPtrs.push_back(MCTruthPtr);

                evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
                index++;
            }
        }
        else if(not fOverlay && fHasGHEP){
            //Starts at 0, evt starts at 1
            fGTreeChain->GetEntry(eventnumber-1);

            genie::NtpMCRecHeader rec_header = fMCRec->hdr;
            genie::EventRecord *event = fMCRec->event;

            simb::MCTruth mctruth;
            simb::GTruth  gtruth;

            evgb::FillMCTruth(event, 0., mctruth);
            evgb::FillGTruth(event, gtruth);

            mctruthcol->push_back(mctruth);
            gtruthcol->push_back(gtruth);

            //Make a vector of mctruth art ptr
            art::Ptr<simb::MCTruth> MCTruthPtr = makeMCTruthPtr(mctruthcol->size() - 1);
            mctPtrs.push_back(MCTruthPtr);

            evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
        } else {

            //Case where things are made from particle gun! no ghep file provided. Need to create MCTruth object / GTruth object
            simb::MCTruth truth;
            truth.SetOrigin(simb::kSingleParticle);

            for (std::vector<TG4PrimaryVertex>::const_iterator t = fEvent->Primaries.begin(); t != fEvent->Primaries.end(); ++t)
            {
                TLorentzVector pos(t->Position.X() / CLHEP::cm, t->Position.Y() / CLHEP::cm, t->Position.Z() / CLHEP::cm, t->Position.T());

                for (std::vector<TG4PrimaryParticle>::const_iterator p = t->Particles.begin(); p != t->Particles.end(); ++p) {
                    int trackid = p->GetTrackId();
                    std::string primary("primary");

                    TLorentzVector pvec(p->Momentum.Px() * CLHEP::MeV / CLHEP::GeV, p->Momentum.Py() * CLHEP::MeV / CLHEP::GeV, p->Momentum.Pz() * CLHEP::MeV / CLHEP::GeV, p->Momentum.E() * CLHEP::MeV / CLHEP::GeV);

                    simb::MCParticle part(trackid, p->GetPDGCode(), primary);
                    part.AddTrajectoryPoint(pos, pvec);

                    LOG_DEBUG("ConvertEdep2Art") << "Adding primary particle with "
                    << " momentum " << part.P()
                    << " position " << part.Vx() << " " << part.Vy() << " " << part.Vz();

                    truth.Add(part);
                }
            }

            LOG_DEBUG("ConvertEdep2Art") << "Adding mctruth with "
            << " nParticles " << truth.NParticles()
            << " Origin " << truth.Origin();

            mctruthcol->push_back(truth);

            //Make a vector of mctruth art ptr
            art::Ptr<simb::MCTruth> MCTruthPtr = makeMCTruthPtr(mctruthcol->size() - 1);
            mctPtrs.push_back(MCTruthPtr);
        }

        //-----------------------------------

        std::unique_ptr< std::vector<simb::MCParticle> > partCol( new std::vector<simb::MCParticle> );
        std::unique_ptr< art::Assns<simb::MCTruth, simb::MCParticle> > tpassn( new art::Assns<simb::MCTruth, simb::MCParticle> );

        fParticleList->clear();
        fTrkIDParent.clear();
        fTrackIDToMCTruthIndex.clear();

        for (std::vector<TG4Trajectory>::const_iterator t = fEvent->Trajectories.begin(); t != fEvent->Trajectories.end(); ++t)
        {
            int trackID = t->GetTrackId();
            int parentID = t->GetParentId();
            int pdg = t->GetPDGCode();
            std::string name = t->GetName();
            double mass = pdglib->Find(pdg)->Mass();//in GeV
            std::string process_name = "unknown";

            if(parentID == -1) {
                process_name = "primary";
                parentID = 0;
            }
            else {
                //Get the first point that created this particle
                process_name = t->Points.at(0).GetProcessName();
            }// end if not a primary particle

            LOG_DEBUG("ConvertEdep2Art")
            << " Particle " << name
            << " with pdg " << pdg
            << " trackID " << trackID
            << " parent id " << parentID
            << " created with process [ " << process_name << " ]"
            << " is EM " << CheckProcess( process_name );

            simb::MCParticle *fParticle = new simb::MCParticle(trackID, pdg, process_name, parentID, mass);

            for (std::vector<TG4TrajectoryPoint>::const_iterator p = t->Points.begin(); p != t->Points.end(); ++p)
            {
                TLorentzVector position = p->GetPosition();

                TLorentzVector fourPos(position.X() / CLHEP::cm, position.Y() / CLHEP::cm, position.Z() / CLHEP::cm, position.T() );
                TVector3 momentum = p->GetMomentum();
                double px = momentum.x() * CLHEP::MeV / CLHEP::GeV;
                double py = momentum.y() * CLHEP::MeV / CLHEP::GeV;
                double pz = momentum.z() * CLHEP::MeV / CLHEP::GeV;
                TLorentzVector fourMom(px, py, pz, std::sqrt( px*px + py*py + pz*pz + mass*mass ));
                std::string process = p->GetProcessName();

                if(p == t->Points.begin()) process = "Start";
                fParticle->AddTrajectoryPoint(fourPos, fourMom, process);
            }

            std::string end_process = t->Points.at(t->Points.size()-1).GetProcessName();
            fParticle->SetEndProcess(end_process);

            fParticleList->Add( fParticle );
        }

        size_t mcTruthIndex = 0;
        int fCurrentTrackID = 0;
        //Work on MCParticle list
        for (std::map<int, simb::MCParticle*>::iterator itPart = fParticleList->begin(); itPart != fParticleList->end(); ++itPart)
        {
            simb::MCParticle &part = *(itPart->second);
            int parentID = part.Mother();
            int trackID = part.TrackId();
            fCurrentTrackID = trackID;
            std::string process_name = part.Process();

            if( process_name == "primary" )
            {
                if(not fOverlay) {
                    mcTruthIndex = 0;
                } else {
                    //Need to get the index correctly... Check the list of mctruth and g4 particles, match them according to type and energy and pos?
                    mcTruthIndex = FindMCTruthIndex(mctruthcol.get(), part);
                }
            }
            else {

                TLorentzVector part_start = part.Trajectory().Position(0);
                TGeoNode *node = fGeo->FindNode(part_start.X(), part_start.Y(), part_start.Z());
                std::string material_name = "";

                if(node)
                material_name = node->GetMedium()->GetMaterial()->GetName();

                //Skip the creation of the mcp if it is part of shower (based on process name) and only if it is not matching the material
                //TODO debug as it does not work for anatree at the moment, problem related to pdg mom exception
                //Maybe need to keep relation trackid, parentid
                if(not fkeepEMShowers){
                    bool isEMShowerProcess = CheckProcess( process_name );
                    std::regex const re_material(fEMShowerDaughterMatRegex);
                    if( isEMShowerProcess && not std::regex_match(material_name, re_material) ) {

                        //link trkid and parent id to be able to go back in history only for these and stop at the parent that created the shower
                        fTrkIDParent[trackID] = parentID;
                        fCurrentTrackID = -1 * GetParentage(trackID);

                        LOG_DEBUG("ConvertEdep2Art")
                        << " Skipping EM shower daughter "
                        << " with trackID " << trackID
                        << " with parent id " << parentID
                        << " Ultimate parentage " << GetParentage(trackID)
                        << " created with process [ " << process_name << " ]";

                        fParticleList->Archive(itPart->second);
                        continue;
                    }
                }//end not keep EM Shower particles

                if( not fParticleList->KnownParticle(parentID) ) {
                    // do add the particle to the parent id map
                    // just in case it makes a daughter that we have to track as well
                    fTrkIDParent[trackID] = parentID;
                    int pid = GetParentage(parentID);

                    // if we still can't find the parent in the particle navigator,
                    // we have to give up
                    if( not fParticleList->KnownParticle(pid) ) {
                        LOG_DEBUG("ConvertEdep2Art")
                        << "can't find parent id: "
                        << parentID << " in the particle list, or fTrkIDParent."
                        << " Make " << parentID << " the mother ID for track ID "
                        << fCurrentTrackID << " in the hope that it will aid debugging.";
                    }
                    else
                    parentID = pid;
                }

                // Attempt to find the MCTruth index corresponding to the
                // current particle.  If the fCurrentTrackID is not in the
                // map try the parent ID, if that is not there, throw an
                // exception
                try {
                    if(fTrackIDToMCTruthIndex.count(fCurrentTrackID) > 0 )
                    mcTruthIndex = fTrackIDToMCTruthIndex.at(fCurrentTrackID);
                    else if(fTrackIDToMCTruthIndex.count(parentID) > 0 )
                    mcTruthIndex = fTrackIDToMCTruthIndex.at(parentID);
                }
                catch (std::exception& e) {
                    LOG_DEBUG("ConvertEdep2Art")
                    << "Cannot find MCTruth index for track id "
                    << fCurrentTrackID << " or " << parentID
                    << " exception " << e.what();
                    throw;
                }
            } //end not primary particle

            fTrackIDToMCTruthIndex[fCurrentTrackID] = mcTruthIndex;
        }

        // Make link between MCTruth and MCParticles
        size_t nGeneratedParticles = 0;
        const std::map<int, size_t> fTrackIDToMCTruthIndex = this->TrackIDToMCTruthIndexMap(); //Need to make trackID to MCTruth index map

        for (std::map<int, simb::MCParticle*>::iterator itPart = fParticleList->begin(); itPart != fParticleList->end(); ++itPart)
        {
            simb::MCParticle& p = *(itPart->second);

            LOG_DEBUG("ConvertEdep2Art")
            << "adding mc particle with track id: "
            << p.TrackId();

            int trackID = p.TrackId();
            partCol->push_back(std::move(p));

            try {
                if( fTrackIDToMCTruthIndex.count(trackID) > 0) {
                    size_t mctidx = fTrackIDToMCTruthIndex.find(trackID)->second;
                    evgb::util::CreateAssn(*this, evt, *partCol, mctPtrs.at(mctidx), *tpassn, nGeneratedParticles);
                }
            }
            catch ( std::exception& e ) {
                LOG_DEBUG("ConvertEdep2Art")
                << "Cannot find MCTruth for Track Id: " << trackID
                << " to create association between Particle and MCTruth"
                << " exception " << e.what();
                throw;
            }

            fParticleList->Archive(itPart->second);
            ++nGeneratedParticles;
        }

        LOG_DEBUG("ConvertEdep2Art") << "Finished linking MCTruth and MCParticles";

        //--------------------------------------------------------------------------
        std::unique_ptr< std::vector< gar::sdp::EnergyDeposit>  > TPCCol(new std::vector<gar::sdp::EnergyDeposit> );
        std::unique_ptr< std::vector< gar::sdp::CaloDeposit > > ECALCol(new std::vector<gar::sdp::CaloDeposit> );
        std::unique_ptr< std::vector< gar::sdp::CaloDeposit > > TrackerCol(new std::vector<gar::sdp::CaloDeposit> );
        std::unique_ptr< std::vector< gar::sdp::CaloDeposit > > MuIDCol(new std::vector<gar::sdp::CaloDeposit> );
        std::unique_ptr< std::vector< gar::sdp::LArDeposit > > LArCol(new std::vector<gar::sdp::LArDeposit> );

        std::unique_ptr< art::Assns<gar::sdp::EnergyDeposit, simb::MCParticle> > ghmcassn(new art::Assns<gar::sdp::EnergyDeposit, simb::MCParticle>);
        std::unique_ptr< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> > ehmcassn(new art::Assns<gar::sdp::CaloDeposit, simb::MCParticle>); //ECAL
        std::unique_ptr< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> > thmcassn(new art::Assns<gar::sdp::CaloDeposit, simb::MCParticle>); //TrackerSc
        std::unique_ptr< art::Assns<gar::sdp::CaloDeposit, simb::MCParticle> > mhmcassn(new art::Assns<gar::sdp::CaloDeposit, simb::MCParticle>); //MuID
        // std::unique_ptr< art::Assns<gar::sdp::LArDeposit, simb::MCParticle> > lhmcassn(new art::Assns<gar::sdp::LArDeposit, simb::MCParticle>); //LAr

        m_ECALDeposits.clear();
        m_TrackerDeposits.clear();
        m_MuIDDeposits.clear();
        fGArDeposits.clear();
        fECALDeposits.clear();

        //Fill simulated hits
        for (auto d = fEvent->SegmentDetectors.begin(); d != fEvent->SegmentDetectors.end(); ++d)
        {
            if( d->first == "TPC_Drift1" || d->first == "TPC_Drift2" )
            {
                //GAr deposits
                for (std::vector<TG4HitSegment>::const_iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    const TG4HitSegment *hit = &(*h);

                    int trackID = hit->GetPrimaryId();
                    double edep = hit->GetEnergyDeposit() * CLHEP::MeV / CLHEP::GeV;
                    double time = (hit->GetStart().T() + hit->GetStop().T())/2 / CLHEP::s;
                    double x = (hit->GetStart().X() + hit->GetStop().X())/2 /CLHEP::cm;
                    double y = (hit->GetStart().Y() + hit->GetStop().Y())/2 /CLHEP::cm;
                    double z = (hit->GetStart().Z() + hit->GetStop().Z())/2 /CLHEP::cm;
                    double stepLength = hit->GetTrackLength() /CLHEP::cm;

                    if(edep == 0 || edep < fEnergyCut)
                    continue;

                    TGeoNode *node = fGeo->FindNode(x, y, z);//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fTPCMaterial)) ) continue;

                    fGArDeposits.emplace_back(trackID, time, edep, x, y, z, stepLength, (trackID > 0));
                }
            }
            else if( d->first == "BarrelECal_vol" || d->first == "EndcapECal_vol"){
                //ECAL deposits
                for (std::vector<TG4HitSegment>::const_iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    const TG4HitSegment *hit = &(*h);

                    int trackID = hit->GetPrimaryId();
                    double edep = VisibleEnergyDeposition(hit, fApplyBirks) * CLHEP::MeV / CLHEP::GeV;
                    double time = (hit->GetStart().T() + hit->GetStop().T())/2 / CLHEP::s;
                    double x = (hit->GetStart().X() + hit->GetStop().X())/2 /CLHEP::cm;
                    double y = (hit->GetStart().Y() + hit->GetStop().Y())/2 /CLHEP::cm;
                    double z = (hit->GetStart().Z() + hit->GetStop().Z())/2 /CLHEP::cm;

                    if(edep == 0 || edep < fEnergyCut)
                    continue;

                    //Check if it is in the active material of the ECAL
                    TGeoNode *node = fGeo->FindNode(x, y, z);//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) continue;

                    unsigned int layer = GetLayerNumber(VolumeName); //get layer number
                    unsigned int slice = GetSliceNumber(VolumeName); // get slice number
                    unsigned int det_id = GetDetNumber(VolumeName); // 1 == Barrel, 2 = Endcap
                    unsigned int stave = GetStaveNumber(VolumeName); //get the stave number
                    unsigned int module = GetModuleNumber(VolumeName); //get the module number

                    std::array<double, 3> GlobalPosCM = {x, y, z};
                    std::array<double, 3> LocalPosCM;
                    gar::geo::LocalTransformation<TGeoHMatrix> trans;
                    fGeo->WorldToLocal(GlobalPosCM, LocalPosCM, trans);

                    LOG_DEBUG("ConvertEdep2Art")
                    << "Sensitive volume " << d->first
                    << " Hit " << hit
                    << " in volume " << VolumeName
                    << " in material " << volmaterial
                    << " det_id " << det_id
                    << " module " << module
                    << " stave " << stave
                    << " layer " << layer
                    << " slice " << slice;

                    gar::raw::CellID_t cellID = fGeo->GetCellID(node, det_id, stave, module, layer, slice, LocalPosCM);//encoding the cellID on 64 bits

                    double G4Pos[3] = {0., 0., 0.}; // in cm
                    G4Pos[0] = GlobalPosCM[0];
                    G4Pos[1] = GlobalPosCM[1];
                    G4Pos[2] = GlobalPosCM[2];

                    gar::sdp::CaloDeposit calohit( trackID, time, edep, G4Pos, cellID );
                    if(m_ECALDeposits.find(cellID) != m_ECALDeposits.end())
                    m_ECALDeposits[cellID].push_back(calohit);
                    else {
                        std::vector<gar::sdp::CaloDeposit> vechit;
                        vechit.push_back(calohit);
                        m_ECALDeposits.emplace(cellID, vechit);
                    }
                }
            }
            else if( d->first == "Tracker_vol" ) {
                //Minerva Style Sc Tracker for temporary det -> triangle of base 4 cm and height 2 cm
                for (std::vector<TG4HitSegment>::const_iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    const TG4HitSegment *hit = &(*h);

                    int trackID = hit->GetPrimaryId();
                    double edep = VisibleEnergyDeposition(hit, fApplyBirks) * CLHEP::MeV / CLHEP::GeV;
                    double time = (hit->GetStart().T() + hit->GetStop().T())/2 / CLHEP::s;
                    double x = (hit->GetStart().X() + hit->GetStop().X())/2 /CLHEP::cm;
                    double y = (hit->GetStart().Y() + hit->GetStop().Y())/2 /CLHEP::cm;
                    double z = (hit->GetStart().Z() + hit->GetStop().Z())/2 /CLHEP::cm;

                    if(edep == 0 || edep < fEnergyCut)
                    continue;

                    //Check if it is in the active material of the ECAL
                    TGeoNode *node = fGeo->FindNode(x, y, z);//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) continue;

                    unsigned int layer = GetLayerNumber(VolumeName); //get layer number
                    unsigned int slice = GetSliceNumber(VolumeName); // get slice number
                    unsigned int det_id = 3;

                    std::array<double, 3> GlobalPosCM = {x, y, z};
                    std::array<double, 3> LocalPosCM;
                    gar::geo::LocalTransformation<TGeoHMatrix> trans;
                    fGeo->WorldToLocal(GlobalPosCM, LocalPosCM, trans);

                    LOG_DEBUG("ConvertEdep2Art")
                    << "Sensitive volume " << d->first
                    << " Hit " << hit
                    << " in volume " << VolumeName
                    << " in material " << volmaterial
                    << " det_id " << det_id
                    << " layer " << layer
                    << " slice " << slice;

                    gar::raw::CellID_t cellID = fGeo->GetCellID(node, det_id, 0, 0, layer, slice, LocalPosCM);//encoding the cellID on 64 bits

                    double G4Pos[3] = {0., 0., 0.}; // in cm
                    G4Pos[0] = GlobalPosCM[0];
                    G4Pos[1] = GlobalPosCM[1];
                    G4Pos[2] = GlobalPosCM[2];

                    gar::sdp::CaloDeposit calohit( trackID, time, edep, G4Pos, cellID );
                    if(m_TrackerDeposits.find(cellID) != m_TrackerDeposits.end())
                    m_TrackerDeposits[cellID].push_back(calohit);
                    else {
                        std::vector<gar::sdp::CaloDeposit> vechit;
                        vechit.push_back(calohit);
                        m_TrackerDeposits.emplace(cellID, vechit);
                    }
                }
            }
            else if( d->first == "MuID_vol" ) {
                //MuonID detector in the SPY
                for (std::vector<TG4HitSegment>::const_iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    const TG4HitSegment *hit = &(*h);

                    int trackID = hit->GetPrimaryId();
                    double edep = VisibleEnergyDeposition(hit, fApplyBirks) * CLHEP::MeV / CLHEP::GeV;
                    double time = (hit->GetStart().T() + hit->GetStop().T())/2 / CLHEP::s;
                    double x = (hit->GetStart().X() + hit->GetStop().X())/2 /CLHEP::cm;
                    double y = (hit->GetStart().Y() + hit->GetStop().Y())/2 /CLHEP::cm;
                    double z = (hit->GetStart().Z() + hit->GetStop().Z())/2 /CLHEP::cm;

                    if(edep == 0 || edep < fEnergyCut)
                    continue;

                    //Check if it is in the active material of the ECAL
                    TGeoNode *node = fGeo->FindNode(x, y, z);//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) continue;

                    unsigned int layer = GetLayerNumber(VolumeName); //get layer number
                    unsigned int slice = GetSliceNumber(VolumeName); // get slice number
                    unsigned int det_id = 4;
                    unsigned int stave = GetStaveNumber(VolumeName);
                    unsigned int module = GetModuleNumber(VolumeName);

                    std::array<double, 3> GlobalPosCM = {x, y, z};
                    std::array<double, 3> LocalPosCM;
                    gar::geo::LocalTransformation<TGeoHMatrix> trans;
                    fGeo->WorldToLocal(GlobalPosCM, LocalPosCM, trans);

                    LOG_DEBUG("ConvertEdep2Art")
                    << "Sensitive volume " << d->first
                    << " Hit " << hit
                    << " in volume " << VolumeName
                    << " in material " << volmaterial
                    << " det_id " << det_id
                    << " layer " << layer
                    << " slice " << slice
                    << " stave " << stave
                    << " module " << module;

                    gar::raw::CellID_t cellID = fGeo->GetCellID(node, det_id, stave, module, layer, slice, LocalPosCM);//encoding the cellID on 64 bits

                    double G4Pos[3] = {0., 0., 0.}; // in cm
                    G4Pos[0] = GlobalPosCM[0];
                    G4Pos[1] = GlobalPosCM[1];
                    G4Pos[2] = GlobalPosCM[2];

                    gar::sdp::CaloDeposit calohit( trackID, time, edep, G4Pos, cellID );
                    if(m_MuIDDeposits.find(cellID) != m_MuIDDeposits.end())
                    m_MuIDDeposits[cellID].push_back(calohit);
                    else {
                        std::vector<gar::sdp::CaloDeposit> vechit;
                        vechit.push_back(calohit);
                        m_MuIDDeposits.emplace(cellID, vechit);
                    }
                }
            }
            else{
                LOG_DEBUG("ConvertEdep2Art")
                << "Ignoring hits for sensitive material: "
                << d->first;
                continue;
            }
        }

        LOG_DEBUG("ConvertEdep2Art") << "Finished collection sensitive hits";

        //--------------------------------------------------------------------------

        bool hasGAr, hasECAL, hasTrackerSc, hasMuID, hasLAr = false;
        if(fGArDeposits.size() > 0) hasGAr = true;
        if(m_ECALDeposits.size() > 0) hasECAL = true;
        if(m_TrackerDeposits.size() > 0) hasTrackerSc = true;
        if(m_MuIDDeposits.size() > 0) hasMuID = true;

        if(hasGAr) {
            std::sort(fGArDeposits.begin(), fGArDeposits.end());

            for(auto const& garhit : fGArDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding GAr deposits for track id: "
                << garhit.TrackID();
                TPCCol->emplace_back(garhit);
            }
        }

        if(hasECAL) {
            this->AddHits(m_ECALDeposits, fECALDeposits);
            std::sort(fECALDeposits.begin(), fECALDeposits.end());

            for(auto const& ecalhit : fECALDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding calo deposits for track id: "
                << ecalhit.TrackID();

                ECALCol->emplace_back(ecalhit);
            }
        }

        if(hasTrackerSc) {
            this->AddHitsMinerva(m_TrackerDeposits, fTrackerDeposits);
            std::sort(fTrackerDeposits.begin(), fTrackerDeposits.end());

            for(auto const& trkhit : fTrackerDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding tracker Sc deposits for track id: "
                << trkhit.TrackID();

                TrackerCol->emplace_back(trkhit);
            }
        }

        if(hasMuID) {
            this->AddHits(m_MuIDDeposits, fMuIDDeposits);
            std::sort(fMuIDDeposits.begin(), fMuIDDeposits.end());

            for(auto const& muidhit : fMuIDDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding muID deposits for track id: "
                << muidhit.TrackID();

                MuIDCol->emplace_back(muidhit);
            }
        }

        //Create assn between hits and mcp
        art::PtrMaker<simb::MCParticle> makeMCPPtr(evt);
        art::PtrMaker<gar::sdp::EnergyDeposit> makeEnergyDepositPtr(evt);
        art::PtrMaker<gar::sdp::CaloDeposit> makeCaloDepositPtr(evt, "ECAL");
        art::PtrMaker<gar::sdp::CaloDeposit> makeTrackerDepositPtr(evt, "TrackerSc");
        art::PtrMaker<gar::sdp::CaloDeposit> makeMuIDDepositPtr(evt, "MuID");

        unsigned int imcp = 0;
        for(auto const &part : *partCol)
        {
            int mpc_trkid = part.TrackId();
            art::Ptr<simb::MCParticle> partPtr = makeMCPPtr(imcp);

            unsigned int igashit = 0;
            unsigned int iecalhit = 0;
            unsigned int itrkhit = 0;
            unsigned int imuidhit = 0;

            if(hasGAr) {
                for(auto const& gashit : *TPCCol)
                {
                    if(mpc_trkid == gashit.TrackID()){
                        art::Ptr<gar::sdp::EnergyDeposit> gashitPtr = makeEnergyDepositPtr(igashit);
                        ghmcassn->addSingle(gashitPtr, partPtr);
                    }
                    igashit++;
                }
            }

            if(hasECAL) {
                for(auto const& ecalhit : *ECALCol)
                {
                    if(mpc_trkid == ecalhit.TrackID()){
                        art::Ptr<gar::sdp::CaloDeposit> ecalhitPtr = makeCaloDepositPtr(iecalhit);
                        ehmcassn->addSingle(ecalhitPtr, partPtr);
                    }
                    iecalhit++;
                }
            }

            if(hasTrackerSc) {
                for(auto const& trkhit : *TrackerCol)
                {
                    if(mpc_trkid == trkhit.TrackID()){
                        art::Ptr<gar::sdp::CaloDeposit> trkhitPtr = makeTrackerDepositPtr(itrkhit);
                        thmcassn->addSingle(trkhitPtr, partPtr);
                    }
                    itrkhit++;
                }
            }

            if(hasMuID) {
                for(auto const& muidhit : *MuIDCol)
                {
                    if(mpc_trkid == muidhit.TrackID()){
                        art::Ptr<gar::sdp::CaloDeposit> muIDhitPtr = makeMuIDDepositPtr(imuidhit);
                        mhmcassn->addSingle(muIDhitPtr, partPtr);
                    }
                    imuidhit++;
                }
            }

            imcp++;
        }

        evt.put(std::move(mctruthcol));
        if(fHasGHEP) {
            evt.put(std::move(gtruthcol));
            evt.put(std::move(tgassn));
        }
        evt.put(std::move(tpassn));
        evt.put(std::move(partCol));

        if(hasGAr) {
            evt.put(std::move(TPCCol));
            evt.put(std::move(ghmcassn));
        }
        if(hasECAL) {
            evt.put(std::move(ECALCol), "ECAL");
            evt.put(std::move(ehmcassn), "ECAL");
        }
        if(hasTrackerSc) {
            evt.put(std::move(TrackerCol), "TrackerSc");
            evt.put(std::move(thmcassn), "TrackerSc");
        }
        if(hasMuID) {
            evt.put(std::move(MuIDCol), "MuID");
            evt.put(std::move(mhmcassn), "MuID");
        }
        if(hasLAr) {
            // evt.put(std::move(LArCol));
            // evt.put(std::move(lhmcassn));
        }

        return;
    }
} // namespace util

namespace util {

    DEFINE_ART_MODULE(ConvertEdep2Art)

} // namespace util

#endif
