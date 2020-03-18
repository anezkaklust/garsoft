#ifndef CONVEDEP2ART_H
#define CONVEDEP2ART_H

// C++ Includes
#include <sstream>
#include <vector>
#include <map>
#include <set>
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
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
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
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/Geometry.h"
#include "Geometry/LocalTransformation.h"
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

        unsigned int GetDetNumber(std::string volname);
        unsigned int GetStaveNumber(std::string volname);
        unsigned int GetModuleNumber(std::string volname);
        unsigned int GetLayerNumber(std::string volname);
        unsigned int GetSliceNumber(std::string volname);

        double VisibleEnergyDeposition(std::vector<TG4HitSegment>::iterator hit, bool applyBirks);

        std::string fEDepSimfile;
        std::string fGhepfile;
        bool fOverlay;
        std::string fECALMaterial;
        bool fApplyBirks;

        TChain* fTreeChain;
        unsigned int nEntries;

        TChain* fGTreeChain;
        unsigned int nEntriesGhep;
        genie::NtpMCEventRecord*                 fMCRec;

        TG4Event* fEvent;
        double fBirksCoeff;

        genie::PDGLibrary *pdglib;
        const gar::geo::GeometryCore*            fGeo;               ///< geometry information
        TGeoManager*                             fGeoManager;
        const gar::detinfo::ECALProperties*      fEcalProp;

        std::vector<simb::MCParticle> fMCParticles;
        std::vector<gar::sdp::CaloDeposit> fECALDeposits;
        std::vector<gar::sdp::EnergyDeposit> fGArDeposits;
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
    fApplyBirks( pset.get< bool >("applyBirks", true) ),
    fTreeChain(new TChain("EDepSimEvents")),
    fGTreeChain(new TChain("gtree")),
    fMCRec(nullptr),
    fEvent(nullptr)
    {
        pdglib = genie::PDGLibrary::Instance();

        produces< std::vector<simb::MCTruth> >();
        produces< std::vector<simb::GTruth>  >();
        produces< art::Assns<simb::MCTruth, simb::GTruth> >();

        produces< std::vector<simb::MCParticle> >();
        produces< std::vector<gar::sdp::EnergyDeposit> >();
        produces< std::vector<gar::sdp::CaloDeposit> >();
        // produces< std::vector<sdp::LArDeposit> >();

        if(fEDepSimfile.empty() || fGhepfile.empty())
        {
            throw cet::exception("ConvertEdep2Art")
            << "Empty edep-sim file or ghep file";
        }

        fTreeChain->Add(fEDepSimfile.c_str());
        nEntries = fTreeChain->GetEntries();
        fTreeChain->SetBranchAddress("Event", &fEvent);

        fGTreeChain->Add(fGhepfile.c_str());
        nEntriesGhep = fTreeChain->GetEntries();
        fGTreeChain->SetBranchAddress("gmcrec", &fMCRec);

        fGeo = gar::providerFrom<gar::geo::Geometry>();
        fGeoManager = fGeo->ROOTGeoManager();
        fEcalProp = gar::providerFrom<gar::detinfo::ECALPropertiesService>();
        fBirksCoeff = fEcalProp->ScintBirksConstant(); //CLHEP::mm/CLHEP::MeV;;
    }

    //----------------------------------------------------------------------
    // Destructor
    ConvertEdep2Art::~ConvertEdep2Art()
    {
        if ( fTreeChain ) delete fTreeChain;
        if ( fGTreeChain ) delete fGTreeChain;
    }

    //----------------------------------------------------------------------
    void ConvertEdep2Art::beginJob()
    {
        return;
    }

    //--------------------------------------------------------------------------
    void ConvertEdep2Art::beginRun(::art::Run& run)
    {

        return;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetLayerNumber(std::string volname)
    {
        return std::atoi( (volname.substr( volname.find("layer_") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetSliceNumber(std::string volname)
    {
        return std::atoi( (volname.substr( volname.find("slice") + 5, 1)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetDetNumber(std::string volname)
    {
        unsigned int det_id = 0;
        if( volname.find("Barrel") !=  std::string::npos )
        det_id = 1;
        if( volname.find("Endcap") !=  std::string::npos )
        det_id = 2;
        return det_id;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetStaveNumber(std::string volname)
    {
        return std::atoi( (volname.substr( volname.find("_stave") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetModuleNumber(std::string volname)
    {
        return std::atoi( (volname.substr( volname.find("_module") + 7, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    double ConvertEdep2Art::VisibleEnergyDeposition(std::vector<TG4HitSegment>::iterator hit, bool applyBirks)
    {
        double edep = hit->GetEnergyDeposit();
        double niel = hit->GetSecondaryDeposit();
        double length = hit->GetTrackLength();
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

        return evis;
    }

    //--------------------------------------------------------------------------
    void ConvertEdep2Art::produce(::art::Event& evt)
    {
        LOG_DEBUG("ConvertEdep2Art") << "produce()";
        art::EventNumber_t eventnumber = evt.id().event();

        if( eventnumber > nEntries ){
            //stop...
            return;
        }

        std::unique_ptr< std::vector<simb::MCTruth> > mctruthcol(new std::vector<simb::MCTruth>);
        std::unique_ptr< std::vector<simb::GTruth> > gtruthcol(new std::vector<simb::GTruth>);
        std::unique_ptr< art::Assns<simb::MCTruth, simb::GTruth> > tgassn(new art::Assns<simb::MCTruth, simb::GTruth>);

        if(fOverlay){
            int fSpillCount = 0;
            int fStartSpill[10000];
            int fStopSpill[10000];

            //need to get where is the end of the spill
            for(int ientry = 0; ientry < fGTreeChain->GetEntries(); ientry++)
            {
                fGTreeChain->GetEntry(ientry);

                genie::NtpMCRecHeader rec_header = fMCRec->hdr;
                genie::EventRecord *event = fMCRec->event;
                // genie::Interaction *interaction = event->Summary();

                // std::cout << rec_header;
                // std::cout << *event;

                //Need to get the Rootino showing the end of spill
                genie::GHepParticle *neutrino = event->Probe();
                if(nullptr == neutrino){
                    if(fSpillCount == 0)
                    fStartSpill[fSpillCount] = 0;
                    else{
                        fStartSpill[fSpillCount] = fStopSpill[fSpillCount-1]+1;
                    }

                    fStopSpill[fSpillCount] = ientry;
                    fSpillCount++;
                }
            }

            mf::LogInfo("ConvertEdep2Art") << "Number of spills in the ghep file " << fSpillCount << std::endl;

            for(int ientry = fStartSpill[eventnumber]; ientry < fStopSpill[eventnumber]; ientry++)
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

                evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
            }
        }
        else{
            fGTreeChain->GetEntry(eventnumber);

            genie::NtpMCRecHeader rec_header = fMCRec->hdr;
            genie::EventRecord *event = fMCRec->event;

            simb::MCTruth mctruth;
            simb::GTruth  gtruth;

            evgb::FillMCTruth(event, 0., mctruth);
            evgb::FillGTruth(event, gtruth);

            mctruthcol->push_back(mctruth);
            gtruthcol->push_back(gtruth);

            evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
        }

        //Get the event
        fTreeChain->GetEntry(eventnumber);

        //--------------------------------------------------------------------------
        std::unique_ptr< std::vector<simb::MCParticle> > partCol(new std::vector<simb::MCParticle> );

        fMCParticles.clear();

        for (std::vector<TG4Trajectory>::iterator t = fEvent->Trajectories.begin(); t != fEvent->Trajectories.end(); ++t)
        {
            int trkID = t->GetTrackId();
            int parentID = t->GetParentId();
            int pdg = t->GetPDGCode();
            std::string name = t->GetName();
            double mass = pdglib->Find(pdg)->Mass();

            std::string process_name = "unknown";
            if(parentID == -1) { process_name = "primary"; parentID = 0; }
            else{
                //Get the first point that created this particle
                process_name = std::to_string( t->Points.at(0).GetProcess() );
            }

            // std::cout << "Part name " << name << " pdg " << pdg << " trkID " << trkID << " parentID " << parentID << " process " << process_name << std::endl;

            simb::MCParticle mcpart(trkID, pdg, process_name, parentID, mass);

            for (std::vector<TG4TrajectoryPoint>::iterator p = t->Points.begin(); p != t->Points.end(); ++p)
            {
                TLorentzVector position = p->GetPosition();
                TLorentzVector fourPos(position.X() / CLHEP::cm, position.Y() / CLHEP::cm, position.Z() / CLHEP::cm, position.T() / CLHEP::s );
                TVector3 momentum = p->GetMomentum();
                TLorentzVector fourMom(momentum.x() * CLHEP::MeV / CLHEP::GeV, momentum.y() * CLHEP::MeV / CLHEP::GeV, momentum.z() * CLHEP::MeV / CLHEP::GeV, 0.);
                std::string process = std::to_string( p->GetProcess() + 1000*p->GetSubprocess() );

                if(p == t->Points.begin()) process = "Start";
                mcpart.AddTrajectoryPoint(fourPos, fourMom, process);
            }

            std::string end_process = std::to_string( t->Points.at(t->Points.size()-1).GetProcess() );
            mcpart.SetEndProcess(end_process);

            fMCParticles.push_back( mcpart );
        }

        for(auto const& mcpart : fMCParticles)
        {
            LOG_DEBUG("ConvertEdep2Art")
            << "adding mc particle with track id: "
            << mcpart.TrackId();
            partCol->emplace_back(mcpart);
        }

        //--------------------------------------------------------------------------
        std::unique_ptr< std::vector< gar::sdp::EnergyDeposit>  > TPCCol(new std::vector<gar::sdp::EnergyDeposit> );
        std::unique_ptr< std::vector< gar::sdp::CaloDeposit > > ECALCol(new std::vector<gar::sdp::CaloDeposit> );

        fGArDeposits.clear();
        fECALDeposits.clear();

        //Fill simulated hits
        for (auto d = fEvent->SegmentDetectors.begin(); d != fEvent->SegmentDetectors.end(); ++d)
        {
            if( d->first == "TPC_Drift1" || d->first == "TPC_Drift2" )
            {
                //GAr deposits
                for (std::vector<TG4HitSegment>::iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    auto trackID = h->GetPrimaryId();
                    double edep = h->GetEnergyDeposit() / CLHEP::MeV;
                    double time = (h->GetStart().T() + h->GetStop().T())/2 / CLHEP::s;
                    double x = (h->GetStart().X() + h->GetStop().X())/2;
                    double y = (h->GetStart().Y() + h->GetStop().Y())/2;
                    double z = (h->GetStart().Z() + h->GetStop().Z())/2;
                    double stepLength = h->GetTrackLength();
                    fGArDeposits.emplace_back(trackID, time, edep, x/CLHEP::cm, y/CLHEP::cm, z/CLHEP::cm, stepLength/CLHEP::cm, (trackID > 0));
                }
            }
            else if( d->first == "BarrelECal_vol" || d->first == "EndcapECal_vol"){
                //ECAL deposits
                for (std::vector<TG4HitSegment>::iterator h = d->second.begin(); h != d->second.end(); ++h)
                {
                    auto trackID = h->GetPrimaryId();
                    double edep = VisibleEnergyDeposition(h, fApplyBirks);
                    double time = (h->GetStart().T() + h->GetStop().T())/2 / CLHEP::s;
                    TVector3 GlobalPosCM( (h->GetStart().X() + h->GetStop().X())/2 /CLHEP::cm,
                    (h->GetStart().Y() + h->GetStop().Y())/2 /CLHEP::cm ,
                    (h->GetStart().Z() + h->GetStop().Z())/2 /CLHEP::cm );

                    TGeoNode *node = fGeoManager->FindNode(GlobalPosCM.x(), GlobalPosCM.y(), GlobalPosCM.z());//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) continue;

                    // std::cout << "Volume " << VolumeName << " material " << volmaterial << std::endl;
                    unsigned int layer = GetLayerNumber(VolumeName); //get layer number
                    unsigned int slice = GetSliceNumber(VolumeName); // get slice number
                    unsigned int det_id = GetDetNumber(VolumeName); // 1 == Barrel, 2 = Endcap
                    unsigned int stave = GetStaveNumber(VolumeName); //get the stave number
                    unsigned int module = GetModuleNumber(VolumeName); //get the module number

                    TVector3 LocalPosCM;
                    fGeo->WorldToLocal(GlobalPosCM, LocalPosCM);
                    G4ThreeVector G4LocalPosCM(LocalPosCM.x(), LocalPosCM.y(), LocalPosCM.z());

                    // std::cout << "layer " << layer;
                    // std::cout << " slice " << slice;
                    // std::cout << " det_id " << det_id;
                    // std::cout << " stave " << stave;
                    // std::cout << " module " << module;
                    // std::cout << std::endl;

                    long long int cellID = fGeo->cellID(node, det_id, stave, module, layer, slice, G4LocalPosCM);//encoding the cellID on 64 bits

                    double G4Pos[3] = {0., 0., 0.}; // in cm
                    if(fGeo->isTile(cellID))
                    {
                        G4ThreeVector G4SegLocalCM = fGeo->position(node, cellID);//in cm
                        TVector3 SegLocalCM(G4SegLocalCM.x(), G4SegLocalCM.y(), G4SegLocalCM.z());
                        TVector3 SegGlobalCM;
                        fGeo->LocalToWorld(SegLocalCM, SegGlobalCM);

                        G4Pos[0] = SegGlobalCM.x();
                        G4Pos[1] = SegGlobalCM.y();
                        G4Pos[2] = SegGlobalCM.z();
                    } else {
                        G4Pos[0] = GlobalPosCM.x();
                        G4Pos[1] = GlobalPosCM.y();
                        G4Pos[2] = GlobalPosCM.z();
                    }

                    fECALDeposits.emplace_back(trackID, time, edep, G4Pos, cellID);
                }
            }
            else{
                continue;
            }
        }

        std::sort(fGArDeposits.begin(), fGArDeposits.end());
        std::sort(fECALDeposits.begin(), fECALDeposits.end());

        for(auto const& garhit : fGArDeposits)
        {
            LOG_DEBUG("ConvertEdep2Art")
            << "adding GAr deposits for track id: "
            << garhit.TrackID();
            TPCCol->emplace_back(garhit);
        }

        for(auto const& ecalhit : fECALDeposits)
        {
            LOG_DEBUG("ConvertEdep2Art")
            << "adding calo deposits for track id: "
            << ecalhit.TrackID();
            ECALCol->emplace_back(ecalhit);
        }

        evt.put(std::move(mctruthcol));
        evt.put(std::move(gtruthcol));
        evt.put(std::move(tgassn));

        evt.put(std::move(partCol));
        evt.put(std::move(TPCCol));
        evt.put(std::move(ECALCol));
        // evt.put(std::move(LArCol));
        // evt.put(std::move(tpassn));
        return;
    }
} // namespace util

namespace util {

    DEFINE_ART_MODULE(ConvertEdep2Art)

} // namespace util

#endif
