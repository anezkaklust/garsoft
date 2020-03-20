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

        unsigned int GetDetNumber(std::string volname) const;
        unsigned int GetStaveNumber(std::string volname) const;
        unsigned int GetModuleNumber(std::string volname) const;
        unsigned int GetLayerNumber(std::string volname) const;
        unsigned int GetSliceNumber(std::string volname) const;

        double VisibleEnergyDeposition(const TG4HitSegment *hit, bool applyBirks) const;
        bool CheckProcess( std::string process_name ) const;

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
        std::string fEMShowerDaughterMatRegex;   ///< keep EM shower daughters only in these materials

        genie::PDGLibrary *pdglib;
        const gar::geo::GeometryCore*            fGeo;               ///< geometry information
        TGeoManager*                             fGeoManager;
        const gar::detinfo::ECALProperties*      fEcalProp;

        std::vector<simb::MCParticle> fMCParticles;
        std::vector<gar::sdp::CaloDeposit> fECALDeposits;
        std::vector<gar::sdp::EnergyDeposit> fGArDeposits;

        int fSpillCount;
        std::vector<int> fStartSpill;
        std::vector<int> fStopSpill;
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
    fEMShowerDaughterMatRegex( pset.get< std::string >("EMShowerDaughterMatRegex", ".*") )
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
        fBirksCoeff = fEcalProp->ScintBirksConstant(); //CLHEP::mm/CLHEP::MeV;

        fSpillCount = 0;
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
        if(fOverlay){
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

            mf::LogInfo("ConvertEdep2Art::beginJob()")
            << "Number of spills in the ghep file "
            << fSpillCount;
        }

        return;
    }

    //--------------------------------------------------------------------------
    void ConvertEdep2Art::beginRun(::art::Run& run)
    {

        return;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetLayerNumber(std::string volname) const
    {
        return std::atoi( (volname.substr( volname.find("layer_") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetSliceNumber(std::string volname) const
    {
        return std::atoi( (volname.substr( volname.find("slice") + 5, 1)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetDetNumber(std::string volname) const
    {
        unsigned int det_id = 0;
        if( volname.find("Barrel") !=  std::string::npos )
        det_id = 1;
        if( volname.find("Endcap") !=  std::string::npos )
        det_id = 2;
        return det_id;
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetStaveNumber(std::string volname) const
    {
        return std::atoi( (volname.substr( volname.find("_stave") + 6, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    unsigned int ConvertEdep2Art::GetModuleNumber(std::string volname) const
    {
        return std::atoi( (volname.substr( volname.find("_module") + 7, 2)).c_str() );
    }

    //------------------------------------------------------------------------------
    double ConvertEdep2Art::VisibleEnergyDeposition(const TG4HitSegment *hit, bool applyBirks) const
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
    bool ConvertEdep2Art::CheckProcess( std::string process_name ) const
    {
        bool isEMShowerProcess = false;

        if( process_name.find("conv")            != std::string::npos  ||
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

        //--------------------------------------------------------------------------
        if(fOverlay){
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

                evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
            }
        }
        else{
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

            evgb::util::CreateAssn(*this, evt, *mctruthcol, *gtruthcol, *tgassn, gtruthcol->size()-1, gtruthcol->size());
        }

        //--------------------------------------------------------------------------
        //Get the event
        //Starts at 0, evt starts at 1
        fTreeChain->GetEntry(eventnumber-1);

        std::unique_ptr< std::vector<simb::MCParticle> > partCol(new std::vector<simb::MCParticle> );

        fMCParticles.clear();

        for (std::vector<TG4Trajectory>::const_iterator t = fEvent->Trajectories.begin(); t != fEvent->Trajectories.end(); ++t)
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
                process_name = t->Points.at(0).GetProcessName();
            }

            LOG_DEBUG("ConvertEdep2Art")
            << "Particle " << name
            << " with pdg " << pdg
            << " trackID " << trkID
            << " parent id " << parentID
            << " created with process [ " << process_name << " ]";

            //Skip the creation of the mcp if it is part of shower (based on process name) and only if it is not matching the material
            //TODO debug as it does not work for anatree at the moment, problem related to pdg mom exception
            //Maybe need to keep relation trackid, parentid

            bool isEMShowerProcess = CheckProcess( process_name );
            TLorentzVector part_start = t->Points.at(0).GetPosition();
            TGeoNode *node = fGeoManager->FindNode(part_start.X() / CLHEP::cm, part_start.Y() / CLHEP::cm, part_start.Z() / CLHEP::cm);
            std::string material_name = "";

            if(node)
            material_name = node->GetMedium()->GetMaterial()->GetName();
            // std::cout << "Material name " << material_name << std::endl;

            std::regex const re_material(fEMShowerDaughterMatRegex);
            if( isEMShowerProcess && !std::regex_match(material_name, re_material) ) {

                LOG_DEBUG("ConvertEdep2Art")
                << "Skipping EM shower daughter " << name
                << " with trackID " << trkID
                << " with parent id " << parentID
                << " created with process [ " << process_name << " ]"
                << " in material " << material_name;

                continue;
            }

            simb::MCParticle mcpart(trkID, pdg, process_name, parentID, mass);

            for (std::vector<TG4TrajectoryPoint>::const_iterator p = t->Points.begin(); p != t->Points.end(); ++p)
            {
                TLorentzVector position = p->GetPosition();
                TLorentzVector fourPos(position.X() / CLHEP::cm, position.Y() / CLHEP::cm, position.Z() / CLHEP::cm, position.T() / CLHEP::s );
                TVector3 momentum = p->GetMomentum();
                TLorentzVector fourMom(momentum.x() * CLHEP::MeV / CLHEP::GeV, momentum.y() * CLHEP::MeV / CLHEP::GeV, momentum.z() * CLHEP::MeV / CLHEP::GeV, 0.);
                std::string process = p->GetProcessName();

                if(p == t->Points.begin()) process = "Start";
                mcpart.AddTrajectoryPoint(fourPos, fourMom, process);
            }

            std::string end_process = t->Points.at(t->Points.size()-1).GetProcessName();
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

                    TGeoNode *node = fGeoManager->FindNode(x, y, z);//Node in cm...
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

                    TVector3 GlobalPosCM(x, y, z);
                    //Check if it is in the active material of the ECAL
                    TGeoNode *node = fGeoManager->FindNode(GlobalPosCM.x(), GlobalPosCM.y(), GlobalPosCM.z());//Node in cm...
                    std::string VolumeName  = node->GetVolume()->GetName();
                    std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
                    if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) continue;

                    unsigned int layer = GetLayerNumber(VolumeName); //get layer number
                    unsigned int slice = GetSliceNumber(VolumeName); // get slice number
                    unsigned int det_id = GetDetNumber(VolumeName); // 1 == Barrel, 2 = Endcap
                    unsigned int stave = GetStaveNumber(VolumeName); //get the stave number
                    unsigned int module = GetModuleNumber(VolumeName); //get the module number

                    TVector3 LocalPosCM;
                    fGeo->WorldToLocal(GlobalPosCM, LocalPosCM);
                    G4ThreeVector G4LocalPosCM(LocalPosCM.x(), LocalPosCM.y(), LocalPosCM.z());

                    LOG_DEBUG("ConvertEdep2Art")
                    << "Hit " << hit
                    << " in volume " << VolumeName
                    << " in material " << volmaterial
                    << " det_id " << det_id
                    << " module " << module
                    << " stave " << stave
                    << " layer " << layer
                    << " slice " << slice;

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
                LOG_DEBUG("ConvertEdep2Art")
                << "Ignoring hits for sensitive material: "
                << d->first;
                continue;
            }
        }

        //--------------------------------------------------------------------------

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
