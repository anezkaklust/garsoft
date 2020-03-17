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

///Geant4 interface
namespace gar {
    namespace garg4 {

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

            std::string fInputfile;
            std::string fECALMaterial;
            bool fAddMCFlux;
            bool fApplyBirks;

            TChain* fTreeChain;
            unsigned int nEntries;
            TG4Event* event;
            double fBirksCoeff;

            const gar::geo::GeometryCore*            fGeo;               ///< geometry information
            TGeoManager*                             fGeoManager;
            const detinfo::ECALProperties*           fEcalProp;

            std::vector<gar::sdp::CaloDeposit> fECALDeposits;
            std::vector<gar::sdp::EnergyDeposit> fGArDeposits;
        };

    } // namespace garg4

    namespace garg4 {

        //----------------------------------------------------------------------
        // Constructor
        ConvertEdep2Art::ConvertEdep2Art(fhicl::ParameterSet const& pset)
        : art::EDProducer{pset},
        fInputfile( pset.get< std::string >("Inputfile", "") ),
        fECALMaterial( pset.get< std::string >("ECALMaterial", "Scintillator") ),
        fAddMCFlux( pset.get< bool >("addMCFlux", false) ),
        fApplyBirks( pset.get< bool >("applyBirks", true) ),
        fTreeChain(new TChain("EDepSimEvents"))
        {
            // produces< std::vector<simb::MCTruth> >();
            // produces< std::vector<simb::GTruth>  >();
            // if ( fAddMCFlux ) {
            //     produces< std::vector<simb::MCFlux>  >();
            // }
            //
            // produces< std::vector<simb::MCParticle> >();
            produces< std::vector<sdp::EnergyDeposit> >();
            produces< std::vector<sdp::CaloDeposit> >();
            // produces< std::vector<sdp::LArDeposit> >();

            fTreeChain->Add(fInputfile.c_str());
            nEntries = fTreeChain->GetEntries();

            event = nullptr;
            fTreeChain->SetBranchAddress("Event", &event);

            fGeo = gar::providerFrom<geo::Geometry>();
            fGeoManager = fGeo->ROOTGeoManager();
            fEcalProp = gar::providerFrom<detinfo::ECALPropertiesService>();
            fBirksCoeff = fEcalProp->ScintBirksConstant(); //CLHEP::mm/CLHEP::MeV;
        }

        //----------------------------------------------------------------------
        // Destructor
        ConvertEdep2Art::~ConvertEdep2Art()
        {

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
            double edep = hit->EnergyDeposit;
            double niel = hit->SecondaryDeposit;
            double length = hit->TrackLength;
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

            // std::unique_ptr< std::vector<simb::MCParticle> > partCol(new std::vector<simb::MCParticle>                );
            // std::unique_ptr< ::art::Assns<simb::MCTruth, simb::MCParticle> > tpassn (new ::art::Assns<simb::MCTruth, simb::MCParticle>);

            // std::unique_ptr< std::vector< sdp::LArDeposit > >           LArCol  (new std::vector<sdp::LArDeposit>           );
            std::unique_ptr< std::vector<sdp::EnergyDeposit>  > TPCCol(new std::vector<sdp::EnergyDeposit> );
            std::unique_ptr< std::vector< sdp::CaloDeposit > > ECALCol(new std::vector<sdp::CaloDeposit> );

            fGArDeposits.clear();
            fECALDeposits.clear();

            if( eventnumber > nEntries ){
                //stop...
                return;
            }

            //Get the event
            fTreeChain->GetEntry(eventnumber);

            //Fill simulated hits
            for (auto d = event->SegmentDetectors.begin(); d != event->SegmentDetectors.end(); ++d)
            {
                if( d->first == "TPC_Drift1" || d->first == "TPC_Drift2" )
                {
                    //GAr deposits
                    for (std::vector<TG4HitSegment>::iterator h = d->second.begin(); h != d->second.end(); ++h)
                    {
                        auto trackID = h->PrimaryId;
                        double edep = h->EnergyDeposit / CLHEP::MeV;
                        double time = (h->Start.T() + h->Stop.T())/2 / CLHEP::ns;
                        double x = (h->Start.X() + h->Stop.X())/2;
                        double y = (h->Start.Y() + h->Stop.Y())/2;
                        double z = (h->Start.Z() + h->Stop.Z())/2;
                        double stepLength = h->TrackLength;
                        fGArDeposits.emplace_back(trackID, time, edep, x/CLHEP::cm, y/CLHEP::cm, z/CLHEP::cm, stepLength/CLHEP::cm, (trackID > 0));
                    }
                }
                else if( d->first == "BarrelECal_vol" || d->first == "EndcapECal_vol"){
                    //ECAL deposits
                    for (std::vector<TG4HitSegment>::iterator h = d->second.begin(); h != d->second.end(); ++h)
                    {
                        auto trackID = h->PrimaryId;
                        double edep = VisibleEnergyDeposition(h, fApplyBirks);
                        double time = (h->Start.T() + h->Stop.T())/2 / CLHEP::ns;
                        TVector3 GlobalPosCM( (h->Start.X() + h->Stop.X())/2 /CLHEP::cm,
                        (h->Start.Y() + h->Stop.Y())/2 /CLHEP::cm ,
                        (h->Start.Z() + h->Stop.Z())/2 /CLHEP::cm );

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

            for(auto hit : fGArDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding GAr deposits for track id: "
                << hit.TrackID();
                TPCCol->emplace_back(hit);
            }
            for(auto hit : fECALDeposits)
            {
                LOG_DEBUG("ConvertEdep2Art")
                << "adding calo deposits for track id: "
                << hit.TrackID();
                ECALCol->emplace_back(hit);
            }
            evt.put(std::move(TPCCol));
            evt.put(std::move(ECALCol));
            // evt.put(std::move(LArCol));
            // evt.put(std::move(partCol));
            // evt.put(std::move(tpassn));
            return;
        }
    } // namespace garg4
    namespace garg4 {
        DEFINE_ART_MODULE(ConvertEdep2Art)
    } // namespace garg4
} // gar

#endif
