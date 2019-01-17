//
//  AuxDetAction.cxx
//
//  Created by Brian Rebel on 11/22/16.
//

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "Math/Vector3D.h"

// G4 includes
#include "Geant4/G4Event.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PrimaryParticle.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4VUserPrimaryParticleInformation.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4VProcess.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "Geant4/G4NavigationHistory.hh"

#include "GArG4/AuxDetAction.h"
#include "GArG4/ParticleListAction.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

namespace gar {

    namespace garg4 {

        //-------------------------------------------------------------
        // Constructor.
        AuxDetAction::AuxDetAction(CLHEP::HepRandomEngine*    engine,
        fhicl::ParameterSet const& pset)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            fGeoManager = fGeo->ROOTGeoManager();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            this->reconfigure(pset);
        }

        //-------------------------------------------------------------
        // Destructor.
        AuxDetAction::~AuxDetAction()
        {
            // Delete anything that we created with "new'.
        }

        //-------------------------------------------------------------
        void AuxDetAction::reconfigure(fhicl::ParameterSet const& pset)
        {
            fECALMaterial = pset.get<std::string>("ECALMaterial");

            fLArEnergyCut = pset.get<double>("LArEnergyCut");
            fLArVolumeName = pset.get<std::vector<std::string>>("LArVolumeName");
            fLArMaterial = pset.get<std::string>("LArMaterial");

            std::cout << "AuxDetAction: Name of the Volumes to track for the LArTPC" << std::endl;
            for(unsigned int i = 0; i < fLArVolumeName.size(); i++) std::cout << fLArVolumeName.at(i) << " ";
            std::cout << std::endl;

            return;
        }

        //-------------------------------------------------------------
        void AuxDetAction::BeginOfEventAction(const G4Event*)
        {
            // Clear any previous information.
            fECALDeposits.clear();
            fLArDeposits.clear();
        }

        //-------------------------------------------------------------
        void AuxDetAction::PreTrackingAction(const G4Track* track)
        {

        }

        //-------------------------------------------------------------
        void AuxDetAction::PostTrackingAction( const G4Track* /*track*/)
        {
        }

        //-------------------------------------------------------------
        // With every step, check and see if we are in an auxdet and save the energy
        void AuxDetAction::SteppingAction(const G4Step* step)
        {
            LOG_DEBUG("AuxDetAction")
            << "AuxDetAction::SteppingAction";

            this->LArSteppingAction(step);

            this->ECALSteppingAction(step);

        }// end of AuxDetAction::SteppingAction

        //------------------------------------------------------------------------------
        void AuxDetAction::EndOfEventAction(const G4Event*)
        {
            //sort per time the hits
            std::sort(fECALDeposits.begin(), fECALDeposits.end());
            std::sort(fLArDeposits.begin(), fLArDeposits.end());
        }

        //------------------------------------------------------------------------------
        void AuxDetAction::LArSteppingAction(const G4Step* step)
        {
            LOG_DEBUG("AuxDetAction")
            << "AuxDetAction::LArSteppingAction";

            // Get the pointer to the track
            G4Track *track = step->GetTrack();

            const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
            const CLHEP::Hep3Vector &stop  = track->GetPosition();

            // If it's a null step, don't use it.
            if(start == stop) return;
            if(step->GetTotalEnergyDeposit() == 0) return;

            // check that we are in the correct material to record a hit
            std::string VolumeName   = this->GetVolumeName(track);

            if( std::find( fLArVolumeName.begin(), fLArVolumeName.end(), VolumeName ) == fLArVolumeName.end() )
            return;

            // check the material
            auto pos = 0.5 * (start + stop);
            TGeoNode *node = fGeoManager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm...

            if(!node){
                LOG_DEBUG("AuxDetAction::LArSteppingAction")
                << "Node not found in "
                << pos.x() << " mm "
                << pos.y() << " mm "
                << pos.z() << " mm";
                return;
            }

            std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
            if ( ! std::regex_match(volmaterial, std::regex(fLArMaterial)) ) return;

            // only worry about energy depositions larger than the minimum required
            if(step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fLArEnergyCut){

                LOG_DEBUG("AuxDetAction::LArSteppingAction")
                << "In volume "
                << VolumeName
                << " Material is "
                << volmaterial
                << " step size is "
                << step->GetStepLength() / CLHEP::cm
                << " cm and deposited "
                << step->GetTotalEnergyDeposit()
                << " MeV of energy with a minimum of "
                << fLArEnergyCut
                << " required.";

                // the step mid point is used for the position of the deposit
                auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
                step->GetPostStepPoint()->GetPosition() );

                // get the track id for this step
                auto trackID  = ParticleListAction::GetCurrentTrackID();
                float time = step->GetPreStepPoint()->GetGlobalTime();

                fLArDeposits.emplace_back(trackID,
                time,//get the time of the first subhit
                step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV,
                midPoint.x() / CLHEP::cm,
                midPoint.y() / CLHEP::cm,
                midPoint.z() / CLHEP::cm,
                step->GetStepLength() / CLHEP::cm,
                (trackID > 0));

            } // end if enough energy to worry about this step
        }

        //------------------------------------------------------------------------------
        void AuxDetAction::ECALSteppingAction(const G4Step* step)
        {
            LOG_DEBUG("AuxDetAction")
            << "AuxDetAction::ECALSteppingAction";

            // Get the pointer to the track
            G4Track *track = step->GetTrack();

            const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
            const CLHEP::Hep3Vector &stop  = track->GetPosition();

            // If it's a null step, don't use it.
            if(start == stop) return;
            if(step->GetTotalEnergyDeposit() == 0) return;

            // check that we are in the correct material to record a hit
            std::string VolumeName   = this->GetVolumeName(track);

            // if( std::find( fECALVolumeName.begin(), fECALVolumeName.end(), VolumeName ) == fECALVolumeName.end() )
            // return;

            // check the material
            auto pos = 0.5 * (start + stop);
            TGeoNode *node = fGeoManager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm...

            if(!node){
                LOG_DEBUG("AuxDetAction::ECALSteppingAction")
                << "Node not found in "
                << pos.x() << " mm "
                << pos.y() << " mm "
                << pos.z() << " mm";
                return;
            }

            std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
            if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) return;

            // only worry about energy depositions larger than the minimum required
            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "In volume "
            << VolumeName
            << " Material is "
            << volmaterial
            << " step size is "
            << step->GetStepLength() / CLHEP::cm
            << " cm and deposited "
            << step->GetTotalEnergyDeposit()
            << " MeV of energy";

            // the step mid point is used for the position of the deposit
            G4ThreeVector G4Global = 0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition() );

            //get layer and slice number from the volume name (easier for segmentation later)
            unsigned int layer = this->GetLayerNumber(VolumeName); //get layer number
            unsigned int slice = this->GetSliceNumber(VolumeName); // get slice number
            unsigned int det_id = this->GetDetNumber(VolumeName); // 1 == Barrel, 2 = Endcap
            unsigned int stave = this->GetStaveNumber(VolumeName); //get the stave number
            unsigned int module = this->GetModuleNumber(VolumeName); //get the module number

            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "det_id " << det_id
            << " stave " << stave
            << " module " << module
            << " layer " << layer
            << " slice " << slice;

            //Transform from global coordinates to local coordinates
            G4ThreeVector G4Local = this->globalToLocal(step, G4Global);

            //Get cellID
            const ROOT::Math::XYZVector G4local(G4Local.x(), G4Local.y(), G4Local.z());
            long long int cellID = fGeo->cellID(det_id, module, stave, layer, slice, G4local);//encoding the cellID on 64 bits

            // std::cout << "layer before cellID " << layer << std::endl;
            // std::cout << "cellID " << cellID << " cellX " << fGeo->getIDbyCellID(cellID, "cellX") << " cellY " << fGeo->getIDbyCellID(cellID, "cellY") << " layer " << fGeo->getIDbyCellID(cellID, "layer") << std::endl;

            //G4 position
            double G4Pos[3] = { G4Global.x() / CLHEP::cm, G4Global.y() / CLHEP::cm, G4Global.z() / CLHEP::cm };

            // //Calculate position based on cellID in the local frame
            // const ROOT::Math::XYZVector SegLocal = fSegmentation->position(cellID);
            // //Transform from local to global coordinates
            // G4ThreeVector SegGlobal = this->localToGlobal(step, SegLocal);
            // //Segmented hit position
            // double SegPos[3] = { SegGlobal.x() / CLHEP::cm, SegGlobal.y() / CLHEP::cm, SegGlobal.z() / CLHEP::cm };

            // get the track id for this step
            auto trackID = ParticleListAction::GetCurrentTrackID();
            float time = step->GetPreStepPoint()->GetGlobalTime();
            float edep = fGArG4EmSaturation.VisibleEnergyDeposition(step) * CLHEP::MeV / CLHEP::GeV;

            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "Energy deposited "
            << step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV
            << " GeV in cellID "
            << cellID
            << " after Birks "
            << edep
            << " GeV";

            fECALDeposits.emplace_back(trackID,
            time,
            edep,
            G4Pos,
            cellID);
        }

        //------------------------------------------------------------------------------
        std::string AuxDetAction::GetVolumeName(const G4Track *track)
        {
            std::string VolName = track->GetVolume()->GetName();
            VolName.erase(VolName.length()-3, 3);
            return VolName;
        }

        //------------------------------------------------------------------------------
        unsigned int AuxDetAction::GetLayerNumber(std::string volname)
        {
            return std::atoi( (volname.substr( volname.find("layer_") + 6, 2)).c_str() );
        }

        //------------------------------------------------------------------------------
        unsigned int AuxDetAction::GetSliceNumber(std::string volname)
        {
            return std::atoi( (volname.substr( volname.find("slice") + 5, 1)).c_str() );
        }

        //------------------------------------------------------------------------------
        unsigned int AuxDetAction::GetDetNumber(std::string volname)
        {
            unsigned int det_id = 0;

            if( volname.find("Barrel") !=  std::string::npos )
            det_id = 1;
            if( volname.find("Endcap") !=  std::string::npos )
            det_id = 2;

            return det_id;
        }

        //------------------------------------------------------------------------------
        unsigned int AuxDetAction::GetStaveNumber(std::string volname)
        {
            return std::atoi( (volname.substr( volname.find("_stave") + 6, 2)).c_str() );
        }

        //------------------------------------------------------------------------------
        unsigned int AuxDetAction::GetModuleNumber(std::string volname)
        {
            return std::atoi( (volname.substr( volname.find("_module") + 7, 2)).c_str() );
        }

        //------------------------------------------------------------------------------
        G4ThreeVector AuxDetAction::globalToLocal(const G4Step* step, const G4ThreeVector& glob)
        {
            return step->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(glob);
        }

        //------------------------------------------------------------------------------
        G4ThreeVector AuxDetAction::localToGlobal(const G4Step* step, const ROOT::Math::XYZVector& loc)
        {
            return localToGlobal( step, G4ThreeVector(loc.X(), loc.Y(), loc.Z()) );
        }

        //------------------------------------------------------------------------------
        G4ThreeVector AuxDetAction::localToGlobal(const G4Step* step, const G4ThreeVector& loc)
        {
            return step->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().Inverse().TransformPoint(loc);
        }

    } // garg4

} // gar
