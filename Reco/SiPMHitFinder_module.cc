////////////////////////////////////////////////////////////////////////
// Class:       SiPMHitFinder
// Plugin Type: producer (art v2_11_02)
// File:        CompressedHitFinder_module.cc
//
// Generated at Wed Jun 13 15:27:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Report the compressed raw digit blocks as hits -- just a threshold
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/LocalTransformation.h"
#include "Geometry/BitFieldCoder.h"

#include "RecoAlg/SiPMUtils.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TGeoNode.h"
#include "TGeoManager.h"

#include <memory>

namespace gar {
    namespace rec {

        class SiPMHitFinder : public art::EDProducer {
        public:
            explicit SiPMHitFinder(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            SiPMHitFinder(SiPMHitFinder const &) = delete;
            SiPMHitFinder(SiPMHitFinder &&) = delete;
            SiPMHitFinder & operator = (SiPMHitFinder const &) = delete;
            SiPMHitFinder & operator = (SiPMHitFinder &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        protected:
            float CalibrateToMIP(unsigned int ADC);

            float CalibratetoMeV(float x, float y, float z, double MIP);

            std::array<double, 3U> CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID);

            float CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID);

        private:

            // Declare member data here.
            std::string fRawDigitLabel;  ///< label to find the right raw digits
            std::string fECALInstanceName; ///< product instance name for the ECAL
            std::string fMuIDInstanceName; ///< product instance name for the MuID

            float fMIPThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
            bool fDesaturation; ///< flag to perform the SiPM desaturation
            float fECALfactorSamplingGeV;
            bool fUseTimePositionReco;

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
            std::unique_ptr<util::SiPMUtils>    fSiPMUtils;    ///< pointer to the Util fcn containing the SiPM desaturation function
            gar::geo::BitFieldCoder const* fFieldDecoder;
        };


        SiPMHitFinder::SiPMHitFinder(fhicl::ParameterSet const & p) : EDProducer{p}
        // :
        {
            fRawDigitLabel = p.get<std::string>("RawDigitLabel", "daqecal");
            fECALInstanceName =  p.get<std::string >("ECALInstanceName", "");

            fMIPThreshold = p.get<float>("MIPThreshold", 0.25);
            fDesaturation = p.get<bool>("Desaturation", false);
            fECALfactorSamplingGeV = p.get<float>("ECALSamplingFactorGeV", 1/0.24);
            fUseTimePositionReco = p.get<bool>("UseTimePositionReco", true);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fSiPMUtils = std::make_unique<util::SiPMUtils>(fDetProp->EffectivePixel());

            std::string fECALEncoding = fGeo->GetECALCellIDEncoding();
            fFieldDecoder = new gar::geo::BitFieldCoder( fECALEncoding );

            art::InputTag ecaltag(fRawDigitLabel, fECALInstanceName);
            consumes< std::vector<raw::CaloRawDigit> >(ecaltag);
            produces< std::vector<rec::CaloHit> >(fECALInstanceName);
            produces< art::Assns<rec::CaloHit, raw::CaloRawDigit>  >(fECALInstanceName);

            //Muon ID
            fMuIDInstanceName =  p.get<std::string >("MuIDInstanceName", "");

            if(fGeo->HasMuonDetector()) {
                art::InputTag muidtag(fRawDigitLabel, fMuIDInstanceName);
                consumes< std::vector<raw::CaloRawDigit> >(muidtag);
                produces< std::vector<rec::CaloHit> >(fMuIDInstanceName);
                produces< art::Assns<rec::CaloHit, raw::CaloRawDigit>  >(fMuIDInstanceName);
            }
        }

        void SiPMHitFinder::produce(art::Event & e)
        {
            // create an emtpy output hit collection -- to add to.
            std::unique_ptr< std::vector<rec::CaloHit> > ecalhitCol ( new std::vector< rec::CaloHit > );
            std::unique_ptr< art::Assns<rec::CaloHit, raw::CaloRawDigit> > ecalRecoDigiHitsAssns( new art::Assns<rec::CaloHit, raw::CaloRawDigit> );

            art::PtrMaker<rec::CaloHit> makeCaloHitPtr(e, fECALInstanceName);

            // the input raw digits
            art::InputTag ecaltag(fRawDigitLabel, fECALInstanceName);
            auto digiCol = e.getValidHandle< std::vector<raw::CaloRawDigit> >(ecaltag);

            for (size_t idigit = 0; idigit < digiCol->size(); ++idigit)
            {
                const raw::CaloRawDigit& digitHit = (*digiCol)[idigit];

                unsigned int hitADC = digitHit.ADC().first;
                std::pair<float, float> hitTime = digitHit.Time();
                float x = digitHit.X();
                float y = digitHit.Y();
                float z = digitHit.Z();
                raw::CellID_t cellID = digitHit.CellID();

                //Do Calibration of the hit in MIPs
                float hitMIP = this->CalibrateToMIP(hitADC);

                if(hitMIP < fMIPThreshold)
                {
                    MF_LOG_DEBUG("SiPMHitFinder") << "Signal under the " << fMIPThreshold << " MIP threshold" << std::endl;
                    continue;
                }

                //Desaturation
                float energy = 0.;
                if(fDesaturation)
                {
                    double sat_px = hitMIP * fDetProp->LightYield();
                    //DeSaturate
                    double unsat_px = fSiPMUtils->DeSaturate(sat_px);

                    //Calibrate to the MeV scale
                    double unsat_energy = unsat_px / fDetProp->LightYield();
                    energy = this->CalibratetoMeV(x, y, z, unsat_energy);
                }
                else{
                    //Calibrate to the MeV scale
                    energy = this->CalibratetoMeV(x, y, z, hitMIP);
                }

                //Position reconstruction based on time for the strips
                float pos[3] = {0., 0., 0.};
                std::pair<float, float> time;
                const std::array<double, 3> point = { x, y, z };

                if(fGeo->isTile(point, cellID))
                {
                    pos[0] = x;
                    pos[1] = y;
                    pos[2] = z;
                    time = hitTime;
                }
                else{
                    if(fUseTimePositionReco){
                        std::array<double, 3> strip_pos = this->CalculateStripHitPosition(x, y, z, hitTime, cellID);
                        pos[0] = strip_pos[0];
                        pos[1] = strip_pos[1];
                        pos[2] = strip_pos[2];
                        time = std::make_pair( this->CorrectStripHitTime(pos[0], pos[1], pos[2], hitTime, cellID), 0. );
                    } else {
                        pos[0] = x;
                        pos[1] = y;
                        pos[2] = z;
                        time = hitTime;
                    }
                }

                //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                unsigned int layer = fFieldDecoder->get(cellID, "layer");
                rec::CaloHit hit(fECALfactorSamplingGeV * energy * CLHEP::MeV / CLHEP::GeV, time, pos, cellID, layer);
                ecalhitCol->emplace_back(hit);

                MF_LOG_DEBUG("SiPMHitFinder") << "recohit " << &hit
                << " with cellID " << cellID
                << " has energy " << fECALfactorSamplingGeV * energy * CLHEP::MeV / CLHEP::GeV
                << " pos (" << pos[0] << ", " <<  pos[1] << ", " << pos[2] << ")";

                //Make association between digi hits and reco hits
                art::Ptr<rec::CaloHit> hitPtr = makeCaloHitPtr(ecalhitCol->size() - 1);
                art::Ptr<raw::CaloRawDigit> digiArtPtr = art::Ptr<raw::CaloRawDigit>(digiCol, idigit);
                ecalRecoDigiHitsAssns->addSingle(hitPtr, digiArtPtr);
            }

            //move the reco hit collection
            e.put(std::move(ecalhitCol), fECALInstanceName);
            e.put(std::move(ecalRecoDigiHitsAssns), fECALInstanceName);

            if(fGeo->HasMuonDetector()) {
                // create an emtpy output hit collection -- to add to.
                std::unique_ptr< std::vector<rec::CaloHit> > muidhitCol ( new std::vector< rec::CaloHit > );
                std::unique_ptr< art::Assns<rec::CaloHit, raw::CaloRawDigit> > muidRecoDigiHitsAssns( new art::Assns<rec::CaloHit, raw::CaloRawDigit> );

                art::PtrMaker<rec::CaloHit> makeMuIDHitPtr(e, fMuIDInstanceName);

                // the input raw digits
                art::InputTag muidtag(fRawDigitLabel, fMuIDInstanceName);
                auto digiCol = e.getValidHandle< std::vector<raw::CaloRawDigit> >(muidtag);

                for (size_t idigit = 0; idigit < digiCol->size(); ++idigit)
                {
                    const raw::CaloRawDigit& digitHit = (*digiCol)[idigit];

                    unsigned int hitADC = digitHit.ADC().first;
                    std::pair<float, float> hitTime = digitHit.Time();
                    float x = digitHit.X();
                    float y = digitHit.Y();
                    float z = digitHit.Z();
                    raw::CellID_t cellID = digitHit.CellID();

                    //Do Calibration of the hit in MIPs
                    float hitMIP = this->CalibrateToMIP(hitADC);

                    if(hitMIP < fMIPThreshold)
                    {
                        MF_LOG_DEBUG("SiPMHitFinder") << "Signal under the " << fMIPThreshold << " MIP threshold" << std::endl;
                        continue;
                    }

                    //Desaturation
                    float energy = 0.;
                    if(fDesaturation)
                    {
                        double sat_px = hitMIP * fDetProp->LightYield();
                        //DeSaturate
                        double unsat_px = fSiPMUtils->DeSaturate(sat_px);

                        //Calibrate to the MeV scale
                        double unsat_energy = unsat_px / fDetProp->LightYield();
                        energy = this->CalibratetoMeV(x, y, z, unsat_energy);
                    }
                    else{
                        //Calibrate to the MeV scale
                        energy = this->CalibratetoMeV(x, y, z, hitMIP);
                    }

                    //Position reconstruction based on time for the strips
                    float pos[3] = {0., 0., 0.};
                    std::pair<float, float> time;

                    std::array<double, 3> strip_pos = this->CalculateStripHitPosition(x, y, z, hitTime, cellID);
                    pos[0] = strip_pos[0];
                    pos[1] = strip_pos[1];
                    pos[2] = strip_pos[2];
                    time = std::make_pair( this->CorrectStripHitTime(pos[0], pos[1], pos[2], hitTime, cellID), 0. );

                    //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                    unsigned int layer = fFieldDecoder->get(cellID, "layer");
                    rec::CaloHit hit(energy * CLHEP::MeV / CLHEP::GeV, time, pos, cellID, layer);
                    muidhitCol->emplace_back(hit);

                    MF_LOG_DEBUG("SiPMHitFinder") << "recohit " << &hit
                    << " with cellID " << cellID
                    << " has energy " << energy * CLHEP::MeV / CLHEP::GeV
                    << " pos (" << pos[0] << ", " <<  pos[1] << ", " << pos[2] << ")";

                    //Make association between digi hits and reco hits
                    art::Ptr<rec::CaloHit> hitPtr = makeMuIDHitPtr(muidhitCol->size() - 1);
                    art::Ptr<raw::CaloRawDigit> digiArtPtr = art::Ptr<raw::CaloRawDigit>(digiCol, idigit);
                    muidRecoDigiHitsAssns->addSingle(hitPtr, digiArtPtr);
                }

                //move the reco hit collection
                e.put(std::move(muidhitCol), fMuIDInstanceName);
                e.put(std::move(muidRecoDigiHitsAssns), fMuIDInstanceName);
            }

            return;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CalibrateToMIP(unsigned int ADC)
        {
            float calibrated_ADC = 0.;
            if(ADC <= 0) {
                return calibrated_ADC;
            }

            calibrated_ADC = ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
            return calibrated_ADC;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CalibratetoMeV(float x, float y, float z, double MIP)
        {
            float calibrated_MIP = 0.;
            if(MIP <= 0) {
                return calibrated_MIP;
            }

            TVector3 point(x, y, z);
            float factor = fGeo->GetSensVolumeThickness(point) / 0.5;

            calibrated_MIP = MIP * fDetProp->MeVtoMIP() * factor;
            return calibrated_MIP; // MeV
        }

        //----------------------------------------------------------------------------
        std::array<double, 3> SiPMHitFinder::CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            std::array<double, 3> pointLocal;
            gar::geo::LocalTransformation<TGeoHMatrix> trans;
            fGeo->WorldToLocal(point, pointLocal, trans);

            //Calculate the position of the hit based on the time of both SiPM along the strip
            // pos along strip is
            // x = c * (t1 - t2) / 2 (0 is the center of the strip)
            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float xlocal = c * ( hitTime.first - hitTime.second ) / 2.;

            std::array<double, 3> local_back = fGeo->ReconstructStripHitPosition(point, pointLocal, xlocal, cID);
            std::array<double, 3> world_back;
            fGeo->LocalToWorld(local_back, world_back, trans);

            return world_back;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            double stripLength = fGeo->getStripLength(point, cID); // in cm

            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float time = (hitTime.first + hitTime.second) / 2. - (stripLength / (2 * c));

            return time;
        }

        DEFINE_ART_MODULE(SiPMHitFinder)

    } // namespace rec
} // namespace gar
