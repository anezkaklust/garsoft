////////////////////////////////////////////////////////////////////////
/// \file  SiPMReadout_module.cc
/// \brief Create the readout from the sipm signals
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_SIPMREADOUT
#define GAR_READOUTSIMULATION_SIPMREADOUT

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
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
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "art/Persistency/Common/PtrMaker.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "ReadoutSimulation/MuIDReadoutSimStandardAlg.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// Forward declarations

///Geant4 interface
namespace gar {
    namespace rosim {

        /**
        * @brief Runs simulation including smearing of energy and SiPM saturation for the ECAL readout
        *
        *
        * Randomness
        * -----------
        *
        * The random number generators used by this process are:
        * - 'GEANT' instance: used by Geant4
        * - 'propagation' instance: used in electron propagation
        *
        */
        class SiPMReadout : public ::art::EDProducer{
        public:

            /// Standard constructor and destructor for an FMWK module.
            explicit SiPMReadout(fhicl::ParameterSet const& pset);
            virtual ~SiPMReadout();

            // Plugins should not be copied or assigned.
            SiPMReadout(SiPMReadout const &) = delete;
            SiPMReadout(SiPMReadout &&) = delete;
            SiPMReadout & operator = (SiPMReadout const &) = delete;
            SiPMReadout & operator = (SiPMReadout &&) = delete;

            void produce (::art::Event& evt) override;

            void reconfigure(fhicl::ParameterSet const& pset);

        private:

            void CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

            std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

            std::string                         fG4Label;    ///< label of G4 module
            std::string                         fECALInstanceName; ///< product instance name for the ECAL
            std::string                         fMuIDInstanceName; ///< product instance name for the MuID

            const gar::geo::GeometryCore*       fGeo;        ///< geometry information
            std::unique_ptr<SiPMReadoutSimAlg>  fROECALSimAlg;   ///< algorithm to simulate the electronics of the ECAL
            std::unique_ptr<SiPMReadoutSimAlg>  fROMuIDSimAlg;   ///< algorithm to simulate the electronics of the MuID

            CLHEP::HepRandomEngine              &fEngine;  ///< random engine
        };

    } // namespace rosim

    namespace rosim {

        //----------------------------------------------------------------------
        // Constructor
        SiPMReadout::SiPMReadout(fhicl::ParameterSet const& pset)
        : art::EDProducer{pset},
        fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,pset,"Seed"))
        {
            fGeo = gar::providerFrom<geo::Geometry>();

            this->reconfigure(pset);

            art::InputTag ecaltag(fG4Label, fECALInstanceName);
            consumes< std::vector<sdp::CaloDeposit> >(ecaltag);
            produces< std::vector<raw::CaloRawDigit> >(fECALInstanceName);
            produces< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>  >(fECALInstanceName);

            if(fGeo->HasMuonDetector()) {
                art::InputTag muidtag(fG4Label, fMuIDInstanceName);
                consumes< std::vector<sdp::CaloDeposit> >(muidtag);
                produces< std::vector<raw::CaloRawDigit> >(fMuIDInstanceName);//for the MuID
                produces< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>  >(fMuIDInstanceName);//for the MuID
            }
            return;
        }

        //----------------------------------------------------------------------
        // Destructor
        SiPMReadout::~SiPMReadout()
        {
        }

        //----------------------------------------------------------------------
        void SiPMReadout::reconfigure(fhicl::ParameterSet const& pset)
        {
            MF_LOG_DEBUG("SiPMReadout") << "Debug: SiPMReadout()";

            fG4Label = pset.get<std::string >("G4ModuleLabel", "geant");
            fECALInstanceName =  pset.get<std::string >("ECALInstanceName", "");

            auto ECALROAlgPars = pset.get<fhicl::ParameterSet>("ECALReadoutSimAlgPars");
            auto ECALROAlgName = ECALROAlgPars.get<std::string>("ECALReadoutSimType");

            if(ECALROAlgName.compare("Standard") == 0) {
                fROECALSimAlg = std::make_unique<gar::rosim::ECALReadoutSimStandardAlg>(fEngine, ECALROAlgPars);
            }
            else {
                throw cet::exception("SiPMReadout")
                << "Unable to determine which ECAL readout simulation algorithm to use, bail";
            }

            //Muon ID
            fMuIDInstanceName =  pset.get<std::string >("MuIDInstanceName", "");

            auto MuIDROAlgPars = pset.get<fhicl::ParameterSet>("MuIDReadoutSimAlgPars");
            auto MuIDROAlgName = MuIDROAlgPars.get<std::string>("MuIDReadoutSimType");

            if(MuIDROAlgName.compare("Standard") == 0) {
                fROMuIDSimAlg = std::make_unique<gar::rosim::MuIDReadoutSimStandardAlg>(fEngine, MuIDROAlgPars);
            }
            else {
                throw cet::exception("SiPMReadout")
                << "Unable to determine which MuID readout simulation algorithm to use, bail";
            }

            return;
        }

        //--------------------------------------------------------------------------
        void SiPMReadout::produce(::art::Event& evt)
        {
            MF_LOG_DEBUG("SiPMReadout") << "produce()";

            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<sdp::CaloDeposit> > artECALHits;
            this->CollectHits(evt, fG4Label, fECALInstanceName, artECALHits);

            //Get contributions per cellID for association to digi hit
            std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > m_cIDMapArtPtrVec = this->MakeCellIDMapArtPtr(artECALHits);

            //Pass the sim hits to the algo
            fROECALSimAlg->PrepareAlgo(artECALHits);

            //Perform the digitization
            fROECALSimAlg->DoDigitization();

            //Get the digitized hits
            std::vector< raw::CaloRawDigit* > digiECALVec = fROECALSimAlg->GetDigitizedHits();

            // loop over the lists and put the particles and voxels into the event as collections
            std::unique_ptr< std::vector<raw::CaloRawDigit> > digitECALCol (new std::vector<raw::CaloRawDigit>);
            std::unique_ptr< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit> > DigiSimECALHitsAssns(new art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>);

            art::PtrMaker<raw::CaloRawDigit> makeDigiECALPtr(evt, fECALInstanceName);

            for(auto const &it : digiECALVec)
            {
                digitECALCol->emplace_back(*it);
                art::Ptr<raw::CaloRawDigit> digiECALPtr = makeDigiECALPtr(digitECALCol->size() - 1);
                //get the associated vector of art ptr based on cellID
                std::vector< art::Ptr<sdp::CaloDeposit> > simPtrVec = m_cIDMapArtPtrVec[it->CellID()];
                for(auto hitpointer : simPtrVec)
                DigiSimECALHitsAssns->addSingle(digiECALPtr, hitpointer);
            }

            evt.put(std::move(digitECALCol), fECALInstanceName);
            evt.put(std::move(DigiSimECALHitsAssns), fECALInstanceName);

            //Treat the MuID hits
            if(fGeo->HasMuonDetector())
            {
                std::vector< art::Ptr<sdp::CaloDeposit> > artMuIDHits;
                this->CollectHits(evt, fG4Label, fMuIDInstanceName, artMuIDHits);

                //Get contributions per cellID for association to digi hit
                m_cIDMapArtPtrVec = this->MakeCellIDMapArtPtr(artMuIDHits);

                //Pass the sim hits to the algo
                fROMuIDSimAlg->PrepareAlgo(artMuIDHits);

                //Perform the digitization
                fROMuIDSimAlg->DoDigitization();

                //Get the digitized hits
                std::vector< raw::CaloRawDigit* > digiMuIDVec = fROMuIDSimAlg->GetDigitizedHits();

                // loop over the lists and put the particles and voxels into the event as collections
                std::unique_ptr< std::vector<raw::CaloRawDigit> > digitMuIDCol (new std::vector<raw::CaloRawDigit>);
                std::unique_ptr< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit> > DigiSimMuIDHitsAssns(new art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>);

                art::PtrMaker<raw::CaloRawDigit> makeDigiMuIDPtr(evt, fMuIDInstanceName);

                for(auto const &it : digiMuIDVec)
                {
                    digitMuIDCol->emplace_back(*it);
                    art::Ptr<raw::CaloRawDigit> digiMuIDPtr = makeDigiMuIDPtr(digitMuIDCol->size() - 1);
                    //get the associated vector of art ptr based on cellID
                    std::vector< art::Ptr<sdp::CaloDeposit> > simPtrVec = m_cIDMapArtPtrVec[it->CellID()];
                    for(auto hitpointer : simPtrVec)
                    DigiSimMuIDHitsAssns->addSingle(digiMuIDPtr, hitpointer);
                }

                evt.put(std::move(digitMuIDCol), fMuIDInstanceName);
                evt.put(std::move(DigiSimMuIDHitsAssns), fMuIDInstanceName);
            }

            return;
        } // CaloReadout::produce()

        //--------------------------------------------------------------------------
        void SiPMReadout::CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
        {
            art::Handle< std::vector<sdp::CaloDeposit> > theHits;
            evt.getByLabel(label, instance, theHits);

            if (!theHits.isValid())
            return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<sdp::CaloDeposit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        //--------------------------------------------------------------------------
        std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > SiPMReadout::MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
        {
            std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > cIDMapArtPtrVec;

            for (std::vector< art::Ptr<sdp::CaloDeposit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<sdp::CaloDeposit> hitPtr = *iter;
                const sdp::CaloDeposit *hit = hitPtr.get();

                if( cIDMapArtPtrVec.count(hit->CellID()) == 0 )
                {
                    std::vector< art::Ptr<sdp::CaloDeposit> > vecArtPtr;
                    vecArtPtr.push_back(hitPtr);
                    cIDMapArtPtrVec.insert( std::make_pair(hit->CellID(), vecArtPtr) );
                }
                else
                {
                    cIDMapArtPtrVec[hit->CellID()].push_back(hitPtr);
                }
            }

            return cIDMapArtPtrVec;
        }


    } // namespace rosim

    namespace rosim {

        DEFINE_ART_MODULE(SiPMReadout)

    } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_SIPMREADOUT
