/* 
Art producer module to run Pandora on GArSoft data 
Author: Eldwan Brianne
Date 15.05.2020
*/

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nutools/MagneticField/MagneticField.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "CaloHitCreator.h"
#include "GeometryCreator.h"
#include "MCParticleCreator.h"
#include "PfoCreator.h"
#include "TrackCreator.h"
#include "BFieldPlugin.h"
#include "PseudoLayerPlugin.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include "ReconstructionDataProducts/Cluster.h"

//Pandora
#include "Api/PandoraApi.h"

//LCContent
#include "LCContent.h"
#include "LCPlugins/LCSoftwareCompensation.h"
#include "LCPlugins/LCBFieldPlugin.h"
#include "LCPlugins/LCEnergyCorrectionPlugins.h"
#include "LCPlugins/LCParticleIdPlugins.h"
#include "LCPlugins/LCPseudoLayerPlugin.h"
#include "LCPlugins/LCShowerProfilePlugin.h"

#include <exception>

namespace gar {
    namespace gar_pandora {

        class PandoraInterface : public art::EDProducer {

        public:

            class Settings
            {
            public:
                /**
                 *  @brief  Default constructor
                 */
                Settings();

                std::string            m_pandoraSettingsXmlFile = "";      ///< The pandora settings xml file
                float                  m_innerBField = 0.0;                ///< The bfield in the main tracker and ecal, units Tesla
                std::vector<float>     m_inputEnergyCorrectionPoints{};    ///< The input energy points for non-linearity energy correction
                std::vector<float>     m_outputEnergyCorrectionPoints{};   ///< The output energy points for non-linearity energy correction
            };

            explicit PandoraInterface(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            PandoraInterface(PandoraInterface const &) = delete;
            PandoraInterface(PandoraInterface &&) = delete;
            PandoraInterface & operator = (PandoraInterface const &) = delete;
            PandoraInterface & operator = (PandoraInterface &&) = delete;

            // Required functions.
            void beginJob() override;
            void produce(art::Event & e) override;
            void endJob() override;

        private:
            void reconfigure(fhicl::ParameterSet const& pset);
            void FinaliseSteeringParameters();

            pandora::StatusCode RegisterUserComponents() const;
            void Reset();
            void printParameters();

            pandora::Pandora                   *m_pPandora = nullptr;                 ///< Address of the pandora instance
            CaloHitCreator                     *m_pCaloHitCreator = nullptr;          ///< The calo hit creator
            GeometryCreator                    *m_pGeometryCreator = nullptr;         ///< The geometry creator
            TrackCreator                       *m_pTrackCreator = nullptr;            ///< The track creator
            MCParticleCreator                  *m_pMCParticleCreator = nullptr;       ///< The mc particle creator
            PfoCreator                         *m_pPfoCreator = nullptr;              ///< The pfo creator

            Settings                           m_settings{};                       ///< The settings for the pandora interface module
            CaloHitCreator::Settings           m_caloHitCreatorSettings{};         ///< The calo hit creator settings
            GeometryCreator::Settings          m_geometryCreatorSettings{};        ///< The geometry creator settings
            MCParticleCreator::Settings        m_mcParticleCreatorSettings{};      ///< The mc particle creator settings
            TrackCreator::Settings             m_trackCreatorSettings{};           ///< The track creator settings
            PfoCreator::Settings               m_pfoCreatorSettings{};             ///< The pfo creator settings
            PseudoLayerPlugin::Settings        m_pseudoLayerParameters{};          ///< The pseudo layer settings

            typedef std::map<const pandora::Pandora *, art::Event *> PandoraToEventMap;
            static PandoraToEventMap          m_pandoraToEventMap;              ///< The pandora to event map

            const geo::GeometryCore* fGeo; ///< pointer to the geometry
        };

        PandoraInterface::PandoraToEventMap PandoraInterface::m_pandoraToEventMap;

        //-----------------------------------------------------------------------------------
        PandoraInterface::PandoraInterface(fhicl::ParameterSet const& pset) :
        EDProducer{pset}
        {
            this->reconfigure(pset);
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::beginJob()
        {
            printParameters();

            try
            {
                LOG_INFO("PandoraInterface - beginJob");

                this->FinaliseSteeringParameters();

                m_pPandora = new pandora::Pandora();
                m_pGeometryCreator = new GeometryCreator(m_geometryCreatorSettings, m_pPandora);
                m_pCaloHitCreator = new CaloHitCreator(m_caloHitCreatorSettings, m_pPandora);
                m_pTrackCreator = new TrackCreator(m_trackCreatorSettings, m_pPandora);
                m_pMCParticleCreator = new MCParticleCreator(m_mcParticleCreatorSettings, m_pPandora);
                m_pPfoCreator = new PfoCreator(m_pfoCreatorSettings, m_pPandora);

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry());
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
            }
            catch (pandora::StatusCodeException &statusCodeException)
            {
                LOG_ERROR("Failed to initialize PandoraInterface: ") << statusCodeException.ToString();
                throw statusCodeException;
            }
            catch (std::exception &exception)
            {
                LOG_ERROR("Failed to initialize PandoraInterface: std exception ") << exception.what();
                throw exception;
            }
            catch (...)
            {
                LOG_ERROR("Failed to initialize PandoraInterface: unrecognized exception");
                throw;
            }
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::produce(art::Event & e)
        {
            try
            {
                LOG_INFO("PandoraInterface - produce");
                (void) m_pandoraToEventMap.insert(PandoraToEventMap::value_type(m_pPandora, &e));

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateMCParticles(&e));
                // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(e));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(&e));
                // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateTrackToMCParticleRelationships(e, m_pTrackCreator->GetTrackVector()));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(&e));
                // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateCaloHitToMCParticleRelationships(e, m_pCaloHitCreator->GetCalorimeterHitVector()));

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pPfoCreator->CreateParticleFlowObjects(&e));

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
                this->Reset();
            }
            catch (pandora::StatusCodeException &statusCodeException)
            {
                LOG_ERROR("Marlin pandora failed to process event: ") << statusCodeException.ToString();
                throw statusCodeException;
            }
            catch (std::exception &exception)
            {
                LOG_ERROR("Marlin pandora failed to process event: std exception ") << exception.what();
                throw exception;
            }
            catch (...)
            {
                LOG_ERROR("Marlin pandora failed to process event: unrecognized exception");
                throw;
            }   
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::endJob()
        {
            delete m_pPandora;
            delete m_pGeometryCreator;
            delete m_pCaloHitCreator;
            delete m_pTrackCreator;
            delete m_pMCParticleCreator;
            delete m_pPfoCreator;
        }

        //-----------------------------------------------------------------------------------
        pandora::StatusCode PandoraInterface::RegisterUserComponents() const
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));

            //Custom Plugins
            art::ServiceHandle<mag::MagneticField> magFieldService;
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldPlugin(*m_pPandora, new BFieldPlugin(magFieldService)));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pPandora, new PseudoLayerPlugin(m_pseudoLayerParameters)));
            //Shower Profile Plugin
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfilePlugin(*m_pPandora, new lc_content::LCShowerProfilePlugin));

            //Energy Correction Plugins
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterEnergyCorrectionPlugin(*m_pPandora, "CleanClusters", pandora::HADRONIC, new lc_content::LCEnergyCorrectionPlugins::CleanCluster));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterEnergyCorrectionPlugin(*m_pPandora, "ScaleHotHadrons", pandora::HADRONIC, new lc_content::LCEnergyCorrectionPlugins::ScaleHotHadrons));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterEnergyCorrectionPlugin(*m_pPandora, "MuonCoilCorrection", pandora::HADRONIC, new lc_content::LCEnergyCorrectionPlugins::MuonCoilCorrection));

            //PID Plugins
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdPlugin(*m_pPandora, "LCEmShowerId" , new lc_content::LCParticleIdPlugins::LCEmShowerId));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdPlugin(*m_pPandora, "LCPhotonId" , new lc_content::LCParticleIdPlugins::LCPhotonId));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdPlugin(*m_pPandora, "LCElectronId" , new lc_content::LCParticleIdPlugins::LCElectronId));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdPlugin(*m_pPandora, "LCMuonId" , new lc_content::LCParticleIdPlugins::LCMuonId));


            return pandora::STATUS_CODE_SUCCESS;
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::reconfigure(fhicl::ParameterSet const& pset)
        {
            //Find pandora xml config file
            cet::search_path sp("FW_SEARCH_PATH");
            std::string fullConfigFileName;
            std::string pandora_file = pset.get<std::string>("PandoraXML");

            if( !sp.find_file(pandora_file, fullConfigFileName) )
                throw cet::exception("PandoraInterface") << " Failed to find xml configuration file " << pandora_file << " in FW search path";
            m_settings.m_pandoraSettingsXmlFile = fullConfigFileName;

            m_trackCreatorSettings.m_trackCollection = pset.get<std::string>("TrackLabel", "track"); //Track hits
            m_trackCreatorSettings.m_minTrackHits = pset.get<unsigned int>("MinTrackHits", 5); //Track quality cut: the minimum number of track hits
            m_trackCreatorSettings.m_maxTrackHits =  pset.get<unsigned int>("MaxTrackHits", 5000); //Track quality cut: the maximum number of track hits
            m_trackCreatorSettings.m_curvatureToMomentumFactor = pset.get<float>("CurvatureToMomentumFactor", 0.3/2000.); //Constant relating track curvature in b field to momentum
            m_trackCreatorSettings.m_maxTrackSigmaPOverP = pset.get<unsigned int>("MaxTrackSigmaPOverP", 3); //Cut on fractional track momentum error
            m_trackCreatorSettings.m_minMomentumForTrackHitChecks = pset.get<float>("MinMomentumForTrackHitChecks", 0.5); //Min track momentum required to perform final quality checks on number of hits

            m_pfoCreatorSettings.m_emStochasticTerm = pset.get<float>("EMStockTerm", 0.05); //The stochastic term for EM shower
            m_pfoCreatorSettings.m_emConstantTerm = pset.get<float>("EMConstTerm", 2.); //The constant term for EM shower

            m_caloHitCreatorSettings.m_CaloHitCollection = pset.get<std::string>("CaloHitLabel", "sscalohit"); //Calo hits
            m_caloHitCreatorSettings.m_eCalToMip = pset.get<float>("ECaltoMipCalibration", 1.); //The calibration from deposited ECal energy to mip
            m_caloHitCreatorSettings.m_eCalMipThreshold = pset.get<float>("ECalMipThreshold", 0.25); //Threshold for creating calo hits in the ECal, units mip
            m_caloHitCreatorSettings.m_eCalToEMGeV = pset.get<float>("ECalToEMGeVCalibration", 1.); //The calibration from deposited ECal energy to EM energy
            m_caloHitCreatorSettings.m_eCalToHadGeVBarrel = pset.get<float>("ECalToHadGeVCalibrationBarrel", 1.); //The calibration from deposited ECal energy to hadronic energy in the barrel
            m_caloHitCreatorSettings.m_eCalToHadGeVEndCap = pset.get<float>("ECalToHadGeVCalibrationEndCap", 1.); //The calibration from deposited ECal energy to hadronic energy in the endcap
            m_caloHitCreatorSettings.m_maxECalHitHadronicEnergy = pset.get<float>("MaxECalHitHadronicEnergy", 10000.); //The maximum hadronic energy allowed for a single ecal hit
            m_caloHitCreatorSettings.m_nOuterSamplingLayers = pset.get<unsigned int>("NOuterSamplingLayers", 3); //Number of layers from edge for hit to be flagged as an outer layer hit
            m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance = pset.get<float>("LayersFromEdgeMaxRearDistance", 250.); //Maximum number of layers from candidate outer layer hit to rear of detector
            m_caloHitCreatorSettings.m_eCalBarrelNormalVector = pset.get<std::vector<float>>("ECalBarrelNormalVector", std::vector<float>{0., 0., 1.}); //Normal vector for the ECal barrel sensitive layers in local coordinates

            // produces< std::vector<gar::rec::Cluster> >();
            // produces< std::vector<gar::rec::PFParticle> >();
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::FinaliseSteeringParameters()
        {
            art::ServiceHandle<mag::MagneticField> magFieldService;
            G4ThreeVector zerovec(0,0,0);
            G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);
            m_settings.m_innerBField = magfield[0]; //x component at (0, 0, 0)

            fGeo = gar::providerFrom<geo::Geometry>();

            m_trackCreatorSettings.m_bField                         = magfield[0];
            m_trackCreatorSettings.m_eCalBarrelInnerSymmetry        = fGeo->GetECALInnerSymmetry();
            m_trackCreatorSettings.m_eCalBarrelInnerPhi0            = 0.;
            m_trackCreatorSettings.m_eCalBarrelInnerR               = fGeo->GetECALInnerBarrelRadius();
            m_trackCreatorSettings.m_eCalEndCapInnerX               = fGeo->GetECALEndcapStartX();

            m_caloHitCreatorSettings.m_eCalBarrelInnerPhi0          = 0.;
            m_caloHitCreatorSettings.m_eCalBarrelInnerSymmetry      = fGeo->GetECALInnerSymmetry();
            m_caloHitCreatorSettings.m_eCalBarrelOuterR             = fGeo->GetECALOuterBarrelRadius();
            m_caloHitCreatorSettings.m_eCalEndCapOuterX             = fGeo->GetECALEndcapOuterX();
            m_caloHitCreatorSettings.m_eCalEndCapOuterR             = fGeo->GetECALOuterEndcapRadius();
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::Reset()
        {
            m_pCaloHitCreator->Reset();
            m_pTrackCreator->Reset();

            PandoraToEventMap::iterator iter = m_pandoraToEventMap.find(m_pPandora);

            if (m_pandoraToEventMap.end() == iter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            m_pandoraToEventMap.erase(iter);
        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::printParameters()
        {

        }

        //-----------------------------------------------------------------------------------
        PandoraInterface::Settings::Settings() :
        m_innerBField(0.5),
        m_inputEnergyCorrectionPoints(0),
        m_outputEnergyCorrectionPoints(0)
        {
        }

    } //namespace gar_pandora
} // namespace gar

namespace gar {
    namespace gar_pandora {
        DEFINE_ART_MODULE(PandoraInterface)
    }
}