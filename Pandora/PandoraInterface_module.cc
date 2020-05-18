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

#include "fhiclcpp/ParameterSet.h"

#include "CaloHitCreator.h"
#include "GeometryCreator.h"
#include "MCParticleCreator.h"
#include "PfoCreator.h"
#include "TrackCreator.h"
#include "BFieldPlugin.h"
#include "PseudoLayerPlugin.h"

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

                std::string     m_pandoraSettingsXmlFile = "";      ///< The pandora settings xml file

                float           m_innerBField = 0.0;                ///< The bfield in the main tracker, ecal and hcal, units Tesla

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
            pandora::StatusCode RegisterUserComponents() const;
            void ProcessSteeringFile();
            void FinaliseSteeringParameters();
            void Reset();
            void printParameters();

            pandora::Pandora                   *m_pPandora = NULL;                  ///< Address of the pandora instance
            CaloHitCreator                     *m_pCaloHitCreator = NULL;           ///< The calo hit creator
            GeometryCreator                    *m_pGeometryCreator = NULL;          ///< The geometry creator
            TrackCreator                       *m_pTrackCreator = NULL;             ///< The track creator
            MCParticleCreator                  *m_pMCParticleCreator = NULL;      ///< The mc particle creator
            PfoCreator                         *m_pPfoCreator = NULL;             ///< The pfo creator

            Settings                            m_settings{};                       ///< The settings for the pandora interface module
            CaloHitCreator::Settings            m_caloHitCreatorSettings{};         ///< The calo hit creator settings
            GeometryCreator::Settings           m_geometryCreatorSettings{};        ///< The geometry creator settings
            MCParticleCreator::Settings         m_mcParticleCreatorSettings{};      ///< The mc particle creator settings
            TrackCreator::Settings              m_trackCreatorSettings{};           ///< The track creator settings
            PfoCreator::Settings                m_pfoCreatorSettings{};             ///< The pfo creator settings
            PseudoLayerPlugin::Settings         m_pseudoLayerParameters{};          ///< The pseudo layer settings

            typedef std::map<const pandora::Pandora *, art::Event *> PandoraToEventMap;
            static PandoraToEventMap          m_pandoraToEventMap;              ///< The pandora to event map
        };

        PandoraInterface::PandoraToEventMap PandoraInterface::m_pandoraToEventMap;

        //-----------------------------------------------------------------------------------
        PandoraInterface::PandoraInterface(fhicl::ParameterSet const& pset) :
        EDProducer{pset}
        {
            cet::search_path sp("FW_SEARCH_PATH");
            std::string fullConfigFileName;

            if( !sp.find_file(pset.get<std::string>("PandoraXML"), fullConfigFileName) )
                throw cet::exception("PandoraInterface") << " Failed to find xml configuration file " << pset.get<std::string>("PandoraXML") << " in FW search path";
            m_settings.m_pandoraSettingsXmlFile = fullConfigFileName;

            this->ProcessSteeringFile();
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
            // PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBasicPlugins(*m_pPandora));

            //Custom Plugins
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pPandora, new PseudoLayerPlugin(m_pseudoLayerParameters)));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldPlugin(*m_pPandora, new BFieldPlugin()));
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
        void PandoraInterface::ProcessSteeringFile()
        {

        }

        //-----------------------------------------------------------------------------------
        void PandoraInterface::FinaliseSteeringParameters()
        {

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
        m_innerBField(3.5f),
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