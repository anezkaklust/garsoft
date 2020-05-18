#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Api/PandoraApi.h"

namespace gar {
    namespace gar_pandora {

        class PfoCreator
        {
        public:

            class Settings
            {
            public:
                Settings();

                std::string     m_clusterCollectionName = "";            ///< The name of the cluster output collection
                std::string     m_pfoCollectionName = "";                ///< The name of the pfo output collection
                std::string     m_startVertexCollectionName = "";        ///< The name of the start vertex output collection
                std::string     m_startVertexAlgName = "";               ///< The name of the algorithm to fill the start vertex output collection
                float           m_emStochasticTerm = 0;                  ///< The stochastic term for EM shower energy resolution
                float           m_hadStochasticTerm = 0;                 ///< The stochastic term for Hadronic shower energy resolution
                float           m_emConstantTerm = 0;                    ///< The constant term for EM shower energy resolution
                float           m_hadConstantTerm = 0;                   ///< The constant term for Hadronic shower energy resolution
            };

            PfoCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~PfoCreator();

            pandora::StatusCode CreateParticleFlowObjects(const art::Event *const pEvent);
            
        private:
            enum Index
            {
                ECAL_INDEX = 0
            };

            void InitialiseSubDetectorNames(pandora::StringVector &subDetectorNames) const;

            const Settings              m_settings;                         ///< The pfo creator settings
            const pandora::Pandora      &m_pandora;                        ///< Reference to the pandora object from which to extract the pfos
        };
    }
}

#endif // #ifndef PFO_CREATOR_H