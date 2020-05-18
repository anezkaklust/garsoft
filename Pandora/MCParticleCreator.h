#ifndef MCPARTICLECREATOR_H
#define MCPARTICLECREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Api/PandoraApi.h"

#include "CaloHitCreator.h"
#include "TrackCreator.h"

namespace gar {
    namespace gar_pandora {

        class MCParticleCreator
        {
        public:
            typedef std::vector<std::string> StringVector;

            class Settings
            {
            public:
                Settings();

                StringVector    m_mcParticleCollections;                ///< The mc particle collections
                StringVector    m_CaloHitRelationCollections;         ///< The SimCaloHit to CaloHit particle relations
                StringVector    m_TrackRelationCollections;           ///< The SimTrackerHit to TrackerHit particle relations
            };

            MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~MCParticleCreator();

            pandora::StatusCode CreateMCParticles(const art::Event *const pEvent) const;

            pandora::StatusCode CreateTrackToMCParticleRelationships(const art::Event *const pEvent, const TrackVector &trackVector) const;
            pandora::StatusCode CreateCaloHitToMCParticleRelationships(const art::Event *const pEvent, const CalorimeterHitVector &calorimeterHitVector) const;

        private:
            const Settings          m_settings;                         ///< The mc particle creator settings
            const pandora::Pandora &m_pandora;                          ///< Reference to the pandora object to create the mc particles
            const float             m_bField;                           ///< The bfield
        };
    }
}

#endif // #ifndef MCPARTICLECREATOR_H