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

            class Settings
            {
            public:
                Settings();

                std::string    m_mcParticleCollection;              ///< The mc particle collection
                std::string    m_CaloHitRelationCollection;         ///< The SimCaloHit to CaloHit particle relation
                std::string    m_TrackRelationCollection;           ///< The SimTrackerHit to TrackerHit particle relation
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