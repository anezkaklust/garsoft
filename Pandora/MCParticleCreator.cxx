#include "CaloHitCreator.h"
#include "MCParticleCreator.h"
#include "TrackCreator.h"

#include <cmath>
#include <limits>

namespace gar {
    namespace gar_pandora {

        MCParticleCreator::MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pandora(*pPandora),
        m_bField(0.5)
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        MCParticleCreator::~MCParticleCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateMCParticles(const art::Event *const pEvent) const
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const art::Event *const pEvent, const TrackVector &trackVector) const
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const art::Event *const pEvent, const CalorimeterHitVector &calorimeterHitVector) const
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        MCParticleCreator::Settings::Settings():
        m_mcParticleCollection( "" ),
        m_CaloHitRelationCollection( "" ),
        m_TrackRelationCollection( "" )

        {
        }
    }
}