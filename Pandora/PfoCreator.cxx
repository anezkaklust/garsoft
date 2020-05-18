#include "PfoCreator.h"

#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

namespace gar {
    namespace gar_pandora {

        PfoCreator::PfoCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pandora(*pPandora)
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        PfoCreator::~PfoCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode PfoCreator::CreateParticleFlowObjects(const art::Event *const pEvent)
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void PfoCreator::InitialiseSubDetectorNames(pandora::StringVector &subDetectorNames) const
        {
            subDetectorNames.push_back("ecal");
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        PfoCreator::Settings::Settings():
        m_clusterCollectionName(""),
        m_pfoCollectionName(""),
        m_startVertexCollectionName(""),
        m_startVertexAlgName(""),
        m_emStochasticTerm(0.17f),
        m_hadStochasticTerm(0.6f),
        m_emConstantTerm(0.01f),
        m_hadConstantTerm(0.03f)
        {
        }
    }
}
