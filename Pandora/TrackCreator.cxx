#include "TrackCreator.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "LCObjects/LCTrack.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
    namespace gar_pandora {

        TrackCreator::TrackCreator(const Settings &settings , const pandora::Pandora *const pPandora)
        : m_settings(settings),
        m_pandora(*pPandora)
        {
            m_TrackFactory = std::make_shared<lc_content::LCTrackFactory>();
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        TrackCreator::~TrackCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::CollectTracks(const art::Event &pEvent)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectTracks(pEvent, m_settings.m_trackCollection, artTrkVector));

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::CreateTrackAssociations(const art::Event &pEvent)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(pEvent));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(pEvent));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(pEvent));

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::ExtractKinks(const art::Event &pEvent)
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::ExtractProngsAndSplits(const art::Event &pEvent)
        {
            return pandora::STATUS_CODE_SUCCESS;
        }
        
        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::ExtractV0s(const art::Event &pEvent)
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::CreateTracks() const
        {
            for (TrackVector::const_iterator iter = artTrkVector.begin(), iterEnd = artTrkVector.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<gar::rec::Track> artPtrTrack = *iter;
                const gar::rec::Track *pTrack = artPtrTrack.get();

                const float *trackParams = pTrack->TrackParEnd(); //y, z, curvature, phi, lambda
                const float omega = 1 / trackParams[2];
                const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(omega));

                const unsigned int nTrackHits(static_cast<int>(pTrack->NHits()));
                if ((nTrackHits < m_settings.m_minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
                continue;

                PandoraApi::Track::Parameters trackParameters;

                pandora::CartesianVector startPosition(pTrack->Vertex()[0], pTrack->Vertex()[1], pTrack->Vertex()[2]), endPosition(pTrack->End()[0], pTrack->End()[1], pTrack->End()[2]);
                pandora::CartesianVector calorimeterPosition(0, 0, 0);

                trackParameters.m_pParentAddress = (void*)pTrack;
                trackParameters.m_d0 = 0.f;
                trackParameters.m_z0 = 0.f;

                trackParameters.m_particleId = (omega > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
                trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

                if (std::numeric_limits<float>::epsilon() < std::fabs(omega))
                trackParameters.m_charge = static_cast<int>(omega / std::fabs(omega));

                pandora::CartesianVector momentum = pandora::CartesianVector(std::cos(trackParams[3]), std::sin(trackParams[3]), std::tan(trackParams[4])) * pt;
                trackParameters.m_momentumAtDca = momentum;
                trackParameters.m_trackStateAtStart = pandora::TrackState(startPosition, momentum);
                trackParameters.m_trackStateAtEnd = pandora::TrackState(endPosition, momentum);
                trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(calorimeterPosition, momentum);
                trackParameters.m_timeAtCalorimeter = 0.f;

                trackParameters.m_isProjectedToEndCap   = true;
                trackParameters.m_reachesCalorimeter    = true;
                trackParameters.m_canFormPfo            = true;
                trackParameters.m_canFormClusterlessPfo = true;

                try
                {
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters, *m_TrackFactory));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    LOG_DEBUG("TrackCreator::CreateTracks")
                    << "Failed to extract a track: " << statusCodeException.ToString();
                }

            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        bool TrackCreator::PassesQualityCuts(const gar::rec::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
        {
            return true;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void TrackCreator::TrackReachesECAL(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
        {
            return;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void TrackCreator::DefineTrackPfoUsage(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
        {
            return;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::CollectTracks(const art::Event &pEvent, const std::string &label, TrackVector &trkVector)
        {
            art::Handle< RawTrackVector > theTrk;
            pEvent.getByLabel(label, theTrk);

            if (!theTrk.isValid())
            {
                LOG_DEBUG("TrackCreator::CreateTracks") << "  Failed to find tracks... " << std::endl;
                return pandora::STATUS_CODE_NOT_FOUND;
            }

            LOG_DEBUG("TrackCreator::CreateTracks") << "  Found: " << theTrk->size() << " tracks " << std::endl;

            for (unsigned int i = 0; i < theTrk->size(); ++i)
            {
                const art::Ptr<gar::rec::Track> trk(theTrk, i);
                trkVector.push_back(trk);
            }

            return pandora::STATUS_CODE_SUCCESS;
        }
    }
}
