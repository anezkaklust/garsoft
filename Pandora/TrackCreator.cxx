#include "TrackCreator.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
    namespace gar_pandora {

        TrackCreator::TrackCreator(const Settings &settings , const pandora::Pandora *const pPandora, const RotationTransformation *const pRotation)
        : m_settings(settings),
        m_pandora(*pPandora),
        m_rotation(*pRotation)
        {

        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        TrackCreator::~TrackCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode TrackCreator::CreateTracks(const art::Event &pEvent)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectTracks(pEvent, m_settings.m_trackCollection, artTrkVector));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateTracks());

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

                const float *trackParams = pTrack->TrackParEnd(); //y, z, omega, phi, lambda
                const float omega = trackParams[2] / CLHEP::cm;
                const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(omega));

                const unsigned int nTrackHits(static_cast<int>(pTrack->NHits()));
                if ((nTrackHits < m_settings.m_minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
                continue;

                PandoraApi::Track::Parameters trackParameters;

                const pandora::CartesianVector startPosition(pTrack->Vertex()[0], pTrack->Vertex()[1], pTrack->Vertex()[2]);
                const pandora::CartesianVector endPosition(pTrack->End()[0], pTrack->End()[1], pTrack->End()[2]);
                const pandora::CartesianVector calorimeterPosition(0, 0, 0);

                const pandora::CartesianVector newstartPosition = m_rotation.MakeRotation(startPosition);
                const pandora::CartesianVector newendPosition = m_rotation.MakeRotation(endPosition);
                const pandora::CartesianVector newcalorimeterPosition = m_rotation.MakeRotation(calorimeterPosition);

                trackParameters.m_pParentAddress = (void*)pTrack;
                trackParameters.m_d0 = trackParams[0];
                trackParameters.m_z0 = trackParams[1];

                trackParameters.m_particleId = (omega > 0) ? pandora::MU_PLUS : pandora::MU_MINUS;
                trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::MU_PLUS);

                if (std::numeric_limits<float>::epsilon() < std::fabs(omega))
                trackParameters.m_charge = static_cast<int>(omega / std::fabs(omega));

                const pandora::CartesianVector momentum = pandora::CartesianVector(std::cos(trackParams[3]), std::sin(trackParams[3]), std::tan(trackParams[4])) * pt;
                const pandora::CartesianVector newmomentum = m_rotation.MakeRotation(momentum);

                LOG_INFO("TrackCreator::CreateTracks()")
                << "Creating Track " << pTrack
                << " starting at " << newstartPosition
                << " ending at " << newendPosition
                << " with pt " << pt
                << " with momentum " << newmomentum;

                trackParameters.m_momentumAtDca = newmomentum;
                trackParameters.m_trackStateAtStart = pandora::TrackState(newstartPosition, momentum);
                trackParameters.m_trackStateAtEnd = pandora::TrackState(newendPosition, momentum);
                trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(newcalorimeterPosition, momentum);
                trackParameters.m_timeAtCalorimeter = 0.f;

                trackParameters.m_isProjectedToEndCap   = false;
                trackParameters.m_reachesCalorimeter    = true;
                trackParameters.m_canFormPfo            = true;
                trackParameters.m_canFormClusterlessPfo = true;

                try
                {
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters));
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
