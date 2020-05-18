#include "TrackCreator.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace gar {
  namespace gar_pandora {

    TrackCreator::TrackCreator(const Settings &settings , const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora)
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    TrackCreator::~TrackCreator()
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::CreateTracks(const art::Event *const pEvent)
    {
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

  }
}