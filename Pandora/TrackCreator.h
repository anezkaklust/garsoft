#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Api/PandoraApi.h"

#include "ReconstructionDataProducts/Track.h"

namespace gar {
  namespace gar_pandora {

    typedef std::vector<gar::rec::Track *> TrackVector;

    class TrackCreator
    {
    public:

      class Settings
      {
      public:

        std::string       m_trackCollection;                     ///< The reconstructed track collection

        unsigned int      m_minTrackHits;                         ///< Track quality cut: the minimum number of track hits
        unsigned int      m_maxTrackHits;                         ///< Track quality cut: the maximum number of track hits

        float             m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum

        float             m_maxTrackSigmaPOverP;                  ///< Track fraction momentum error cut
        float             m_minMomentumForTrackHitChecks;         ///< Min track momentum required to perform final quality checks on number of hits

        float             m_bField;                       ///< The bfield
        int               m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
        float             m_eCalBarrelInnerPhi0;          ///< ECal barrel inner phi 0
        float             m_eCalBarrelInnerR;             ///< ECal barrel inner radius
        float             m_eCalEndCapInnerZ;             ///< ECal endcap inner z
      };

      TrackCreator(const Settings &settings, const pandora::Pandora *const pPandora);

      virtual ~TrackCreator();

      pandora::StatusCode CreateTracks(const art::Event *const pEvent);

      void Reset();

    protected:

      const Settings          m_settings;

      const pandora::Pandora  &m_pandora;

      virtual bool PassesQualityCuts(const gar::rec::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const;

      virtual void TrackReachesECAL(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

      virtual void DefineTrackPfoUsage(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;
    };

    inline void TrackCreator::Reset()
    {

    }
  }
}

#endif // #ifndef TRACK_CREATOR_H
