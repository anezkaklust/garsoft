#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Api/PandoraApi.h"

#include "ReconstructionDataProducts/Track.h"

namespace gar {
  namespace gar_pandora {

    typedef std::vector<std::string> StringVector;
    typedef std::vector<gar::rec::Track *> TrackVector;

    class TrackCreator
    {
    public:

      class Settings
      {
      public:

        StringVector    m_trackCollections;                     ///< The reconstructed track collections

        int             m_shouldFormTrackRelationships;         ///< Whether to form pandora track relationships using v0 and kink info

        int             m_minTrackHits;                         ///< Track quality cut: the minimum number of track hits
        int             m_minFtdTrackHits;                      ///< Track quality cut: the minimum number of FTD track hits for FTD only tracks
        int             m_maxTrackHits;                         ///< Track quality cut: the maximum number of track hits

        float           m_d0TrackCut;                           ///< Track d0 cut used to determine whether track can be used to form pfo
        float           m_z0TrackCut;                           ///< Track z0 cut used to determine whether track can be used to form pfo

        int             m_usingNonVertexTracks;                 ///< Whether can form pfos from tracks that don't start at vertex
        int             m_usingUnmatchedNonVertexTracks;        ///< Whether can form pfos from unmatched tracks that don't start at vertex
        int             m_usingUnmatchedVertexTracks;           ///< Whether can form pfos from unmatched tracks that start at vertex
        float           m_unmatchedVertexTrackMaxEnergy;        ///< Maximum energy for unmatched vertex track

        float           m_d0UnmatchedVertexTrackCut;            ///< d0 cut used to determine whether unmatched vertex track can form pfo
        float           m_z0UnmatchedVertexTrackCut;            ///< z0 cut used to determine whether unmatched vertex track can form pfo
        float           m_zCutForNonVertexTracks;               ///< Non vtx track z cut to determine whether track can be used to form pfo

        int             m_reachesECalNBarrelTrackerHits;                  ///< Minimum number of barrel tracker hits to consider track as reaching ecal
        int             m_reachesECalNFtdHits;                  ///< Minimum number of ftd hits to consider track as reaching ecal
        float           m_reachesECalBarrelTrackerOuterDistance;          ///< Max distance from track to barrel tracker r max to id whether track reaches ecal
        int             m_reachesECalMinFtdLayer;               ///< Min layer in Ftd for tracks to be considered to have reached decal
        float           m_reachesECalBarrelTrackerZMaxDistance;           ///< Max distance from track to barrel tracker z max to id whether track reaches ecal
        float           m_reachesECalFtdZMaxDistance;           ///< Max distance from track hit to ftd z position to identify ftd hits
        float           m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum

        float           m_minTrackECalDistanceFromIp;           ///< Sanity check on separation between ip and track projected ecal position
        float           m_maxTrackSigmaPOverP;                  ///< Track fraction momentum error cut
        float           m_minMomentumForTrackHitChecks;         ///< Min track momentum required to perform final quality checks on number of hits

        float           m_maxBarrelTrackerInnerRDistance;                 ///< Track cut on distance from barrel tracker inner r to id whether track can form pfo
        float           m_minBarrelTrackerHitFractionOfExpected;          ///< Minimum fraction of TPC hits compared to expected
        int             m_minFtdHitsForBarrelTrackerHitFraction;          ///< Minimum number of FTD hits to ignore TPC hit fraction
        float           m_trackStateTolerance; ///< distance below tracker ecal radius the second trackstate in the ecal endcap is still passed to pandora
        std::string     m_trackingSystemName;  ///< name of the tracking system used for getting new track states

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