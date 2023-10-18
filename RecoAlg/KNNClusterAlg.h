//
//  KNNClusterAlg.h
//
//  Created by Eldwan Brianne on 17.04.2019
//  Based on MarlinReco simple clusetring algo
//  https://github.com/iLCSoft/MarlinReco/blob/master/Clustering/NNClustering/include/NNClusterProcessor.h
//
//  Modification:  Addition of usesTracks, Leo Bellantoni 21 Jun2019
//      If we have more clustering algorithms someday, this should be derived from
//      some class with a bunch of virtual methods or some sort of thing like that.


#ifndef GAR_RECOALG_KNNClusterAlg_h
#define GAR_RECOALG_KNNClusterAlg_h

#include "art/Framework/Core/ModuleMacros.h"
#include "Geometry/GeometryCore.h"

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/Track.h"

#include <list>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace fhicl{
  class ParameterSet;
}

namespace gar{
  namespace rec{
    namespace alg{

      typedef std::vector<gar::rec::CaloHit*> CaloHitVec;
      typedef std::vector<gar::rec::Track*> TrackVec;
      typedef std::vector<gar::rec::Cluster*> ClusterVec;

      class KNNClusterAlg {
      public:

        KNNClusterAlg(fhicl::ParameterSet const& pset);

        virtual ~KNNClusterAlg();

        void reconfigure(fhicl::ParameterSet const& pset);

        bool usesTracks() {return false;}

        void PrepareAlgo(const std::vector< art::Ptr<gar::rec::Track> > &trkVector, const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector, std::unordered_map< const gar::rec::Track*, art::Ptr<gar::rec::Track> > &trkMaptoArtPtr, std::unordered_map< const gar::rec::CaloHit*, art::Ptr<gar::rec::CaloHit> > &hitMaptoArtPtr);

        void DoClustering();

        ClusterVec GetFoundClusters() const { return clusterVector; }

      private:

        void ClearLists();

        gar::geo::GeometryCore const* fGeo;        ///< geometry information

        std::string fClusterAlgName;
        float m_EnergyCut;
        float m_DistanceCut;

        CaloHitVec m_CaloHitVec;
        TrackVec m_TrackVec;
        ClusterVec clusterVector;
      };

    } // namespace alg
  } // namespace rec
} //namespace gar

#endif /* GAR_RECOALG_KNNClusterAlg_h */
