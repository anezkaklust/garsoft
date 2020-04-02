//
//  StripSplitterAlg.h
//
//  Created by Eldwan Brianne on 18.04.2019
//  Based on DDStripSplitter
//  https://github.com/iLCSoft/MarlinReco/blob/master/Clustering/hybridEcalSplitter/include/DDStripSplitter.h

#ifndef GAR_RECOALG_StripSplitterAlg_h
#define GAR_RECOALG_StripSplitterAlg_h

#include "art/Framework/Core/ModuleMacros.h"
#include "Geometry/GeometryCore.h"

#include "ReconstructionDataProducts/CaloHit.h"

#include "TVector3.h"

namespace fhicl{
    class ParameterSet;
}

namespace gar{
    namespace rec{
        namespace alg{

            class StripSplitterAlg {
            public:

                StripSplitterAlg(fhicl::ParameterSet const& pset);

                virtual ~StripSplitterAlg();

                void reconfigure(fhicl::ParameterSet const& pset);

                void ClearLists();

                void PrepareAlgo(const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

                void DoStripSplitting();

                bool GetSaveStripEndsFlag() const { return fSaveStripEnds; }

                std::vector<const gar::rec::CaloHit*> getStripEndsHits() const { return fStripEndsHits; }

                std::vector <const gar::rec::CaloHit*> getSplitHits() const { return splitStripHits; }

                std::vector <const gar::rec::CaloHit*> getUnSplitHits() const { return unSplitStripHits; }

            private:

                void getVirtualHits(const gar::rec::CaloHit *hit, int orientation, bool isBarrel, std::vector <const gar::rec::CaloHit*> &virtualhits);

                TVector3 stripIntersect(const gar::rec::CaloHit *hit0, const TVector3& dir0, const gar::rec::CaloHit *hit1, const TVector3& dir1);

                gar::geo::GeometryCore const* fGeo; ///< geometry information

                std::string fSSAAlgName;
                bool fSaveStripEnds;

                int fInnerSymmetry;
                double fStripWidth;
                double fStripLength;
                int fnVirtual;

                std::vector<const gar::rec::CaloHit*> m_CaloHitVecOdd;
                std::vector<const gar::rec::CaloHit*> m_CaloHitVecEven;

                std::vector<const gar::rec::CaloHit*> fStripEndsHits;
                std::vector <const gar::rec::CaloHit*> unSplitStripHits;
                std::vector <const gar::rec::CaloHit*> splitStripHits;

                enum StripOrientation {
                    TRANSVERSE = 0, LONGITUDINAL
                };
            };

        } // namespace alg
    } // namespace rec
} //namespace gar

#endif /* GAR_RECOALG_KNNClusterAlg_h */
