////////////////////////////////////////////////////////////////////////
// Class:       TPCECALAssociation
// Plugin Type: producer (art v2_11_02)
// File:        TPCECALAssociation_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "RecoAlg/TrackPropagator.h"

#include "art_root_io/TFileService.h"
#include "TH2D.h"



namespace gar {
    namespace rec {



        class TPCECALAssociation : public art::EDProducer {
        public:
            explicit TPCECALAssociation(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            TPCECALAssociation(TPCECALAssociation const &) = delete;
            TPCECALAssociation(TPCECALAssociation &&) = delete;
            TPCECALAssociation & operator = (TPCECALAssociation const &) = delete;
            TPCECALAssociation & operator = (TPCECALAssociation &&) = delete;

            // Required functions.
            void beginJob() override;
            void produce(art::Event & e) override;

        private:

            // Declare member data here.
            std::string fTrackLabel;    ///< label to find the reco tracks
            std::string fClusterLabel;  ///< label to find the right reco caloclusters
            int fVerbosity;
            float fTrackEndXCut;        ///< Extrapolate only track ends outside central drift
            float fTrackEndRCut;        ///< Extrapolate only track ends outside central region
            float fPerpRCut;            ///< Max dist cluster center to circle of track (z,y) only
            float fBarrelXCut;          ///< Max dist cluster center (drift time compensated) track in x
            float fEndcapRphiCut;       ///< Max dist cluster center (drift time compensated) track in r*phi

            const detinfo::DetectorProperties*  fDetProp;    ///< detector properties
            const detinfo::DetectorClocks*      fClocks;     ///< Detector clock information
            const geo::GeometryCore*            fGeo;        ///< pointer to the geometry
            float                               maxXdisplacement;

            TH1F* radClusTrack;                              ///< Cluster to track in transverse plane
            TH2F* xClusTrack;                                ///< Cluster to track in barrel, x only, vs x of trackend
            TH1F* phiClusTrack;                              ///< Cluster to track in endcap
            TH1F* rPhiClusTrack;
        };



        TPCECALAssociation::TPCECALAssociation(fhicl::ParameterSet const & p) : EDProducer{p} {

            fTrackLabel    = p.get<std::string>("TrackLabel", "track");
            fClusterLabel  = p.get<std::string>("ClusterLabel","calocluster");
            fVerbosity     = p.get<int>("Verbosity", 0);
            fTrackEndXCut  = p.get<float>("TrackEndXCut",     240.0);
            fTrackEndRCut  = p.get<float>("TrackEndRCut",     238.0);
            fPerpRCut       = p.get<float>("fPerpRCut",        10.0);
            fBarrelXCut    = p.get<float>("fBarrelXCut",      20.0);
            fEndcapRphiCut = p.get<float>("fEndXCut",          32.0);
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fClocks  = gar::providerFrom<detinfo::DetectorClocksService>();
            fGeo     = gar::providerFrom<geo::Geometry>();

            produces< art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd > >();
        }



        void TPCECALAssociation::beginJob() {
            if (fVerbosity>0) {
                art::ServiceHandle< art::TFileService > tfs;
                radClusTrack = tfs->make<TH1F>("radClusTrack",
                    "(z,y) distance from cluster to track", 100, 0.0,60.0);
                xClusTrack   = tfs->make<TH2F>("xClusTrack",
                    "x distance from cluster to track vs x position of track end in barrel",
                    100,-250.0,+250.0, 100,-200.0,+200.0);
                phiClusTrack = tfs->make<TH1F>("phiClusTrack",
                    "angular distance from cluster to track in endcap", 90,-M_PI,+M_PI);
                rPhiClusTrack = tfs->make<TH1F>("rPhiClusTrack",
                    "distance from cluster to track in endcap", 100,-200.0,+200.0);
            }
        }



        void TPCECALAssociation::produce(art::Event & e) {
        
            // Get tracks and clusters.  If either is missing, just skip this event
            // processing.  That's not an exception
            art::Handle< std::vector<gar::rec::Cluster> > theClusterHandle;
            e.getByLabel(fClusterLabel, theClusterHandle);
            if (!theClusterHandle.isValid()) return;
            art::Handle< std::vector<gar::rec::Track> >   theTrackHandle;
            e.getByLabel(fTrackLabel, theTrackHandle);
            if (!theTrackHandle.isValid()) return;

            // fClocks service provides this info on event-by-event basis not at beginJob.
            maxXdisplacement = fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature())
                              *fClocks->SpillLength();



            std::unique_ptr<art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd>>
                            ClusterTrackAssns(new art::Assns<gar::rec::Cluster,gar::rec::Track, gar::rec::TrackEnd>);

            for (size_t iCluster=0; iCluster<theClusterHandle->size(); ++iCluster) {
                gar::rec::Cluster cluster = (*theClusterHandle)[iCluster];
                float yClus = cluster.Position()[1];    float zClus = cluster.Position()[2];
                float rClus = std::hypot(zClus,yClus);  float xClus = cluster.Position()[0];

                for (size_t iTrack=0; iTrack<theTrackHandle->size(); ++iTrack) {
                    gar::rec::Track track = (*theTrackHandle)[iTrack];

                    // Which if any ends of the track are near the outer edges of the TPC?
                    // These specific cuts could perhaps be increased
                    bool outside[2]; outside[0] = outside[1] = false;
                    
                    if ( abs(track.Vertex()[0]) > fTrackEndXCut ||
                        std::hypot(track.Vertex()[1],track.Vertex()[2]) > fTrackEndRCut ) {
                        outside[TrackEndBeg] = true;
                    }
                    if ( abs(track.End()[0]) > fTrackEndXCut ||
                        std::hypot(track.End()[1],track.End()[2]) > fTrackEndRCut ) {
                        outside[TrackEndEnd] = true;
                    }

                    // Plot and cut on the distance of cluster to helix in transverse plane.
                    for (TrackEnd tEnd = TrackEndBeg; tEnd >= TrackEndEnd; --tEnd) {
                        if (outside[tEnd]) {
                            float trackPar[5];    float trackEnd[3];
                            if (tEnd==TrackEndBeg) {
                                for (int i=0; i<5; ++i) trackPar[i] = track.TrackParBeg()[i];
                                for (int i=0; i<3; ++i) trackEnd[i] = track.Vertex()[i];
                            } else {
                                for (int i=0; i<5; ++i) trackPar[i] = track.TrackParEnd()[i];
                                for (int i=0; i<3; ++i) trackEnd[i] = track.End()[i];
                            }
                            float radius = 1.0/trackPar[2];
                            float zCent = trackPar[1] - radius*sin(trackPar[3]);                        
                            float yCent = trackPar[0] + radius*cos(trackPar[3]);

                            float distRadially = std::hypot(zClus-zCent,yClus-yCent) -abs(radius);
                            distRadially = abs(distRadially);
                            radClusTrack->Fill( distRadially );
                            
                            if ( distRadially > fPerpRCut ) {
                                // This track-cluster match fails
                                outside[tEnd] = false;
                                continue;
                            }

                            // Require plausible extrapolation in x as well
                            TVector3 clusterCenter(cluster.Position());
                            float retXYZ[3];
                            if (fGeo->PointInECALBarrel(clusterCenter)) {
                                // Extrapolate track to that radius
                                int errcode = util::TrackPropagator::PropagateToCylinder(
                                    trackPar,trackEnd,rClus, 0.0, 0.0, retXYZ);
                                if ( errcode!=0 ) {
                                    // This track-cluster match fails.  Error code 1, there is
                                    // no intersection at all, is possible
                                    outside[tEnd] = false;
                                    continue;
                                }
                                float extrapXerr = retXYZ[0] -clusterCenter.X();
                                xClusTrack->Fill(trackEnd[0], extrapXerr);
                                // extrapXerr is roughly maxXdisplacement/2 for trackend at
                                // x < -25cm and -maxXdisplacement/2 for trackend at x > 25cm.
                                float expected_mean = 0;
                                if (trackEnd[0]<-25) expected_mean = +maxXdisplacement/2.0;
                                if (trackEnd[0]>+25) expected_mean = -maxXdisplacement/2.0;
                                float cutQuantity = extrapXerr -expected_mean;
                                if ( abs(cutQuantity) > maxXdisplacement+fBarrelXCut ) {
                                    outside[tEnd] = false;
                                    continue;
                                }
                            } else {
                                // In an endcap.  How many radians in a maxXdisplacement?
                                float radiansInDrift = trackPar[2]*maxXdisplacement
                                                      / tan(trackPar[4]);
                                if ( abs(radiansInDrift) >= 2.0*M_PI ) {
                                    // Drat!  No distinguishing power here.
                                    outside[tEnd] = false;
                                    continue;
                                }
                                int errcode = util::TrackPropagator::PropagateToX(
                                    trackPar,trackEnd, xClus, retXYZ);
                                if ( errcode!=0 ) {
                                    // This track-cluster match fails.
                                    outside[tEnd] = false;
                                    continue;
                                }
                                // Find how many radians the closest point on the propagated
                                // track is from the cluster along the helix, projected onto
                                // to the perpendicular plane.  (bound to -pi,+pi range)
                                float angClus  = std::atan2(clusterCenter.Y()-yCent,clusterCenter.Z()-zCent);
                                float angXtrap = std::atan2(       retXYZ[1] -yCent,       retXYZ[2] -zCent);
                                float extrapPhiErr = angXtrap -angClus;
                                if (extrapPhiErr > +M_PI) extrapPhiErr -= 2.0*M_PI;
                                if (extrapPhiErr < -M_PI) extrapPhiErr += 2.0*M_PI;
                                phiClusTrack->Fill(extrapPhiErr);
                                float extrapRphiErr = abs(radius)*extrapPhiErr;
                                rPhiClusTrack->Fill( extrapRphiErr );
                                if (abs(extrapRphiErr) > fEndcapRphiCut) {
                                    outside[tEnd] = false;
                                    continue;
                                }
                            }    
                        }
                    }
                }
            }

            e.put(std::move(ClusterTrackAssns));
            return;
        }



        DEFINE_ART_MODULE(TPCECALAssociation)

    } // namespace rec
} // namespace gar
