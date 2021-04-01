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

#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"
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
            int   fVerbosity;
            float fTrackEndXCut;        ///< Extrapolate only track ends outside central drift
            float fTrackEndRCut;        ///< Extrapolate only track ends outside central region
            float fPerpRCut;            ///< Max dist cluster center to circle of track (z,y) only
            float fBarrelXCut;          ///< Max dist cluster center (drift time compensated) track in x
            float fEndcapRphiCut;       ///< Max dist cluster center (drift time compensated) track in r*phi

            const detinfo::DetectorProperties*  fDetProp;    ///< detector properties
            const detinfo::DetectorClocks*      fClocks;     ///< detector clock information
            const geo::GeometryCore*            fGeo;        ///< pointer to the geometry
            float                               maxXdisplacement;

            // Position of TPC from geometry service; 1 S Boston Ave.
            float ItsInTulsa[3];

            TH1F* radClusTrack;                              ///< Cluster to track in transverse plane
            TH2F* xClusTrack;                                ///< Cluster to track in barrel, x only, vs x of trackend
            TH1F* phiClusTrack;                              ///< Cluster to track in endcap
            TH1F* rPhiClusTrack;
            TH1F* nEntryVsCut;                               ///< Number passing each cut
        };



        TPCECALAssociation::TPCECALAssociation(fhicl::ParameterSet const & p) : EDProducer{p} {

            fTrackLabel    = p.get<std::string>("TrackLabel", "track");
            fClusterLabel  = p.get<std::string>("ClusterLabel","calocluster");
            fVerbosity     = p.get<int>("Verbosity", 0);
            // Needs to be computed in produce; so will be 0.0 to do that 
            // otherwise of course the fcl file value appears
            fTrackEndXCut  = p.get<float>("TrackEndXCut",     215.0);
            fTrackEndRCut  = p.get<float>("TrackEndRCut",     230.0);
            fPerpRCut      = p.get<float>("PerpRCut",         10.0);
            fBarrelXCut    = p.get<float>("BarrelXCut",       35.0);
            fEndcapRphiCut = p.get<float>("EndXCut",          32.0);
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fClocks  = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();
            fGeo     = gar::providerFrom<geo::GeometryGAr>();

            produces< art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd > >();
        }



        void TPCECALAssociation::beginJob() {

            ItsInTulsa[0] = fGeo->TPCXCent();
            ItsInTulsa[1] = fGeo->TPCYCent();
            ItsInTulsa[2] = fGeo->TPCZCent();

            if (fVerbosity>0) {
                art::ServiceHandle< art::TFileService > tfs;
                radClusTrack  = tfs->make<TH1F>("radClusTrack",
                    "(z,y) distance from cluster to track", 100, 0.0,60.0);
                xClusTrack    = tfs->make<TH2F>("xClusTrack",
                    "x of track minus x of cluster vs x position of track end in barrel",
                    100,-250.0,+250.0, 100,-200.0,+200.0);
                phiClusTrack  = tfs->make<TH1F>("phiClusTrack",
                    "angular distance from cluster to track in endcap", 90,-M_PI,+M_PI);
                rPhiClusTrack = tfs->make<TH1F>("rPhiClusTrack",
                    "distance from cluster to track in endcap", 100,-200.0,+200.0);
                nEntryVsCut   = tfs->make<TH1F>("nEntryVsCut",
                    "Number track end & cluster combos vs sequentially applied cut", 30, -0.5, 29.5);
            }
        }



        void TPCECALAssociation::produce(art::Event & e) {
        
            // Get tracks and clusters.  If either is missing, just skip this event
            // processing.  That's not an exception
            art::Handle< std::vector<gar::rec::Cluster> > ClusterHandle;
            e.getByLabel(fClusterLabel, ClusterHandle);
            if (!ClusterHandle.isValid()) return;
            art::Handle< std::vector<gar::rec::Track> >   TrackHandle;
            e.getByLabel(fTrackLabel, TrackHandle);
            if (!TrackHandle.isValid()) return;

            // fClocks service provides this info on event-by-event basis not at beginJob.
            maxXdisplacement = fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature())
                              *fClocks->SpillLength();

            std::unique_ptr<art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd>>
                            ClusterTrackAssns(new art::Assns<gar::rec::Cluster,gar::rec::Track, gar::rec::TrackEnd>);
            auto const clusterPtrMaker = art::PtrMaker<rec::Cluster>(e, ClusterHandle.id());
            auto const   trackPtrMaker = art::PtrMaker<rec::Track>  (e, TrackHandle.id());



            for (size_t iCluster=0; iCluster<ClusterHandle->size(); ++iCluster) {
                gar::rec::Cluster cluster = (*ClusterHandle)[iCluster];
                TVector3 clusterCenter(cluster.Position());
                bool inECALBarrel = fGeo->PointInECALBarrel(clusterCenter);
                // fGeo uses one coordinate system, this code uses another.
                clusterCenter -=ItsInTulsa;
                float yClus = clusterCenter[1];
                float zClus = clusterCenter[2];
                float rClus = std::hypot(zClus,yClus);  
                float xClus = clusterCenter[0];

                for (size_t iTrack=0; iTrack<TrackHandle->size(); ++iTrack) {
                    gar::rec::Track track = (*TrackHandle)[iTrack];

                    // Which if any ends of the track are near the outer edges of the TPC?
                    // These specific cuts could perhaps be increased
                    bool outside[2]; outside[0] = outside[1] = false;

                    if ( abs(track.Vertex()[0] -ItsInTulsa[0]) > fTrackEndXCut ||
                        std::hypot(track.Vertex()[1] -ItsInTulsa[1],
                                   track.Vertex()[2] -ItsInTulsa[2]) > fTrackEndRCut ) {
                        outside[TrackEndBeg] = true;
                    }
                    if ( abs(track.End()[0] -ItsInTulsa[0])    > fTrackEndXCut ||
                        std::hypot(track.End()[1] -ItsInTulsa[1],
                                   track.End()[2] -ItsInTulsa[2])    > fTrackEndRCut ) {
                        outside[TrackEndEnd] = true;
                    }

                    // Plot and cut on the distance of cluster to helix in transverse plane.
                    for (TrackEnd iEnd = TrackEndBeg; iEnd >= TrackEndEnd; --iEnd) {
                        if (fVerbosity>0) nEntryVsCut->Fill(0);

                        if (outside[iEnd]) {
                            if (fVerbosity>0) nEntryVsCut->Fill(1);

                            float trackPar[5];    float trackEnd[3];
                            if (iEnd==TrackEndBeg) {
                                for (int i=0; i<5; ++i) trackPar[i] = track.TrackParBeg()[i];
                                for (int i=0; i<3; ++i) trackEnd[i] = track.Vertex()[i];
                            } else {
                                for (int i=0; i<5; ++i) trackPar[i] = track.TrackParEnd()[i];
                                for (int i=0; i<3; ++i) trackEnd[i] = track.End()[i];
                            }
                            // Translate to the center-of-MPD coordinate system
                            trackPar[0] -=ItsInTulsa[1];        trackPar[1] -=ItsInTulsa[2];
                            for (int i=0; i<3; ++i) trackEnd[i] -= ItsInTulsa[i];
                            
                            float radius = 1.0/trackPar[2];
                            float zCent = trackPar[1] - radius*sin(trackPar[3]);
                            float yCent = trackPar[0] + radius*cos(trackPar[3]);

                            float distRadially = std::hypot(zClus-zCent,yClus-yCent) -abs(radius);
                            distRadially = abs(distRadially);
                            if (fVerbosity>0) radClusTrack->Fill( distRadially );
                            
                            if ( distRadially > fPerpRCut ) {
                                // This track-cluster match fails
                                continue;
                            }
                            if (fVerbosity>0) nEntryVsCut->Fill(2);

                            // Require plausible extrapolation in x as well
                            float retXYZ1[3];    float retXYZ2[3];
                            if (inECALBarrel) {
                                if (fVerbosity>0) nEntryVsCut->Fill(10);

                                // Extrapolate track to that radius.  Using MPD center coords.
                                int errcode = util::TrackPropagator::PropagateToCylinder(
                                    trackPar,trackEnd,rClus, 0.0, 0.0, retXYZ1,retXYZ2);
                                if ( errcode!=0 ) {
                                    // This track-cluster match fails.  Error code 1, there is
                                    // no intersection at all, is possible
                                    continue;
                                }
                                if (fVerbosity>0) nEntryVsCut->Fill(11);
                                float extrapXerr;
                                float transDist1 = std::hypot(retXYZ1[2]-zClus,retXYZ1[1]-yClus);
                                float transDist2 = std::hypot(retXYZ2[2]-zClus,retXYZ2[1]-yClus);
                                bool looneyExtrap;
                                if (transDist1<transDist2) {
                                    TVector3 trackXYZ1(retXYZ1);
                                    trackXYZ1 += TVector3(ItsInTulsa);
                                    looneyExtrap =  !fGeo->PointInECALBarrel(trackXYZ1)
                                                 && !fGeo->PointInECALEndcap(trackXYZ1);
                                    extrapXerr = retXYZ1[0] -clusterCenter.X();
                                } else {
                                    TVector3 trackXYZ2(retXYZ2);
                                    trackXYZ2 += TVector3(ItsInTulsa);
                                    looneyExtrap =  !fGeo->PointInECALBarrel(trackXYZ2)
                                                 && !fGeo->PointInECALEndcap(trackXYZ2);
                                    extrapXerr = retXYZ2[0] -clusterCenter.X();
                                }
                                if (looneyExtrap) continue;
                                if (fVerbosity>0) nEntryVsCut->Fill(12);
                                if (fVerbosity>0) xClusTrack->Fill(trackEnd[0], extrapXerr);
                                // extrapXerr is roughly maxXdisplacement/2 for trackend at
                                // x < -25cm and -maxXdisplacement/2 for trackend at x > 25cm.
                                float expected_mean = 0;
                                if (trackEnd[0]<-25) expected_mean = +maxXdisplacement/2.0;
                                if (trackEnd[0]>+25) expected_mean = -maxXdisplacement/2.0;
                                float cutQuantity = extrapXerr -expected_mean;
                                if ( abs(cutQuantity) > maxXdisplacement/2.0+fBarrelXCut ) {
                                    continue;
                                }
                                if (fVerbosity>0) nEntryVsCut->Fill(13);

                            } else {
                                if (fVerbosity>0) nEntryVsCut->Fill(20);

                                // In an endcap.  How many radians in a maxXdisplacement?
                                float radiansInDrift = trackPar[2]*maxXdisplacement
                                                      / tan(trackPar[4]);
                                if ( abs(radiansInDrift) >= 2.0*M_PI ) {
                                    // Drat!  No distinguishing power here.
                                    continue;
                                }
                                if (fVerbosity>0) nEntryVsCut->Fill(21);

                                int errcode = util::TrackPropagator::PropagateToX(
                                    trackPar,trackEnd, xClus, retXYZ1);
                                if ( errcode!=0 ) {
                                    // This track-cluster match fails.
                                    continue;
                                }
                                if (fVerbosity>0) nEntryVsCut->Fill(22);
                                TVector3 trackXYZ(retXYZ1);
                                trackXYZ += TVector3(ItsInTulsa);
                                bool looneyExtrap =  !fGeo->PointInECALBarrel(trackXYZ)
                                                  && !fGeo->PointInECALEndcap(trackXYZ);
                                if (looneyExtrap) continue;
                                if (fVerbosity>0) nEntryVsCut->Fill(23);

                                // Find how many radians the closest point on the propagated
                                // track is from the cluster along the helix, projected onto
                                // to the perpendicular plane.  (bound to -pi,+pi range)
                                float angClus  = std::atan2(clusterCenter.Y()-yCent,clusterCenter.Z()-zCent);
                                float angXtrap = std::atan2(retXYZ1[1] -yCent,retXYZ1[2] -zCent) -angClus;
                                // angXtrap can indeed be outside of -PI to +PI
                                if (angXtrap > +M_PI) angXtrap -= 2.0*M_PI;
                                if (angXtrap < -M_PI) angXtrap += 2.0*M_PI;
                                if (fVerbosity>0) phiClusTrack->Fill(angXtrap);
                                float extrapRphiErr = abs(radius)*angXtrap;
                                if (fVerbosity>0) rPhiClusTrack->Fill( extrapRphiErr );
                                if (abs(extrapRphiErr) > fEndcapRphiCut) {
                                    continue;
                                }
                                if (fVerbosity>0) nEntryVsCut->Fill(24);
                            }
                            // Make the cluster-track association
                            art::Ptr<gar::rec::Cluster> const clusterPtr = clusterPtrMaker(iCluster);
                            art::Ptr<gar::rec::Track>   const   trackPtr = trackPtrMaker(iTrack);
                            ClusterTrackAssns->addSingle(clusterPtr,trackPtr,iEnd);
                        }
                    } // end loop over 2 ends of track
                }
            } // end loop over clusters

            e.put(std::move(ClusterTrackAssns));
            return;
        }



        DEFINE_ART_MODULE(TPCECALAssociation)

    } // namespace rec
} // namespace gar
