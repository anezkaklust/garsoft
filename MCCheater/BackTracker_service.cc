////////////////////////////////////////////////////////////////////////
//
//
// \file: BackTracker_service.cc
//
// brebel@fnal.gov
// modified Feb - Mar 2020 Leo Bellantoni, bellanto@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// GArSoft includes
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nug4/ParticleNavigation/EmEveIdCalculator.h"

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"





namespace gar {
  namespace cheat {
    
    //----------------------------------------------------------------------
    BackTracker::BackTracker(fhicl::ParameterSet const& pset,
                             art::ActivityRegistry    & reg)
      : BackTrackerCore(pset) {
      reg.sPreProcessEvent.watch(this, &BackTracker::Rebuild);
      fHasMC = fHasHits = fHasCalHits = fHasTracks = fHasClusters = false;
      fSTFU = 0;
    }

    //----------------------------------------------------------------------
    BackTracker::~BackTracker() {}

    //----------------------------------------------------------------------
    void BackTracker::beginJob() {}



    //----------------------------------------------------------------------
    void BackTracker::Rebuild(art::Event const& evt, art::ScheduleContext) {
      // Just drop the ScheduleContext that art requires. 
      RebuildNoSC(evt);
    }

    void BackTracker::RebuildNoSC(art::Event const& evt) {
      // do nothing if this is data
      if (evt.isRealData()) return;
      
      // do nothing if we are asked to do nothing
      if (fDisableRebuild) return;



      // We have a new event, so clear the collections from the previous events
      // Just the MC collections now, others as needed
      fParticleList.clear();
      fMCTruthList .clear();
      fTrackIDToMCTruthIndex.clear();



      // get the MCParticles from the event

      try
        {
          auto partCol = evt.getValidHandle<std::vector<simb::MCParticle> >(fG4ModuleLabel);
          if (partCol.failedToGet()) {
            ++fSTFU;
            if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
                << "failed to get handle to simb::MCParticle from " << fG4ModuleLabel
                << ", return";
            }
            return;
          }

          // get mapping from MCParticles to MCTruths; loop to fill fTrackIDToMCTruthIndex
          art::FindOneP<simb::MCTruth> fo(partCol, evt, fG4ModuleLabel);
          if( !fo.isValid() ){
            ++fSTFU;
            if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
                << "failed to associate MCTruth to MCParticle with " << fG4ModuleLabel
                << "; all attempts to backtrack will fail\n";
            }
            return;        

          } else {
            for (size_t iPart = 0; iPart < partCol->size(); ++iPart) {

              simb::MCParticle *part = new simb::MCParticle(partCol->at(iPart));
              fParticleList.Add(part);
          
              // get the simb::MCTruth associated to this MCParticle
              try {
                art::Ptr<simb::MCTruth> mct = fo.at(iPart);
                if (fMCTruthList.size() < 1) {
                  fMCTruthList.push_back(mct);
                } else {
                  // check that we are not adding a simb::MCTruth twice to 
                  // the collection we know that all the particles for a 
                  // given simb::MCTruth are put into the collection of
                  // particles at the same time, so we can just check that the
                  // current ::art::Ptr has a different id than the last one put
                  if (!(mct == fMCTruthList.back())) fMCTruthList.push_back(mct);
                }

                // Fill the track id to mctruth index map.  partCol will not contain
                // negative track ids.
                fTrackIDToMCTruthIndex[partCol->at(iPart).TrackId()] = fMCTruthList.size() - 1;
              }
              catch (cet::exception &ex) {
                ++fSTFU;
                if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
                    << "unable to find MCTruth from ParticleList created in "
                    << fG4ModuleLabel
                    << "; all attempts to backtrack will fail\n"
                    << "message from caught exception:\n" << ex;
                }
              }
            } // end loop over MCParticles
          }

        }


      catch (cet::exception &ex)
        {
          return;
        }

      // Register the electromagnetic Eve calculator
      fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);

      fHasMC = true;



      // Inverse of drift velocity and longitudinal diffusion constant will
      // be needed by the backtracker when fHasHits, and also the geometry
      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      fInverseVelocity = 1.0/detProp->DriftVelocity(detProp->Efield(),
                                                    detProp->Temperature());

      auto garProp = gar::providerFrom<detinfo::GArPropertiesService>();
      fLongDiffConst = std::sqrt(2.0 * garProp->LongitudinalDiffusion());

      fGeo = gar::providerFrom<geo::Geometry>();

      // Create the channel to energy deposit collection.  Start by 
      // looking for the RawDigit collection from the event and create
      // a FindMany mapping between the digits and energy deposits if it exists.
      for (auto vec : fChannelToEDepCol) vec.clear();
      fChannelToEDepCol.clear();
      art::Handle<std::vector<raw::RawDigit>> digCol;
      evt.getByLabel(fRawTPCDataLabel, digCol);
      if (!digCol.isValid()) {
        ++fSTFU;
        if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
          MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
            << "Unable to find RawDigits in " << fRawTPCDataLabel <<
            "; no backtracking in TPC will be possible";
        }

      } else {
        art::FindMany<sdp::EnergyDeposit> fmEnergyDep(digCol, evt, fRawTPCDataLabel);
        if (!fmEnergyDep.isValid()) {
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
              << "Unable to find valid association between RawDigits and "
              << "energy deposits in " << fRawTPCDataLabel <<
              "; no backtracking in TPC will be possible";
          }
 
        } else {
          fChannelToEDepCol.resize(fGeo->NChannels());
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_DEBUG("BackTracker_service::RebuildNoSC")
              << "There are " << fChannelToEDepCol.size()
              << " channels in the geometry";
          }

          std::vector<const sdp::EnergyDeposit*> eDeps;
          for (size_t iDig = 0; iDig < digCol->size(); ++iDig) {
            eDeps.clear();
            fmEnergyDep.get(iDig, eDeps);    // uses vector::insert
            if (eDeps.size() < 1) continue;
            fChannelToEDepCol[ (*digCol)[iDig].Channel() ].swap(eDeps);
          }
          fHasHits = true;
        }
      }



      // Create the CellID to energy deposit collection.  Start by 
      // looking for the CaloRawDigit collection from the event and create
      // a FindMany mapping between these digits and CaloDeposits if it exists.
      fECALTrackToTPCTrack->clear();
      fCellIDToEDepCol.clear();
      art::Handle<std::vector<raw::CaloRawDigit>> caloDigCol;
      evt.getByLabel(fRawCaloDataLabel, caloDigCol);
      if (!caloDigCol.isValid()) {
        ++fSTFU;
        if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
          MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
            << "Unable to find CaloRawDigits in " << fRawCaloDataLabel <<
            "; no backtracking in ECAL will be possible";
        }            

      } else {
        art::FindMany<sdp::CaloDeposit> fmCaloDep(caloDigCol, evt, fRawCaloDataLabel);
        if (!fmCaloDep.isValid()) {
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
              << "Unable to find valid association between CaloRawDigits and "
              << "calorimeter energy deposits in " << fRawCaloDataLabel <<
              "; no backtracking in ECAL will be possible";
          }

        } else {
          std::vector<const sdp::CaloDeposit*> CeDeps;
          for (size_t iCalDig = 0; iCalDig < caloDigCol->size(); ++iCalDig) {
            CeDeps.clear();
            fmCaloDep.get(iCalDig, CeDeps);    // uses vector::insert
            if (CeDeps.size() < 1) continue;
            fCellIDToEDepCol[ (*caloDigCol)[iCalDig].CellID() ] = CeDeps;
          }
          fHasCalHits = true;
        }
      }



      // Create an unordered map of IDNumbers of reco'd Tracks to the hits in those tracks
      // Just to keep the code comprehensible, first map Tracks to TPCClusters, then 
      // TPCClusters to Hits, then build the combined map.
      std::unordered_map< rec::IDNumber, std::vector<const rec::TPCCluster*> > lTrackIDToTPCClusters;
      art::Handle<std::vector<rec::Track>> recoTrackCol;
      evt.getByLabel(fTrackLabel, recoTrackCol);
      if (!recoTrackCol.isValid()) {
        ++fSTFU;
        if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
          MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
            << "Unable to find rec::Tracks in " << fTrackLabel <<
            "; no backtracking of reconstructed tracks will be possible";
        }        

      } else {
        art::FindMany<rec::TPCCluster> fmTPCClusts(recoTrackCol, evt, fTrackLabel);
        if (!fmTPCClusts.isValid()) {
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
              << "Unable to find valid association between Tracks and "
              << "TPCClusters in " << fTrackLabel <<
              "; no backtracking of reconstructed tracks will be possible";
          }

        } else {
          std::vector<const rec::TPCCluster*> tpclusThisTrack;
          for (size_t iTrack = 0; iTrack < recoTrackCol->size(); ++iTrack) {
            tpclusThisTrack.clear();
            fmTPCClusts.get(iTrack, tpclusThisTrack);    // uses vector::insert
            if (tpclusThisTrack.size() < 1) continue;
            lTrackIDToTPCClusters[ (*recoTrackCol)[iTrack].getIDNumber() ] = tpclusThisTrack;
          }
        }
      }

      // Construct map of TPCCluster to Hits
      std::unordered_map< rec::IDNumber,std::vector<art::Ptr<rec::Hit>> > lTPClusIDsToHits;
      art::Handle<std::vector<rec::TPCCluster>> tpclusCol;
      evt.getByLabel(fTPCClusterLabel, tpclusCol);
      if (!tpclusCol.isValid()) {
        ++fSTFU;
        if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
          MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
            << "Unable to find rec::TPCClusters in " << fTPCClusterLabel <<
            "; no backtracking of reconstructed tracks will be possible";
        }

      } else {
        art::FindManyP<rec::Hit> fmHits(tpclusCol, evt, fTPCClusterLabel);
        if (!fmHits.isValid()) {
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
              << "Unable to find valid association between TPCClusters and "
              << "Hits in " << fTPCClusterLabel <<
              "; no backtracking of reconstructed tracks will be possible";
          }

        } else {
          std::vector<art::Ptr<rec::Hit>> hitsThisTPCCluster;
          for (size_t iTPClus = 0; iTPClus < tpclusCol->size(); ++iTPClus) {
            hitsThisTPCCluster.clear();
            fmHits.get(iTPClus, hitsThisTPCCluster);    // uses vector::insert
            if (hitsThisTPCCluster.size() < 1) continue;
            lTPClusIDsToHits[(*tpclusCol)[iTPClus].getIDNumber()] = hitsThisTPCCluster;
          }
        }
      }

      // Now compound the lists to make fTrackIDToHits
      fTrackIDToHits.clear();
      if ( !lTrackIDToTPCClusters.empty() && !lTPClusIDsToHits.empty()) {
        for (auto aTrack : *recoTrackCol) {
          std::vector<const rec::TPCCluster*> tpclusThisTrack =
            lTrackIDToTPCClusters[ aTrack.getIDNumber() ];
          std::vector<art::Ptr<rec::Hit>> hitsThisTrack;
          for (const rec::TPCCluster* aTPClus : tpclusThisTrack) {
            std::vector<art::Ptr<rec::Hit>> hitsThisClus = 
              lTPClusIDsToHits[aTPClus->getIDNumber()];
            hitsThisTrack.insert( hitsThisTrack.end(),
                                  hitsThisClus.begin(),hitsThisClus.end() );
          }
          // Explicit assumption is that each hit is in at most 1 cluster so
          // there are no duplicates in hitsThisTrack
          fTrackIDToHits[aTrack.getIDNumber()] = hitsThisTrack;
        }
      }
      fHasTracks = true;



      // Create an unordered map of IDNumbers of reco'd Clusters to the CaloHits in those tracks
      art::Handle<std::vector<rec::Cluster>> recoClusterCol;
      evt.getByLabel(fClusterLabel, recoClusterCol);
      if (!recoClusterCol.isValid()) {
        ++fSTFU;
        if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
          MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
            << "Unable to find rec::Clusters in " << fClusterLabel <<
            "; no backtracking of reconstructed ECAL clusters will be possible";
        }

      } else {
        art::FindManyP<rec::CaloHit> fmCaloHits(recoClusterCol, evt, fClusterLabel);
        if (!fmCaloHits.isValid()) {
          ++fSTFU;
          if (fSTFU<=10) {    // Ye who comprehend messagelogger, doeth ye better.
            MF_LOG_WARNING("BackTracker_service::RebuildNoSC")
              << "Unable to find valid association between Clusters and "
              << "CaloHits in " << fTrackLabel <<
              "; no backtracking of reconstructed tracks will be possible";
          }

        } else {
          std::vector<art::Ptr<rec::CaloHit>> caloHitsThisCluster;
          for (size_t iCluster = 0; iCluster < recoClusterCol->size(); ++iCluster) {
            caloHitsThisCluster.clear();
            fmCaloHits.get(iCluster, caloHitsThisCluster);    // uses vector::insert
            if (caloHitsThisCluster.size() < 1) continue;
            fClusterIDToCaloHits[ (*recoClusterCol)[iCluster].getIDNumber() ] = caloHitsThisCluster;
          }
        }
      }
      fHasClusters = true;

      return;
    }





    DEFINE_ART_SERVICE(BackTracker)
  } // namespace cheat
} // namespace gar
