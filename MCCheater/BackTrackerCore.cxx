//
//  BackTrackerCore.cxx
//  Created by Brian Rebel on 3/13/17.
//  Major rewrite, Leo Bellantoni Mar 2020
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "SimulationDataProducts/sim.h"
#include "CoreUtils/ServiceUtil.h"
#include "MCCheater/BackTrackerCore.h"
#include "Geometry/Geometry.h"

#include "TMath.h"

// For timing studies
//#include <chrono>
//using namespace std::chrono;


namespace gar{
    namespace cheat{
    
        //--------------------------------------------------------------------------
        BackTrackerCore::BackTrackerCore(fhicl::ParameterSet const& pset)
          : fClocks(nullptr), fGeo(nullptr) {

            fClocks = gar::providerFrom<detinfo::DetectorClocksService>();
            fGeo    = gar::providerFrom<geo::Geometry>();

            fG4ModuleLabel        = pset.get<std::string>("G4ModuleLabel",           "geant");
            fRawTPCDataLabel      = pset.get<std::string>("RawTPCDataLabel",         "daq");
            fRawCaloDataLabel     = pset.get<std::string>("RawCaloDataLabel",        "daqecal");
            fECALtimeResolution   = pset.get<double     >("ECALtimeResolution",       1.0);
            fMinHitEnergyFraction = pset.get<double     >("MinHitEnergyFraction",     0.1);
            fMinCaloHitEnergyFrac = pset.get<double     >("MinCaloHitEnergyFrac",     0.1);
            fDisableRebuild       = pset.get<bool       >("DisableRebuild",           false);

            fECALTrackToTPCTrack  = new std::unordered_map<int, int>;
            return;
        }

        //--------------------------------------------------------------------------
        BackTrackerCore::~BackTrackerCore() {}
    



    
        //----------------------------------------------------------------------
        simb::MCParticle* const BackTrackerCore::TrackIDToParticle(int const& id) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::TrackIDToParticle")
                    << "Attempting to backtrack without MC truth information";
            }

            // What to do for negative trackID?  Seems they are in fParticleList too.

            auto part_it = fParticleList.find(id);

            if (part_it == fParticleList.end()) {
                LOG_WARNING("BackTrackerCore::TrackIDToParticle")
                    << "can't find particle with track id " << id
                    << " in sim::ParticleList returning null pointer";
                    return nullptr;
            }

            return part_it->second;
        }



        //----------------------------------------------------------------------
        simb::MCParticle* const BackTrackerCore::FindEve(simb::MCParticle* const p) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindEve")
                    << "Attempting to backtrack without MC truth information";
            }

            int EveTrackId = fParticleList.EveId( p->TrackId() );
            return TrackIDToParticle(EveTrackId);
        }



        //----------------------------------------------------------------------
        simb::MCParticle* const BackTrackerCore::FindTPCEve(simb::MCParticle* const p) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindTPCEve")
                    << "Attempting to backtrack without MC truth information";
            }

            // If I had ever been here before I would probably know just what to do / Don't you?
            int dejaVu = p->TrackId();
            if ( fECALTrackToTPCTrack->find(dejaVu) != fECALTrackToTPCTrack->end() ) {
                return TrackIDToParticle( fECALTrackToTPCTrack->at(dejaVu) );
            }

            simb::MCParticle* walker  = p;
            if (walker==nullptr) return nullptr;
            while ( walker != nullptr &&
                !fGeo->PointInGArTPC( TVector3(walker->Vx(),walker->Vy(),walker->Vz()) ) ) {
                // Always walk up the ParticleList in case someday somebody changes it
                int mommaTID = fParticleList[walker->TrackId()]->Mother();
                if (mommaTID == 0) {
                    // you are at the top of the tree
                    break;
                }
                simb::MCParticle* momma = TrackIDToParticle( mommaTID );
                walker = momma;
            }
            (*fECALTrackToTPCTrack)[dejaVu] = walker->TrackId();
            return walker;
         }



        //----------------------------------------------------------------------
        simb::MCParticle* const BackTrackerCore::FindTPCEve(int trackID) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindTPCEve")
                    << "Attempting to backtrack without MC truth information";
            }

            // If I had ever been here before I would probably know just what to do / Don't you?
            if ( fECALTrackToTPCTrack->find(trackID) != fECALTrackToTPCTrack->end() ) {
                return TrackIDToParticle( fECALTrackToTPCTrack->at(trackID) );
            }

            simb::MCParticle* walker  = TrackIDToParticle(trackID);
            while ( walker != nullptr &&
                !fGeo->PointInGArTPC( TVector3(walker->Vx(),walker->Vy(),walker->Vz()) ) ) {
                // Always walk up the ParticleList in case someday somebody changes it
                int mommaTID = fParticleList[walker->TrackId()]->Mother();
                if (mommaTID == 0) {
                    // you are at the top of the tree
                    break;
                }
                simb::MCParticle* momma = TrackIDToParticle( mommaTID );
                walker = momma;
            }

           (*fECALTrackToTPCTrack)[trackID] = walker->TrackId();
            return walker;
         }



        //----------------------------------------------------------------------
        bool BackTrackerCore::IsDescendedFrom (simb::MCParticle* const forebear,
                                               simb::MCParticle* const afterbear) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::IsDescendedFrom")
                    << "Attempting to backtrack without MC truth information";
            }

            if (forebear==nullptr || afterbear==nullptr) return false;

            int stopper = forebear->TrackId();
            int walker  = afterbear->TrackId();
            while ( walker!=stopper ) {
                // Always walk up the ParticleList in case someday somebody changes it
                int momma = fParticleList[walker]->Mother();
                if (momma == 0) {
                    return false;
                } else {
                    walker = momma;
                }
            }
            return true;
         }



        //----------------------------------------------------------------------
        art::Ptr<simb::MCTruth> const
        BackTrackerCore::ParticleToMCTruth(simb::MCParticle* const p) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::ParticleToMCTruth")
                    << "Attempting to backtrack without MC truth information";
            }

            int id = p->TrackId();
            // find the entry in the MCTruth collection for this track id
            size_t mct = fTrackIDToMCTruthIndex.find(std::abs(id))->second;

            if( mct > fMCTruthList.size() ) {
                throw cet::exception("BackTrackerCore::ParticleToMCTruth")
                    << "attempting to find MCTruth index for out-of-range value: "
                    << mct << "/" << fMCTruthList.size();
            }
            return fMCTruthList[mct];
        }





        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::HitToHitIDEs(art::Ptr<gar::rec::Hit> const& hit) const {
            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitToHitIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            // By cutting on the StartTime and EndTime for this hit, we get only the
            // IDEs for this hit; ides.size() should basically be the number of 
            // MCParticles that contributed to this hit.
            std::vector<HitIDE> ides;
            ides = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());
            return ides;
        }

        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::HitToHitIDEs(gar::rec::Hit const& hit) const {
            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitToHitIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            // By cutting on the StartTime and EndTime for this hit, we get only the
            // IDEs for this hit; ides.size() should basically be the number of 
            // MCParticles that contributed to this hit.
            std::vector<HitIDE> ides;
            ides = this->ChannelToHitIDEs(hit.Channel(),hit.StartTime(),hit.EndTime());
            return ides;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<gar::rec::Hit>>
        BackTrackerCore::ParticleToHits(simb::MCParticle* const p,
                                        std::vector<art::Ptr<gar::rec::Hit>> const& allhits,
                                        bool checkNeutrals) const {
            // returns a subset of the hits in the allhits collection that match
            // the MC particle passed as an argument

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::ParticleToHits")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            std::vector<art::Ptr<gar::rec::Hit>> hitList;
            int PDGpid = p->PdgCode();
            bool obviousNeutral = PDGpid==22 || PDGpid==2112 || PDGpid==111 ||
                                  abs(PDGpid)==14 || abs(PDGpid)==12;
            if (!checkNeutrals && obviousNeutral) return hitList;

            int tkID = p->TrackId();
            std::vector<HitIDE> hids;
            
            for (auto hit : allhits) {
                hids.clear();

                hids = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());

                for (auto const& hid : hids) {
                    if ( hid.trackID==tkID && hid.energyFrac>fMinHitEnergyFraction ) {
                        hitList.push_back(hit);
                    }
                }
            }

            return hitList;
        }



        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::ChannelToHitIDEs(raw::Channel_t const& channel,
                                          double const start, double const stop) const {

            // loop over the energy deposits in the channel and grab those in time window

            std::vector<HitIDE> hitIDEs;

            if(fChannelToEDepCol.size() < channel){
                LOG_WARNING("BackTrackerCore::ChannelToHitIDEs")
                    << "Attempting to find energy deposits for channel " << channel
                    << " while there are only " << fChannelToEDepCol.size()
                    << " channels available, returning an empty vector.";
                return hitIDEs;
            }

            auto chanEDeps = fChannelToEDepCol[channel];

            if(chanEDeps.size() < 1){
                LOG_WARNING("BackTrackerCore::ChannelToHitIDEs")
                    << "No sdp::EnergyDeposits for selected channel: " << channel
                    << ". There is no way to backtrack the given hit, returning"
                    << " an empty vector.";
                return hitIDEs;
            }

            std::map<int,size_t> tidmap;  // Maps track IDs to index into hitIDEs

            // first get the total energy represented by all track ids for
            // this channel and range of tdc values
            unsigned short tdc = 0;
            float chanpos[3]={0.,0.,0};
            fGeo->ChannelToPosition(channel, chanpos);

            double totalEallTracks = 0.0;
            for (auto const& edep : chanEDeps) {

                if (edep->TrackID() == gar::sdp::NoParticleId) continue;

                tdc = fClocks->TPCG4Time2TDC 
                        (edep->Time() +TMath::Abs(edep->X()-chanpos[0]) *fInverseVelocity);
                if ( (unsigned int)start <= tdc && tdc <= (unsigned int)stop) {
                    float energy = edep->Energy();
                    totalEallTracks += energy;
                    int tid = edep->TrackID();
                    auto tidmapiter = tidmap.find(tid);
                    if (tidmapiter == tidmap.end()) {
                        hitIDEs.emplace_back(tid, 0.0, energy);
                        tidmap[tid] = hitIDEs.size()-1;
                    } else {
                        hitIDEs.at(tidmapiter->second).energyTot += energy;
                    }
                }
            }

            // Why was this here? 
            //if (totalE < 1.0e-5) totalE = 1.0;

            // loop over the hitIDEs to set the fractional energy for each TrackID
            for (size_t i = 0; i < hitIDEs.size(); ++i) {
                hitIDEs[i].energyFrac = hitIDEs[i].energyTot / totalEallTracks;
            }
            return hitIDEs;
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::HitPurity(simb::MCParticle* const p,
            std::vector<art::Ptr<gar::rec::Hit> > const& hits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitCollectionPurity")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }


            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            float total   = 0.0;
            for (auto hit : hits) {
                std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                weight = (weightByCharge) ? hit->Signal() : 1.0;
                total += weight;

                for (auto iHitIDE : hitIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT trackIDs
                    if (iHitIDE.trackID    == trackIDin &&
                        iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            double purity      = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                purity      = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(purity,uncertainty);
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::HitEfficiency(simb::MCParticle* const p,
            std::vector<art::Ptr<gar::rec::Hit> > const& hits,
            std::vector<art::Ptr<gar::rec::Hit> > const& allhits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitCollectionEfficiency")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            for (auto hit : hits) {
                std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                weight = (weightByCharge) ? hit->Signal() : 1.0;

                for (auto iHitIDE : hitIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT track IDs
                    if (iHitIDE.trackID    == trackIDin &&
                        iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            // Next figure out how many hits in the whole collection are associated with this id
            float total   = 0.0;
            if (desired>0) {
                for (auto hit : allhits) {
                    std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                    weight = (weightByCharge) ? hit->Signal() : 1.0;

                    for (auto iHitIDE : hitIDEs) {
                        // Don't double count if this hit has > 1 of the desired GEANT track IDs
                        if (iHitIDE.trackID    == trackIDin &&
                           iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                           total += weight;
                            break;
                        }
                    }
                } // end loop over all hits
            }

            double efficiency  = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                efficiency  = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(efficiency,uncertainty);
        }





        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CaloHitToCalIDEs(art::Ptr<gar::rec::CaloHit> const& hit) const {
            if (!fHasMC || !fHasCalHits) {
                throw cet::exception("BackTrackerCore::CaloHitToCalIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            // Cut on Time for this hit +/- ECALtimeResolution givenin BackTracker.fcl
            std::vector<CalIDE> ides;
            ides = this->CellIDToCalIDEs(hit->CellID(),hit->Time());
            return ides;
        }

        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CaloHitToCalIDEs(gar::rec::CaloHit const& hit) const {
            if (!fHasMC || !fHasCalHits) {
                throw cet::exception("BackTrackerCore::CaloHitToCalIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            // Cut on Time for this hit +/- ECALtimeResolution given in BackTracker.fcl
            std::vector<CalIDE> ides;
            ides = this->CellIDToCalIDEs(hit.CellID(),hit.Time());
            return ides;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<gar::rec::CaloHit>>
        BackTrackerCore::ParticleToCaloHits(simb::MCParticle* const p,
                                            std::vector<art::Ptr<gar::rec::CaloHit>> const& allhits) const {
            // returns a subset of the CaloHits in the allhits collection that match
            // the MC particle passed as an argument

            // Timing
            //time_point start(std::chrono::high_resolution_clock::now());
            //auto totalThere = start -start;

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::ParticleToCaloHits")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            std::vector<art::Ptr<gar::rec::CaloHit>> calHitList;
            int tkID = p->TrackId();
            std::vector<CalIDE> calhids;
            
            for (auto hit : allhits) {
                calhids.clear();

                // Timing
                //start   = high_resolution_clock::now();           
                calhids = this->CellIDToCalIDEs(hit->CellID(),hit->Time());
                // Timing
                //totalThere += high_resolution_clock::now() -start;     

                for (auto const& hid : calhids) {
                    if ( hid.trackID==tkID && hid.energyFrac>fMinCaloHitEnergyFrac ) {
                        calHitList.push_back(hit);
                   }
                }
            }
            // Timing
            //duration<double, std::micro> dA = totalThere;
            //std::cout << "Time in ChannelToHitIDES: " <<
            //    dA.count() << " for " << allhits.size() << " hits" << std::endl;
            return calHitList;
        }



        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CellIDToCalIDEs(raw::CellID_t const& cell,
                                         float const time) const {

            // loop over the energy deposits in the channel and grab those in time window

            std::vector<CalIDE> calIDEs;

            std::vector<const sdp::CaloDeposit*>  cellEDeps;
            try {
                cellEDeps = fCellIDToEDepCol.at(cell);
            }
            catch (std::out_of_range) {
                LOG_WARNING("BackTrackerCore::CellIDToCalIDEs")
                    << "No sdp::EnergyDeposits for selected ECAL cell: " << cell
                    << ". There is no way to backtrack the given calorimeter hit,"
                    << " returning an empty vector.";
                return calIDEs;
            }

            std::map<int,size_t> tidmap;  // Maps track IDs to index into calIDEs

            // Get the total energy produced by all track ids for this cell & time
            float start = time -fECALtimeResolution;
            float stop  = time +fECALtimeResolution;

            double totalEallTracks = 0.0;
            for (auto const& edep : cellEDeps) {
                if (edep->TrackID() == gar::sdp::NoParticleId) continue;

                if ( start<=edep->Time() && edep->Time()<=stop) {
                    float energy = edep->Energy();
                    totalEallTracks += energy;
                    int tid = edep->TrackID();
                    // Chase the tid up to the parent which started in the gas.
                    simb::MCParticle* TPCEve = FindTPCEve(tid);
                    int tidTPC;
                    if (TPCEve!=nullptr) {
                        tidTPC = TPCEve->TrackId();
                    } else {
                        tidTPC = tid;
                    }
                    auto tidmapiter = tidmap.find(tidTPC);
                    if (tidmapiter == tidmap.end()) {
                        calIDEs.emplace_back(tidTPC, 0.0, energy);
                        tidmap[tidTPC] = calIDEs.size()-1;
                    } else {
                        calIDEs.at(tidmapiter->second).energyTot += energy;
                    }
                }
            }

            // loop over the hitIDEs to set the fractional energy for each TrackID
            for (size_t i = 0; i < calIDEs.size(); ++i) {
                calIDEs[i].energyFrac = calIDEs[i].energyTot / totalEallTracks;
            }
            return calIDEs;
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::CaloHitPurity(simb::MCParticle* const p,
            std::vector<art::Ptr<gar::rec::CaloHit> > const& hits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::CaloHitPurity")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }


            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            float total   = 0.0;
            for (auto hit : hits) {
                std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                weight = (weightByCharge) ? hit->Energy() : 1.0;
                total += weight;

                for (auto iCalIDE : calIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT trackIDs
                    if (iCalIDE.trackID    == trackIDin &&
                        iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            double purity      = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                purity      = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(purity,uncertainty);
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::CaloHitEfficiency(simb::MCParticle* const p,
            std::vector<art::Ptr<gar::rec::CaloHit> > const& hits,
            std::vector<art::Ptr<gar::rec::CaloHit> > const& allhits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::CaloHitEfficiency")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            for (auto hit : hits) {
                std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                weight = (weightByCharge) ? hit->Energy() : 1.0;

                for (auto iCalIDE : calIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT track IDs
                    if (iCalIDE.trackID    == trackIDin &&
                        iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            // Next figure out how many hits in the whole collection are associated with this id
            float total   = 0.0;
            if (desired>0.0) {
                for (auto hit : allhits) {
                    std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                    weight = (weightByCharge) ? hit->Energy() : 1.0;

                    for (auto iCalIDE : calIDEs) {
                        // Don't double count if this hit has > 1 of the desired GEANT track IDs
                        if (iCalIDE.trackID    == trackIDin &&
                            iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                            total += weight;
                             break;
                        }
                    }
                 } // end loop over all hits
            }

            double efficiency  = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                efficiency  = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(efficiency,uncertainty);
        }




    
    } // cheat
} // gar
