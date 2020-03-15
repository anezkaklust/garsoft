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

// Timing
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
            fRawDataLabel         = pset.get<std::string>("RawDataLabel",            "daq");
            fMinHitEnergyFraction = pset.get<double     >("MinimumHitEnergyFraction", 0.1);
            fDisableRebuild       = pset.get<bool       >("DisableRebuild",           false);

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

            // Timing
            //time_point start(std::chrono::high_resolution_clock::now());
            //auto totalThere = start -start;

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

                // Timing
                //start       = high_resolution_clock::now();           
                hids = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());
                // Timing
                //totalThere += high_resolution_clock::now() -start;     

                for (auto const& hid : hids) {
                    if ( hid.trackID==tkID && hid.energyFrac>fMinHitEnergyFraction ) {
                        hitList.push_back(hit);
                    }
                }
            }
            // Timing
            //duration<double, std::micro> dA = totalThere;
            //std::cout << "Time in ChannelToHitIDES: " <<
            //    dA.count() << " for " << allhits.size() << " hits" << std::endl;
            return hitList;
        }



        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::ChannelToHitIDEs(raw::Channel_t const& channel,
                                          double const start, double const stop) const {

            // loop over the electrons in the channel and grab those that are in time
            // with the identified hit start and stop times

            std::vector<HitIDE> hitIDEs;

            double totalEallTracks = 0.0;

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
        BackTrackerCore::HitCollectionPurity(simb::MCParticle* const p,
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
        BackTrackerCore::HitCollectionEfficiency(simb::MCParticle* const p,
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

            double efficiency  = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                efficiency  = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(efficiency,uncertainty);
        }






    
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    // plist is assumed to have adopted the appropriate EveIdCalculator prior to
    // having been passed to this method. It is likely that the EmEveIdCalculator is
    // the one you always want to use
    std::vector<HitIDE> BackTrackerCore::HitToEveID(::art::Ptr<gar::rec::Hit> const& hit) const
    {
      std::vector<HitIDE> ides = this->HitToHitIDEs(hit);
      
      // make a map of evd ID values and fraction of energy represented by
      // that eve id in this hit
      std::map<int, HitIDE> eveToE;
      
      for(auto const& hid : ides){
        eveToE[fParticleList.EveId( hid.trackID )] = hid;
      }
      
      // now fill the eveides vector from the map
      std::vector<HitIDE> eveides;
      eveides.reserve(eveToE.size());
      for(auto itr : eveToE){
        eveides.emplace_back(itr.first,
                             itr.second.energyFrac,
                             itr.second.energyTot);
      }
      
      return eveides;
    }
    
    
    //----------------------------------------------------------------------
    std::set<int> BackTrackerCore::GetSetOfTrackIDs() const
    {
      // fParticleList::value_type is a pair (track, particle pointer)
      std::set<int> trackIDs;
      for (auto const& pl : fParticleList)
        trackIDs.insert(pl.first);
      
      return trackIDs;
    }    
        
    //--------------------------------------------------------------------------
    std::map<size_t, const sdp::EnergyDeposit*> BackTrackerCore::ChannelTDCToEnergyDeposit(raw::Channel_t channel) const
    {
      if(channel > fChannelToEDepCol.size())
        throw cet::exception("BackTrackerCore")
        << "Requested channel "
        << channel
        << " for backtracking TDCs is out of range, bail";
      
      LOG_DEBUG("BackTracker")
      << "There are "
      << fChannelToEDepCol[channel].size()
      << " energy depositions for channel "
      << channel;
      
      float chanpos[3]={0.,0.,0};
      fGeo->ChannelToPosition(channel, chanpos);

      std::map<size_t, const sdp::EnergyDeposit*> edeps;
      
      // loop over all the energy deposits for this channel
      size_t centralTDC = 0.;

      for(auto const* edep : fChannelToEDepCol[channel]){
        centralTDC = (size_t)fClocks->TPCG4Time2TDC(edep->Time() + TMath::Abs(edep->X()-chanpos[0]) * fInverseVelocity);
        
        LOG_DEBUG("BackTrackerCore")
        << "centralTDC: "
        << centralTDC
        << " time "
        << edep->Time()
        << " energy "
        << edep->Energy();
        
        edeps[centralTDC] = edep;
      }
      
      return edeps;
    }
    
    
    //----------------------------------------------------------------------
    std::vector<float> BackTrackerCore::HitToXYZ(::art::Ptr<gar::rec::Hit> const& hit) const
    {
      std::vector<float> xyz;
      
      // There is a data member which is a vector of all energy deposits for each channel
      // in the event.  Let's look to see if the vector is non-empty for this one
      // If it is empty, it could be a noise hit
      if(fChannelToEDepCol.size() < hit->Channel()){
        LOG_WARNING("BackTracker")
        << "Attempting to back track a hit from a channel without any corresponding EnergyDeposits, "
        << "return empty vector: "
        << fChannelToEDepCol.size()
        << "/"
        << hit->Channel();
        
        return xyz;
      }
      
      if(fChannelToEDepCol[hit->Channel()].size() < 1){
        LOG_WARNING("BackTracker")
        << "Attempting to back track a hit without any corresponding EnergyDeposits, "
        << "return empty vector";
        return xyz;
      }
      
      // We know time extent of the hit, so that should allow us to narrow down
      // the EnergyDeposits
      // Since multiple times correspond to a single hit, we need to
      // average over all contributing positions
      float avX    = 0.;
      float avY    = 0.;
      float avZ    = 0.;
      float wgtSum = 0.;
      
      unsigned int start = (unsigned int)hit->StartTime();
      unsigned int tdc   = 0;
      unsigned int end   = (unsigned int)hit->EndTime();
      
      for(auto const& edep : fChannelToEDepCol[hit->Channel()]){
        
        // check the TDC value of the deposit to make sure it is in the
        // correct range
        // TODO: strictly speaking, the TDC retrieved this way may not be correct
        // because we are getting a single TDC value for the EnergyDeposit, but
        // the simulation breaks the deposits into clusters of electrons to
        // account for diffusion.  This method is approximately correct
        tdc = fClocks->TPCG4Time2TDC( edep->Time() );
        if(tdc < start || tdc > end) continue;
        
        avX += edep->X() * edep->Energy();
        avY += edep->Y() * edep->Energy();
        avZ += edep->Z() * edep->Energy();
        
        wgtSum += edep->Energy();
      }
      
      if(wgtSum == 0.){
        LOG_WARNING("BackTracker")
        << "no energy recorded for any EnergyDeposit belonging to this hit,"
        << " return an empty vector";
        return xyz;
      }
      
      xyz.push_back(avX / wgtSum);
      xyz.push_back(avY / wgtSum);
      xyz.push_back(avZ / wgtSum);
      
      return xyz;
    }

  } // cheat
} // gar
