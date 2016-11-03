////////////////////////////////////////////////////////////////////////
// $Id: BackTracker_service.cc,v 1.3 2011/12/13 05:57:02 bckhouse Exp $
//
//
// \file: BackTracker_service.cc
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include <map>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "SimulationDataProducts/sim.h"

namespace gar{
  namespace cheat{
    
    //----------------------------------------------------------------------
    BackTracker::BackTracker(fhicl::ParameterSet   const& pset,
                             art::ActivityRegistry      & reg)
    {
      reconfigure(pset);
      
      reg.sPreProcessEvent.watch(this, &BackTracker::Rebuild);
    }
    
    //----------------------------------------------------------------------
    BackTracker::~BackTracker()
    {
    }
    
    //----------------------------------------------------------------------
    void BackTracker::reconfigure(const fhicl::ParameterSet& pset)
    {
      fG4ModuleLabel        = pset.get<std::string>("G4ModuleLabel",           "gargeant");
      fMinHitEnergyFraction = pset.get<double     >("MinimumHitEnergyFraction", 0.1);
    }
    
    //----------------------------------------------------------------------
    void BackTracker::Rebuild(art::Event const& evt)
    {
      // do nothing if this is data
      if(evt.isRealData()) return;
      
      // get the particles from the event
      auto pHandle = evt.getValidHandle<std::vector<simb::MCParticle> >(fG4ModuleLabel);
      
      // first check to see if we got called to early, that is the particles are
      // and sim channels are not made yet in a MC production job
      // if that is the case, we'll take care of it later
      if(pHandle.failedToGet()){
        LOG_WARNING("BackTracker")
        << "failed to get handle to simb::MCParticle from "
        << fG4ModuleLabel
        << ", return";
        return;
      }
      
      // Clear out anything remaining from previous calls to Rebuild
      fParticleList.clear();
      fMCTruthList .clear();
      fSimChannels .clear();
      
      art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
      
      if( fo.isValid() ){
        for(size_t p = 0; p < pHandle->size(); ++p){
          
          simb::MCParticle *part = new simb::MCParticle(pHandle->at(p));
          fParticleList.Add(part);
          
          // get the simb::MCTruth associated to this sim::ParticleList
          try{
            art::Ptr<simb::MCTruth> mct = fo.at(p);
            if(fMCTruthList.size() < 1) fMCTruthList.push_back(mct);
            else{
                // check that we are not adding a simb::MCTruth twice to the collection
                // we know that all the particles for a given simb::MCTruth are put into the
                // collection of particles at the same time, so we can just check that the
                // current art::Ptr has a different id than the last one put
              if(!(mct == fMCTruthList.back())) fMCTruthList.push_back(mct);
            }
            // fill the track id to mctruth index map
            fTrackIDToMCTruthIndex[pHandle->at(p).TrackId()] = fMCTruthList.size() - 1;
          }
          catch(cet::exception &ex){
            LOG_WARNING("BackTracker")
            << "unable to find MCTruth from ParticleList "
            << "created in "
            << fG4ModuleLabel
            << " any attempt to get the MCTruth objects from "
            << "the backtracker will fail\n"
            << "message from caught exception:\n"
            << ex;
          }
        }// end loop over particles to get MCTruthList
      }// end if fo.isValid()
      
      // grab the sim::SimChannels for this event
      evt.getView(fG4ModuleLabel, fSimChannels);
      
      // grab the voxel list for this event
      //fVoxelList = sim::SimListUtils::GetLArVoxelList(evt, fG4ModuleLabel);
      
      fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
      
      LOG_DEBUG("BackTracker")
      << "BackTracker has "
      << fSimChannels.size()
      << " sim::SimChannels and "
      << GetSetOfTrackIDs().size()
      << " tracks.  The particles are:\n"
      << fParticleList
      << "\n the MCTruth information is\n";
      
      for(auto mc : fMCTruthList)
        LOG_DEBUG("BackTracker") << *(mc);
      
      return;
    }
    
    //----------------------------------------------------------------------
    const simb::MCParticle* BackTracker::TrackIDToParticle(int const& id) const
    {
      auto part_it = fParticleList.find(id);
      
      if(part_it == fParticleList.end()){
        LOG_WARNING("BackTracker")
        << "can't find particle with track id "
        << id
        << " in sim::ParticleList returning null pointer";
        
        return 0;
      }
      
      return part_it->second;
    }
    
    //----------------------------------------------------------------------
    const simb::MCParticle* BackTracker::TrackIDToMotherParticle(int const& id) const
    {
      // get the mother id from the particle navigator
      // the EveId was adopted in the Rebuild method
      
      return this->TrackIDToParticle(fParticleList.EveId(abs(id)));
    }
    
    //----------------------------------------------------------------------
    const art::Ptr<simb::MCTruth>& BackTracker::TrackIDToMCTruth(int const& id) const
    {
        // find the entry in the MCTruth collection for this track id
      size_t mct = fTrackIDToMCTruthIndex.find(std::abs(id))->second;
      
      if(/* mct < 0 || */ mct > fMCTruthList.size() )
        throw cet::exception("BackTracker")
        << "attempting to find MCTruth index for "
        << "out of range value: "
        << mct
        << "/" << fMCTruthList.size();
      
      return fMCTruthList[mct];
    }
    
    //----------------------------------------------------------------------
    std::vector<sdp::IDE> BackTracker::TrackIDToSimIDE(int const& id) const
    {
      std::vector<sdp::IDE> ides;
      
      // loop over all sim::SimChannels and fill a vector
      // of sim::IDE objects for the given track id
      for(auto sc : fSimChannels){
        const auto & tdcidemap = sc->TDCIDEs();
        
        // loop over the IDEMAP
        for(auto const& mapitr : tdcidemap){
          
          // loop over the vector of IDE objects.
          auto const& idevec = mapitr.fIDEs;
          for(auto const& ide : idevec){
            if( std::abs(ide.trackID) == id) ides.push_back(ide);
          }
          
        } // end loop over map from sim::SimChannel
      } // end loop over sim::SimChannels
      
      return ides;
    }
    
    //----------------------------------------------------------------------
    const art::Ptr<simb::MCTruth>& BackTracker::ParticleToMCTruth(const simb::MCParticle* p) const
    {
      return this->TrackIDToMCTruth(p->TrackId());
    }
    
    //----------------------------------------------------------------------
    std::vector<const simb::MCParticle*> BackTracker::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
    {
      std::vector<const simb::MCParticle*> ret;
      
        // sim::ParticleList::value_type is a pair (track ID, particle pointer)
      for(auto const& TrackIDpair: fParticleList) {
        if( TrackIDToMCTruth(TrackIDpair.first) == mct )
          ret.push_back(TrackIDpair.second);
      }
      
      return ret;
    }
    
    //----------------------------------------------------------------------
    std::vector<sdp::TrackIDE> BackTracker::HitToTrackID(art::Ptr<gar::rec::Hit> const& hit)
    {
      std::vector<sdp::TrackIDE> trackIDEs;
      
      const double start = hit->StartTime();
      const double end   = hit->EndTime();
      
      this->ChannelToTrackID(trackIDEs, hit->Channel(), start, end);
      
      return trackIDEs;
    }
    
    //----------------------------------------------------------------------
    const std::vector<std::vector<art::Ptr<gar::rec::Hit>>> BackTracker::TrackIDsToHits(std::vector<art::Ptr<gar::rec::Hit>> const& allhits,
                                                                                     std::vector<int>                  const& tkIDs)
    {
      // returns a subset of the hits in the allhits collection that are matched
      // to MC particles listed in tkIDs
      
      // temporary vector of TrackIDs and Ptrs to hits so only one
      // loop through the (possibly large) allhits collection is needed
      std::vector<std::pair<int, art::Ptr<gar::rec::Hit>>> hitList;
      std::vector<sdp::TrackIDE> tids;
      for(auto hit : allhits){
        tids.clear();

        this->ChannelToTrackID(tids,
                               hit->Channel(),
                               hit->StartTime(),
                               hit->EndTime());
        
        for(auto const& tid : tids) {
          for(auto const& itkid : tkIDs) {
            if(tid.trackID == itkid) {
              if(tid.energyFrac > fMinHitEnergyFraction)
                hitList.push_back(std::make_pair(itkid, hit));
            }
          } // itkid
        } // itid
      } // itr
      
      // now build the truHits vector that will be returned to the caller
      std::vector<std::vector<art::Ptr<gar::rec::Hit>>> truHits;
      
      // temporary vector containing hits assigned to one MC particle
      std::vector<art::Ptr<gar::rec::Hit>> tmpHits;
      for(auto const& itkid : tkIDs) {
        tmpHits.clear();
        
        for(auto itr : hitList) {
          if(itkid == itr.first) tmpHits.push_back(itr.second);
        }
        truHits.push_back(tmpHits);
      }
      
      return truHits;
    }
    
    //----------------------------------------------------------------------
    // plist is assumed to have adopted the appropriate EveIdCalculator prior to
    // having been passed to this method. It is likely that the EmEveIdCalculator is
    // the one you always want to use
    std::vector<sdp::TrackIDE> BackTracker::HitToEveID(art::Ptr<gar::rec::Hit> const& hit)
    {
      std::vector<sdp::TrackIDE> trackides = this->HitToTrackID(hit);
      
      // make a map of evd ID values and fraction of energy represented by
      // that eve id in this hit
      std::map<int, float> eveToE;
      
      double totalE = 0.;
      for(auto const& tid : trackides){
        eveToE[fParticleList.EveId( tid.trackID )] += tid.energy;
        totalE += tid.energy;
      }
      
      // now fill the eveides vector from the map
      std::vector<sdp::TrackIDE> eveides;
      eveides.reserve(eveToE.size());
      for(auto itr : eveToE){
        sdp::TrackIDE temp;
        temp.trackID    = itr.first;
        temp.energyFrac = itr.second/totalE;
        temp.energy     = itr.second;
        eveides.push_back(std::move(temp));
      }
      
      return eveides;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTracker::GetSetOfEveIDs()
    {
      std::set<int> eveIDs;
      
      sim::ParticleList::const_iterator plitr = fParticleList.begin();
      while(plitr != fParticleList.end() ){
        int eveID = fParticleList.EveId((*plitr).first);
          // look to see if this eveID is already in the set
        if( eveIDs.find(eveID) == eveIDs.end() ) eveIDs.insert(eveID);
        plitr++;
      }
      
      return eveIDs;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTracker::GetSetOfTrackIDs()
    {
      // fParticleList::value_type is a pair (track, particle pointer)
      std::set<int> trackIDs;
      for (auto const& pl : fParticleList)
        trackIDs.insert(pl.first);
      
      return trackIDs;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTracker::GetSetOfEveIDs(std::vector< art::Ptr<gar::rec::Hit> > const& hits)
    {
      std::set<int> eveIDs;
      
      std::vector< art::Ptr<gar::rec::Hit> >::const_iterator itr = hits.begin();
      for(auto itr : hits){
        
        // get the eve ids corresponding to this hit
        const std::vector<sdp::TrackIDE> ides = HitToEveID(itr);
        
        // loop over the ides and extract the track ids
        for(auto const& ide : ides) eveIDs.insert(ide.trackID);
        
      }
      
      return eveIDs;
    }
    
      //----------------------------------------------------------------------
    std::set<int> BackTracker::GetSetOfTrackIDs(std::vector< art::Ptr<gar::rec::Hit> > const& hits)
    {
      std::set<int> trackIDs;
      std::vector<sdp::TrackIDE> trackIDEs;
      
      for(auto itr : hits ){
        
        trackIDEs.clear();
        
        // get the track ids corresponding to this hit
        const double start = itr->StartTime();
        const double end   = itr->EndTime();
        
        this->ChannelToTrackID(trackIDEs,
                               itr->Channel(),
                               start,
                               end);
        
        // loop over the ides and extract the track ids
        for(auto const& tid : trackIDEs) {
          trackIDs.insert(tid.trackID);
        }
      }
      
      return trackIDs;
    }
    
    //----------------------------------------------------------------------
    double BackTracker::HitCollectionPurity(std::set<int>                              trackIDs,
                                            std::vector< art::Ptr<gar::rec::Hit> > const& hits,
                                            bool                                       weightByCharge)
    {
      // get the list of EveIDs that correspond to the hits in this collection
      // if the EveID shows up in the input list of trackIDs, then it counts
      float total   = 0;
      float desired = 0.;
      float weight  = 1.;
      
      // don't have to check the view in the hits collection because
      // those are assumed to be from the object we are testing and will
      // the correct view by definition then.
      for(auto hit : hits){
        
        weight = (weightByCharge) ? hit->Signal() : 1.;
        
        std::vector<sdp::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        total += weight;
        
        // don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(auto const& tid : hitTrackIDs){
          if(trackIDs.find(tid.trackID) != trackIDs.end()){
            desired += weight;
            break;
          }
        }
        
      }// end loop over hits
      
      double purity = 0.;
      
      if(total > 0) purity = desired/total;
      
      return purity;
    }
    
    
    //----------------------------------------------------------------------
    double BackTracker::HitCollectionEfficiency(std::set<int>                              trackIDs,
                                                std::vector< art::Ptr<gar::rec::Hit> > const& hits,
                                                std::vector< art::Ptr<gar::rec::Hit> > const& allhits,
                                                bool                                          weightByCharge)
    {
      // get the list of EveIDs that correspond to the hits in this collection
      // and the energy associated with the desired trackID
      float desired = 0.;
      float total   = 0.;
      float weight  = 0.;
      
      // don't have to check the view in the hits collection because
      // those are assumed to be from the object we are testing and will
      // the correct view by definition then.
      for(auto hit : hits){
        
        std::vector<sdp::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);

        weight = (weightByCharge) ? hit->Signal() : 1.;
        
        // don't worry about hits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(auto tid : hitTrackIDs){
          if(trackIDs.find(tid.trackID) != trackIDs.end() &&
             tid.energyFrac             >= fMinHitEnergyFraction){
            desired += weight;
            break;
          }
        }
      }// end loop over hits
      
      // now figure out how many hits in the whole collection are associated with this id
      for(auto hit : allhits){
        
        weight = (weightByCharge) ? hit->Signal() : 1.;

        std::vector<sdp::TrackIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        for(auto const& tid : hitTrackIDs){
          // don't worry about hits where the energy fraction for the chosen
          // trackID is < 0.1
          // also don't double count if this hit has more than one of the
          // desired track IDs associated with it
          if(trackIDs.find(tid.trackID) != trackIDs.end() &&
             tid.energyFrac             >= fMinHitEnergyFraction){
            total += weight;
            break;
          }
        }
        
      }// end loop over all hits
      
      double efficiency = 0.;
      if(total > 0.) efficiency = desired/total;
      
      return efficiency;
    }
    
    //----------------------------------------------------------------------
    const sdp::SimChannel* BackTracker::FindSimChannel(raw::Channel_t const& channel) const
    {
      const sdp::SimChannel* chan = nullptr;
      
      for(auto sc : fSimChannels){
        if(sc->Channel() == channel) chan = sc;
      }
      
      if(!chan)
        throw cet::exception("BackTracker")
        << "No sim::SimChannel corresponding "
        << "to channel: "
        << channel;
      
      return chan;
    }
    
    //----------------------------------------------------------------------
    void BackTracker::ChannelToTrackID(std::vector<sdp::TrackIDE>      & trackIDEs,
                                       raw::Channel_t             const& channel,
                                       double                            hit_start_time,
                                       double                            hit_end_time)
    {
      trackIDEs.clear();
      
      double totalE = 0.;
      
      try{
        const sdp::SimChannel* schannel = this->FindSimChannel(channel);
        
        // loop over the electrons in the channel and grab those that are in time
        // with the identified hit start and stop times
        const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
        int start_tdc = ts->TPCTick2TDC( hit_start_time );
        int end_tdc   = ts->TPCTick2TDC( hit_end_time   );
        if(start_tdc<0) start_tdc = 0;
        if(end_tdc<0) end_tdc = 0;
        std::vector<sdp::IDE> simides = schannel->TrackIDsAndEnergies(start_tdc, end_tdc);
        
        // first get the total energy represented by all track ids for
        // this channel and range of tdc values
        for(auto const& ide : simides) totalE += ide.energy;
        
        // protect against a divide by zero below
        if(totalE < 1.e-5) totalE = 1.;
        
        // loop over the entries in the map and fill the input vectors
        
        for(auto const& ide : simides){
          
          if(ide.trackID == gar::sdp::NoParticleId) continue;
          
          sdp::TrackIDE info;
          info.trackID    = ide.trackID;
          info.energyFrac = ide.energy/totalE;
          info.energy     = ide.energy;
          
          trackIDEs.push_back(info);
          
        }
      }// end try
      catch(cet::exception e){
        LOG_WARNING("BackTracker")
        << "caught exception \n"
        << e;
      }
      
      return;
    }
    
    //----------------------------------------------------------------------
    void BackTracker::HitToSimIDEs(gar::rec::Hit           const& hit,
                                   std::vector<sdp::IDE>     & ides) const
    {
      // Get services.
      const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
      
      int start_tdc = ts->TPCTick2TDC( hit.StartTime() );
      int end_tdc   = ts->TPCTick2TDC( hit.EndTime()   );
      if(start_tdc<0) start_tdc = 0;
      if(end_tdc<0) end_tdc = 0;
      
      ides = FindSimChannel(hit.Channel())->TrackIDsAndEnergies(start_tdc, end_tdc);
    }
    
    //----------------------------------------------------------------------
    std::vector<double> BackTracker::SimIDEsToXYZ(std::vector<sdp::IDE> const& ides)
    {
      std::vector<double> xyz(3, std::numeric_limits<double>::min());
      
      double x = 0.;
      double y = 0.;
      double z = 0.;
      double w = 0.;
      
      // loop over electrons.
      
      for(auto const& ide : ides) {
        
        double weight = ide.numElectrons;
        
        w += weight;
        x += weight * ide.x;
        y += weight * ide.y;
        z += weight * ide.z;
        
      }// end loop over sim::IDEs
      
      // if the sum of the weights is still 0, then return
      // the obviously stupid default values
      if(w < 1.e-5)
        throw cet::exception("BackTracker")
        << "No sim::IDEs providing non-zero number of electrons"
        << " can't determine originating location from truth\n";
      
      xyz[0] = x/w;
      xyz[1] = y/w;
      xyz[2] = z/w;
      
      return xyz;
    }
    
    //----------------------------------------------------------------------
    std::vector<double> BackTracker::HitToXYZ(art::Ptr<gar::rec::Hit> const& hit)
    {
      std::vector<sdp::IDE> ides;
      HitToSimIDEs(hit, ides);
      return SimIDEsToXYZ(ides);
    }
    
    DEFINE_ART_SERVICE(BackTracker)

  } // namespace cheat
} // namespace gar