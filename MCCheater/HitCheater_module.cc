////////////////////////////////////////////////////////////////////////
/// \file  HitCheater_module.cc
/// \brief Create rec::Hit products using MC truth information
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_MCCHEATER_HITCHEATER
#define GAR_MCCHEATER_HITCHEATER

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "cetlib_except/exception.h"

#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft Includes
#include "MCCheater/BackTracker.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"
#include "ReconstructionDataProducts/Hit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// ROOT Includes

#include "TMath.h"

// Forward declarations

namespace gar {
  namespace cheat {
    
    class HitCheater : public ::art::EDProducer{
    public:
      
      // Standard constructor and destructor for an FMWK module.
      explicit HitCheater(fhicl::ParameterSet const& pset);
      virtual ~HitCheater();
      
      void produce (::art::Event& evt);
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      void CreateHitsOnChannel(gar::raw::RawDigit                         const& rawDigit,
                               std::vector<std::pair<int, gar::rec::Hit> >     & chanHits,
                               std::vector<const gar::sdp::EnergyDeposit*>     & edeps);

      std::string                         fReadoutLabel; ///< label of module creating raw digits
      std::string                         fG4Label;      ///< label of module creating mc particles
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock
      unsigned int                        fTDCGap;       ///< gap allowed between energy depositions
                                                         ///< from the same G4 track for making hits
      const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
      const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
      
    };
    
  } // cheat
  
  namespace cheat {

    //--------------------------------------------------------------------------
    HitCheater::HitCheater(fhicl::ParameterSet const& pset)
    {
      fTime    = gar::providerFrom<detinfo::DetectorClocksService>();
      fGeo     = gar::providerFrom<geo::Geometry>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      
      this->reconfigure(pset);
      
      produces< std::vector<rec::Hit>                    >();
      produces< ::art::Assns<raw::RawDigit,    rec::Hit> >();
      produces< ::art::Assns<simb::MCParticle, rec::Hit> >();
      
      return;
    }

    //--------------------------------------------------------------------------
    HitCheater::~HitCheater()
    {
      return;
    }

    //--------------------------------------------------------------------------
    void HitCheater::reconfigure(fhicl::ParameterSet const& pset)
    {
      fReadoutLabel = pset.get<std::string >("ReadoutModuleLabel", "daq"  );
      fG4Label      = pset.get<std::string >("G4ModuleLabel",      "geant");
      fTDCGap       = pset.get<unsigned int>("TDCGapAllowed"              );
      return;
    }
    
    //--------------------------------------------------------------------------
    void HitCheater::produce(::art::Event & evt)
    {
      // Hits know about the xyz reconstructed position and the start and end
      // times of the hits, as well as the channel number and signal deposited
      
      // Grab the RawDigits from the event and their associations to EnergyDeposits
      // Then sort the EnergyDeposits for each RawDigit by TrackID and time
      // We will create a hit for each clump of EnergyDeposits in time for a give
      // channel number
      
      // We will want to associate the hits to their raw digits, but also to
      // the MCParticle they correspond to
      
      if( evt.isRealData() ){
        throw cet::exception("HitCheater")
        << "Attempting to cheat the hit reconstruction for real data, "
        << "that will never work";
        return;
      }
      
      std::unique_ptr< std::vector<rec::Hit>                    > hitCol    (new std::vector<rec::Hit>                    );
      std::unique_ptr< ::art::Assns<raw::RawDigit,    rec::Hit> > rdHitAssns(new ::art::Assns<raw::RawDigit,    rec::Hit> );
      std::unique_ptr< ::art::Assns<simb::MCParticle, rec::Hit> > pHitAssns (new ::art::Assns<simb::MCParticle, rec::Hit> );
      
      // get the raw digits from the event
      auto digCol = evt.getValidHandle< std::vector<gar::raw::RawDigit> >(fReadoutLabel);
      
      // get the FindMany for the digit to EnergyDeposit association
      ::art::FindMany<sdp::EnergyDeposit> fmed(digCol, evt, fReadoutLabel);
      
      // test that the FindMany is valid - don't worry about the digCol, the
      // getValidHandle throws if it can't find a valid handle
      if(!fmed.isValid() ){
        throw cet::exception("HitCheater")
        << "Unable to find valid FindMany<EnergyDeposit> "
        << fmed.isValid()
        << " this is a problem for cheating";
      }

      // get the MCParticles from the event and make a vector of Ptrs of them
      auto partCol = evt.getValidHandle< std::vector<simb::MCParticle> >(fG4Label);
      std::vector<art::Ptr<simb::MCParticle> > partVec;
      art::fill_ptr_vector(partVec, partCol);
      std::map<int, size_t> idToPart;
      for(size_t p = 0; p < partVec.size(); ++p) idToPart[partVec[p]->TrackId()] = p;
      
      // loop over the raw digits, ignore any that have no associated energy
      // deposits those digits are just noise
      std::vector<std::pair<int, rec::Hit> > hits;
      std::vector<const sdp::EnergyDeposit*> edeps;
      
      for(size_t rd = 0; rd < digCol->size(); ++rd){
        edeps.clear();
        fmed.get(rd, edeps);
        if(edeps.size() < 1) continue;
        
        auto const& dig = (*digCol)[rd];
        
        LOG_DEBUG("HitCheater")
        << "There are "
        << edeps.size()
        << " energy deposits for channel "
        << dig.Channel();
        
        this->CreateHitsOnChannel(dig, hits, edeps);
        
        if(hits.size() > 0)
          LOG_DEBUG("HitCheater")
          << "There are "
          << hits.size()
          << " hits found for channel "
          << dig.Channel();
        
        // add the hits to the output collection and make the necessary associations
        for(auto hit : hits){
          
          LOG_DEBUG("HitCheater")
          << hit.second;
          
          hitCol->push_back(hit.second);

          // make the digit to hit association
          auto const digPtr = art::Ptr<raw::RawDigit>(digCol, rd);
          util::CreateAssn(*this, evt, *hitCol, digPtr, *rdHitAssns);
          
          // make the hit to MCParticle association
          util::CreateAssn(*this, evt, *hitCol, partVec[idToPart[hit.first]], *pHitAssns);

        }
        
      } // end loop over raw digits
      
      evt.put(std::move(hitCol    ));
      evt.put(std::move(rdHitAssns));
      evt.put(std::move(pHitAssns ));
      
      return;
    }
    
    //--------------------------------------------------------------------------
    void HitCheater::CreateHitsOnChannel(gar::raw::RawDigit                         const& rawDigit,
                                         std::vector<std::pair<int, gar::rec::Hit> >     & chanHits,
                                         std::vector<const gar::sdp::EnergyDeposit*>     & edeps)
    {
      // make sure there are no left overs in the hit vector
      chanHits.clear();
      
      // use the back tracker to get the energy deposits for each TDC value in
      // this channel.  Create a map of TDC and TrackID to energy deposits.
      // Then loop over the map to make hits out of deposits from the same
      // TrackID in contiguous TDC values
      
      auto bt = gar::providerFrom<cheat::BackTracker>();
      
      std::map<int, std::map<size_t, float> > trackIDToTDCEDeps;
      std::vector<float>                      totalTDCEDep(rawDigit.Samples(), 0.);
      int                                     nonZeroCtr = 0;

      auto tdcEDeps = bt->ChannelTDCToEnergyDeposit(rawDigit.Channel());
      LOG_DEBUG("HitCheater")
      << "There are "
      << tdcEDeps.size()
      << " tdcEDeps for channel "
      << rawDigit.Channel();

      // uncompress Raw Digits as needed, since samples are in ticks

      gar::raw::ADCvector_t uncompressed;
      uncompressed.resize(rawDigit.Samples());
      short ped = rawDigit.Pedestal();
      gar::raw::Uncompress(rawDigit.ADCs(),uncompressed,ped,rawDigit.Compression());

      for(size_t tdc = 0; tdc < rawDigit.Samples(); ++tdc){

        if(uncompressed[tdc] != 0){
          ++nonZeroCtr;
          LOG_DEBUG("HitCheater")
          << "RawDigit on channel "
          << rawDigit.Channel()
          << " has nonzero signal "
          << uncompressed[tdc]
          << " at "
          << tdc;
        }

        auto itr = tdcEDeps.find(tdc);
        if(itr == tdcEDeps.end() ) continue;
        
        totalTDCEDep[tdc] += itr->second->Energy();
        trackIDToTDCEDeps[std::abs(itr->second->TrackID())][tdc] += itr->second->Energy();
        
      } // end loop to sort the energy deposits
      
      LOG_DEBUG("HitCheater")
      << "There are "
      << trackIDToTDCEDeps.size()
      << " entries in the track ID to TDCEDeps map";
      
      // no signal on the channel, so return empty handed
      if(nonZeroCtr == 0){
        LOG_DEBUG("HitCheater")
        << "The input raw digit on channel "
        << rawDigit.Channel()
        << " has no signal for any TDC tick but "
        << trackIDToTDCEDeps.size()
        << " energy depositions";
      
        return;
      }
      
      // now for each trackID, find the associated hits
      // Make a hit for each TDC with signal for the track ID to start with
      // and then later combine those which are close enough together
      unsigned int channel = rawDigit.Channel();
      
      // get the YZ position for this hit, will calculate x based on the
      // time below
      float  pos[3]   = {0.};
      fGeo->ChannelToPosition(channel, pos);
      
      LOG_DEBUG("HitCheater")
      << "Channel location is ("
      << pos[0]
      << ", "
      << pos[1]
      << ", "
      << pos[2]
      << ")";
      
      // search for gaps in the TDC ticks with signal for each track id
      // that means the outter loop has to be over track id values
      // then loop over the tdcs for each track id and look for gaps to
      // create hits.  It is possible there is only one hit on a given
      // channel for a given track, but there could be more than one too
      // if the volume is magnetized and the particle momentum is low enough
      // to allow it to curl around.
      int          trkId  = 0;
      unsigned int curTDC = 0;
      unsigned int begTDC = 0;
      unsigned int endTDC = 0;
      float        hitSig = 0.;
      float        hitTime = 0.;
      float        hitSumSq = 0;

      for(auto trkIdItr : trackIDToTDCEDeps){
        
        trkId = trkIdItr.first;
        
        // loop over the tdcs for this track id
        auto const& tdcToEDepMap = trkIdItr.second;

        LOG_DEBUG("HitCheater")
        << "TrackID "
        << trkId
        << " has the following TDCs: ";
//        for(auto const& itr : tdcToEDepMap)
//          LOG_VERBATIM("HitCheater")
//          << itr.first;
        
        std::vector<gar::rec::Hit> trkIDHits;
        std::set<unsigned int>     usedTDCs;
        
        for(auto itr : tdcToEDepMap){
          
          // the x position is given by the product of the drift velocity
          // and the time, given by the weighted TDC average converted to
          // a time
          pos[0] = fDetProp->DriftVelocity() * fTime->TPCTick2Time(curTDC);

          curTDC  = itr.first;
          
          begTDC = (curTDC > fTDCGap            ) ? curTDC - fTDCGap : 0;
          endTDC = (curTDC < rawDigit.Samples() ) ? curTDC + fTDCGap : rawDigit.Samples();
          
          LOG_DEBUG("HitCheater")
          << "begTDC: "
          << begTDC
          << " endTDC: "
          << endTDC;

          // check with ADC values from the raw digit are close enough in
          // time to assign to the current true energy deposit and make hits out
          // of those.  Make sure to not double count

          hitSig = 0.;
	  hitTime = 0;
	  hitSumSq = 0;

          for(size_t tdc = begTDC; tdc < endTDC + 1; ++tdc){
            // only add signal from tdcs which have not already been used for
            // this track id
            if(usedTDCs.count(tdc) < 1){
	      float cursig =  1. * uncompressed[tdc] * itr.second / totalTDCEDep[curTDC];
	      hitSig += cursig;
              usedTDCs.insert(tdc);
	      hitTime += hitSig*tdc;
              hitSumSq = hitSig*tdc*tdc;

              LOG_DEBUG("HitCheater")
              << " tdc: "
              << tdc
              << " hitSig: "
              << hitSig
              << " "
              << uncompressed[tdc]
              << " "
              << itr.second / totalTDCEDep[curTDC];
            }
          }

	  float hitRMS = 0;
	  if (hitSig>0)
	    {
	      hitTime = hitTime/hitSig;
	      hitRMS = TMath::Sqrt( TMath::Max( (float) 0.0,hitSumSq/hitSig - hitTime*hitTime) );
	    }
	  else
	    {
	      hitTime = 0.5*(begTDC+endTDC);
	      hitRMS = 0;
	    }

          // the second arugment is the signal from the current TrackID
          // in this TDC
          trkIDHits.emplace_back(channel,
                                 hitSig,
                                 pos,
                                 begTDC,
                                 endTDC,
				 hitTime,
				 hitRMS);
          
        } // end loop over TDC to energy deposit map
        
        // now add the hits for this track ID to the output vector
        for(auto const& hit : trkIDHits)
          chanHits.emplace_back(std::make_pair(trkId, hit));
        
      } // end loop over track ids
      
      return;
    }
    
  } // cheat
  
  namespace cheat {
    
    DEFINE_ART_MODULE(HitCheater)
    
  } // cheat
} // gar
#endif // GAR_MCCHEATER_HITCHEATER
