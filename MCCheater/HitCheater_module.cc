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
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// GArSoft Includes
#include "MCCheater/BackTracker.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "ReconstructionDataProducts/Hit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

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
      
      void CreateHitsOnChannel(unsigned int                                  channel,
                               std::vector<rec::Hit>                       & chanHits,
                               std::vector<const sdp::EnergyDeposit*> const& edeps);
      
      std::string                         fReadoutLabel; ///< label of module creating raw digits
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock
      
    };
    
  } // cheat
  
  namespace cheat {

    //--------------------------------------------------------------------------
    HitCheater::HitCheater(fhicl::ParameterSet const& pset)
    {
      fTime  = gar::providerFrom<detinfo::DetectorClocksService>();

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
      fReadoutLabel = pset.get<std::string>("ReadoutModuleLabel", "rawdata");
      
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
      
      auto digCol = evt.getValidHandle< std::vector<raw::RawDigit> >(fReadoutLabel);
      
      // get the FindMany for the digit to EnergyDeposit association
      ::art::FindMany<sdp::EnergyDeposit>    fmed(hitCol, evt, fReadoutLabel);
      std::vector<const sdp::EnergyDeposit*> edeps;
      
      if(!digCol.isValid() ||
         !fmed  .isValid() ){
        throw cet::exception("BackTracker")
        << "Unable to find valid collection of RawDigits "
        << digCol.isValid()
        << " or FindMany<EnergyDeposit> "
        << fmed.isValid()
        << " this is a problem for cheating";
      }

      // make a container to store the hits for each channel
      std::vector< std::vector<rec::Hit> > hitsByChannel;
      hitsByChannel.resize(digCol->size());
      
      // loop over the raw digits, ignore any that have no associated energy deposits
      // those digits are just noise
      unsigned int channel = 0;
      std::vector<rec::Hit> hits;
      for(size_t rd = 0; rd < digCol->size(); ++rd){
        fmed.get(d, edeps);
        if(edeps.size() < 1) continue;

        channel = (*digCol)[rd].Channel();
        this->CreateHitsOnChannel(channel, hitsByChannel[channel], edeps);
        
        // add the hits to the output collection and make the necessary associations
        for(auto hit : hits){
          hitCol->push_back(hit);
          
          auto const ptr = art::Ptr<raw::RawDigit>(digCol, rd);
          util::CreateAssn(*this, evt, *hitCol, ptr, *rdHitAssns);

        }
        
      } // end loop over raw digits
      
      evt.put(std::move(hitCol    ));
      evt.put(std::move(rdHitAssns));
      evt.put(std::move(pHitAssns ));
      
      return;
    }
    
    //--------------------------------------------------------------------------
    bool sortEDepByTime(const sdp::EnergyDeposit* a, const sdp::EnergyDeposit* b)
    {
      return a->Time() < b->Time();
    }
    
    //--------------------------------------------------------------------------
    void HitCheater::CreateHitsOnChannel(unsigned int                                  channel,
                                         std::vector<rec::Hit>                       & chanHits,
                                         std::vector<const sdp::EnergyDeposit*> const& edeps)
    {
      // make sure there are no left overs in the hit vector
      chanHits.clear();
      
      // loop over the energy deposits and create a list for each track id
      std::map<int, std::vector<const sdp::EnergyDeposit*> > trackIDToEDep;
      
      for(auto const* edep : edeps) trackIDToEDep[std::abs(edep.TrackID())] = edep;

      // now for each track ID identify hits as being those energy depositions
      // that are within a contiguous block of TDC values
      unsigned short prevTDC = 0;
      unsigned short tdc     = 0;
      
      std::vector<const sdp::EnergyDeposit*> subSet;
      
      for(auto const itr : trackIDToEDep){

        auto edepCol = itr.second;
        
        std::sort(edepCol.begin(), edepCol.end(), sortEDepByTime);
        prevTDC = fTime->TPCG4Time2TDC(edepCol.first()->Time());
        subSet.clear();
        
        for(auto const* edep : edepCol){
          
          tdc = fTime->TPCG4Time2TDC(edep->Time());
          
          if(tdc > prevTDC + 1){
            // make a new hit - use the weighted average position and summed
            // energy
            chanHits.emplace
            
            subSet.clear();
          }

          subSet.push_back(edep);
        } // end loop over the energy deposits for this track ID
        
      } // end loop to make hits
      
      
      return;
    }

    //--------------------------------------------------------------------------

    
  } // cheat
  
  namespace cheat {
      
    DEFINE_ART_MODULE(IonizationReadout)
    
  } // cheat
} // gar
#endif // GAR_MCCHEATER_HITCHEATER
