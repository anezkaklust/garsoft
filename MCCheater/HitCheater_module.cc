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
#include "DetectorInfo/DetectorPropertiesService.h"
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
      
      void CreateHitsOnChannel(gar::raw::RawDigit                         const& rawDigit,
                               std::vector<std::pair<int, gar::rec::Hit> >     & chanHits,
                               std::vector<const gar::sdp::EnergyDeposit*>     & edeps);
      
      std::string                         fReadoutLabel; ///< label of module creating raw digits
      std::string                         fG4Label;      ///< label of module creating mc particles
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock
      unsigned int                        fTDCWindow;    ///< gap allowed between energy depositions
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
      fTDCWindow    = pset.get<unsigned int>("TDCGapAllowed"              );
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
      ::art::FindMany<sdp::EnergyDeposit>    fmed(digCol, evt, fReadoutLabel);
      
      // test that the FindMany is valid - don't worry about the digCol, the
      // getValidHandle throws if it can't find a valid handle
      if(!fmed.isValid() ){
        throw cet::exception("BackTracker")
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
      
      // loop over the raw digits, ignore any that have no associated energy deposits
      // those digits are just noise
      std::vector<std::pair<int, rec::Hit> > hits;
      std::vector<const sdp::EnergyDeposit*> edeps;
      
      for(size_t rd = 0; rd < digCol->size(); ++rd){
        fmed.get(rd, edeps);
        if(edeps.size() < 1) continue;

        auto const& dig = (*digCol)[rd];
        this->CreateHitsOnChannel(dig, hits, edeps);
        
        // add the hits to the output collection and make the necessary associations
        for(auto hit : hits){
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
      
      std::map<int, std::vector<float> > trackIDToTDCEDeps;
      std::vector<float>                 totalTDCEDep(rawDigit.Samples(), 0.);
      
      for(size_t tdc = 0; tdc < rawDigit.Samples(); ++tdc){
        auto tdcEDeps = bt->ChannelTDCToEnergyDeposit(rawDigit.Channel(), tdc);
        
        for(auto const* itr : tdcEDeps){
          totalTDCEDep[tdc] += itr->Energy();
          trackIDToTDCEDeps[std::abs(itr->TrackID())][tdc] += itr->Energy();
        }

      } // end loop to sort the energy deposits
      
      // now for each trackID, find the associated hits
      size_t startTDC = 0;
      float  hitSig   = 0.;
      float  totEDep  = 0.;
      float  trkEDep  = 0.;
      float  wgtTDC   = 0.;
      
      // get the YZ position for this hit, will calculate x based on the
      // time below
      float  pos[3]   = {0.};
      fGeo->ChannelToPosition(rawDigit.Channel(), pos);
      
      for(auto itr : trackIDToTDCEDeps){
       
        auto const& tdcVec = itr.second;
        // Loop over the vector of tdcs, make hits out of contiguous tdcs with
        // energy deposits.
        // Use the raw signal for the tdcs to sum up the total hit signal, then
        // multiply that number by the ratio of the energy deposited for the track
        // to the energy deposited for all tracks in the tdc range to get the
        // signal for just this track/hit
        startTDC = 0;
        wgtTDC   = 0.;
        trkEDep  = 0.;
        totEDep  = 0.;
        hitSig   = 0.;
        for(size_t tdc = 0; tdc < tdcVec.size(); ++tdc){

          if(tdcVec[tdc] == 0.){
            for(size_t t = startTDC; t < tdc; ++t){
              wgtTDC  += tdcVec[t] * t;
              hitSig  += 1. * rawDigit.ADC(t);
              trkEDep += tdcVec[t];
              totEDep += totalTDCEDep[t];
            } // end loop over TDCs with signal
            
            // make the hit
            if(totEDep == 0.){
              LOG_WARNING("HitCheater")
              << "total energy deposited for tdc range "
              << startTDC
              << " : "
              << tdc
              << " is zero, that shouldn't be, do nothing with this range";
              continue;
            }
            
            hitSig *= trkEDep / totEDep;
            wgtTDC /= totEDep;
            
            // the x position is given by the product of the drift velocity
            // and the time, given by the weighted TDC average converted to
            // a time
            pos[0] = fDetProp->DriftVelocity() * fTime->TPCTick2Time(wgtTDC);
            
            chanHits.push_back(std::make_pair(itr.first,
                                              gar::rec::Hit(rawDigit.Channel(),
                                                            hitSig,
                                                            pos,
                                                            startTDC,
                                                            tdc - 1)
                                              )
                               );
            
            startTDC = tdc;
            hitSig   = 0.;
            trkEDep  = 0.;
            totEDep  = 0.;
            wgtTDC   = 0.;
          } // end if we have a gap in TDCs with signal
          
        } // end loop over TDCs
          
      }
      
      return;
    }

    //--------------------------------------------------------------------------

    
  } // cheat
  
  namespace cheat {
      
    DEFINE_ART_MODULE(HitCheater)
    
  } // cheat
} // gar
#endif // GAR_MCCHEATER_HITCHEATER
