/// $Id: SimChannel.cxx,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  SimulationDataProducts/SimChannel.cxx
///
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#include <limits> // std::numeric_limits
#include <utility>
#include <stdexcept>

#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/sim.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace sim{
    
    //--------------------------------------------------------------------------
    bool TDCIDE::operator<(TDCIDE const& other) const
    {
      return this->fTDC < other.fTDC;
    }
    
    //-------------------------------------------------
    IDE::IDE()
    : trackID     (std::numeric_limits<int>::min())
    , numElectrons(std::numeric_limits<float>::max())
    , energy      (std::numeric_limits<float>::max())
    , x           (std::numeric_limits<float>::max())
    , y		        (std::numeric_limits<float>::max())
    , z		        (std::numeric_limits<float>::max())
    {}
    
    //-------------------------------------------------
    IDE::IDE(sim::IDE const& ide,
             int             offset)
    : trackID     (ide.trackID+offset)
    , numElectrons(ide.numElectrons)
    , energy      (ide.energy)
    , x           (ide.x)
    , y           (ide.y)
    , z           (ide.z)
    {}
    
    // Default constructor
    //-------------------------------------------------
    SimChannel::SimChannel()
    : fChannel(std::numeric_limits<unsigned int>::max())
    {}
    
    //-------------------------------------------------
    SimChannel::SimChannel(unsigned int channel)
    : fChannel(channel)
    {}
    
    //-------------------------------------------------
    void SimChannel::AddIonizationElectrons(int          trackID,
                                            unsigned int tdc,
                                            double       numberElectrons,
                                            double      *xyz,
                                            double       energy)
    {
      // look at the map to see if the current TDC already
      // exists, if not, add it, if so, just add a new track id to the
      // vector
      
      // no electrons? no energy? no good!
      if((numberElectrons < std::numeric_limits<double>::epsilon()) ||
         (energy         <= std::numeric_limits<double>::epsilon())){
        // will throw
        LOG_ERROR("SimChannel")
        << "AddIonizationElectrons() trying to add to TDC #"
        << tdc
        << " "
        << numberElectrons
        << " electrons with "
        << energy
        << " MeV of energy from track ID = "
        << trackID;
        return;
      } // if no energy or no electrons
      
      gar::sim::TDCIDE temp;
      temp.fTDC = tdc;
      
      auto const& itr = fTDCIDEs.find(temp);
      if(itr != fTDCIDEs.end()){
        
        // loop over the IDE vector for this tdc and add the electrons
        // to the entry with the same track id
        for(auto & ide : itr->fIDEs){
          
          if( ide.trackID == trackID ){
            // make a weighted average for the location information
            double weight    = ide.numElectrons + numberElectrons;
            temp.fIDEs.emplace_back(trackID,
                                    weight,
                                    ide.energy + energy,
                                    (ide.x * ide.numElectrons + xyz[0] * numberElectrons) / weight,
                                    (ide.y * ide.numElectrons + xyz[1] * numberElectrons) / weight,
                                    (ide.z * ide.numElectrons + xyz[2] * numberElectrons) / weight);
          }
          else temp.fIDEs.push_back(ide);
        } // loop over the IDEs for this tdc

        // now remove the TDCIDE from the set and replace it with the updated one
        fTDCIDEs.erase(itr);
        fTDCIDEs.insert(std::move(temp));
        
      } // if this tdc value is in the set
      else{
        // if we never found the tdc, add a TDCIDE to the vector
        sim::TDCIDE tdcide;
        tdcide.fTDC = tdc;
        tdcide.fIDEs.emplace_back(trackID,
                                  numberElectrons,
                                  energy,
                                  xyz[0],
                                  xyz[1],
                                  xyz[2]);
        
        fTDCIDEs.insert(tdcide);
      }
      
      return;
    }
    
    //-------------------------------------------------
    gar::sim::TDCIDE const& SimChannel::FindTDCIDE(unsigned short tdc) const
    {
      gar::sim::TDCIDE test;
      test.fTDC = tdc;
      
      auto const& found = fTDCIDEs.find(test);
      if(found != fTDCIDEs.end()){
        if(found->fTDC == tdc){
          return *found;
        }
      }
      
      throw cet::exception("SimChannel")
      << "no TDCIDE found for "
      << tdc;
      
      return *(fTDCIDEs.begin());
    }
    
    //-------------------------------------------------
    double SimChannel::Charge(unsigned int tdc) const
    {
      double charge = 0.;

      try{
        auto const& tdcide = this->FindTDCIDE(tdc);

        // loop over the list for this tdc value and add up
        // the total number of electrons
        for(auto const& ide : tdcide.fIDEs){
          charge += ide.numElectrons;
        }
      }
      catch(cet::exception & e){}
      
      return charge;
    }
    
    //-------------------------------------------------
    double SimChannel::Energy(unsigned int tdc) const
    {
      double energy = 0.;

      try{
        auto const& tdcide = this->FindTDCIDE(tdc);
        
        // loop over the list for this tdc value and add up
        // the total energy
        for(auto const& ide : tdcide.fIDEs){
          energy += ide.energy;
        }
      }
      catch(cet::exception & e){}
      
      return energy;
    }
    
    //-----------------------------------------------------------------------
    // the start and end tdc values are assumed to be inclusive
    std::vector<sim::IDE> SimChannel::TrackIDsAndEnergies(unsigned int startTDC,
                                                          unsigned int endTDC) const
    {
      // make a map of track ID values to sim::IDE objects
      std::map<int, sim::IDE> idToIDE;
      
      std::vector<sim::IDE> ides;
      
      if(startTDC > endTDC ){
        mf::LogWarning("SimChannel") << "requested tdc range is bogus: "
        << startTDC << " " << endTDC
        << " return empty vector";
        return ides;
      }
      
      gar::sim::TDCIDE start;
      gar::sim::TDCIDE end;
      start.fTDC = startTDC;
      end  .fTDC = endTDC;
      
      auto mitr  = fTDCIDEs.lower_bound(start);
      auto stItr = fTDCIDEs.lower_bound(start);
      auto edItr = fTDCIDEs.upper_bound(end);
      
      for(mitr = stItr; mitr != edItr; mitr++){
        
        // grab the vector of IDEs for this tdc
        auto const& idelist = (*mitr).fIDEs;

        // now loop over them and add their content to the map
        for(auto const& ide : idelist){
          
          if( idToIDE.find(ide.trackID) != idToIDE.end() ){
            double nel1   = idToIDE[ide.trackID].numElectrons;
            double nel2   = ide.numElectrons;
            double en1    = idToIDE[ide.trackID].energy;
            double en2	  = ide.energy;
            double energy = en1+en2;
            double weight = nel1 + nel2;
            // make a weighted average for the location information
            idToIDE[ide.trackID].x            = (ide.x * nel2 + idToIDE[ide.trackID].x * nel1) / weight;
            idToIDE[ide.trackID].y            = (ide.y * nel2 + idToIDE[ide.trackID].y * nel1) / weight;
            idToIDE[ide.trackID].z            = (ide.z * nel2 + idToIDE[ide.trackID].z * nel1) / weight;
            idToIDE[ide.trackID].numElectrons = weight;
            idToIDE[ide.trackID].energy       = energy;
          } // end if the track id for this one is found
          else{
            sim::IDE temp(ide);
            idToIDE[ide.trackID] = temp;
          }
        } // end loop over vector
      } // end loop over tdc values
      
      // now fill the vector with the ides from the map
      for(auto const& itr : idToIDE){
        ides.push_back(itr.second);
      }
      
      return ides;
    }
    
    //-----------------------------------------------------------------------
    // the start and end tdc values are assumed to be inclusive
    std::vector<sim::TrackIDE>  SimChannel::TrackIDEs(unsigned int startTDC,
                                                      unsigned int endTDC) const
    {
      
      std::vector<sim::TrackIDE> trackIDEs;
      
      if(startTDC > endTDC ){
        LOG_WARNING("SimChannel::TrackIDEs")
        << "requested tdc range is bogus: "
        << startTDC << " " << endTDC
        << " return empty vector";
        return trackIDEs;
      }
      
      double totalE = 0.;
      const std::vector<sim::IDE> ides = TrackIDsAndEnergies(startTDC,endTDC);
      for (auto const& ide : ides)
        totalE += ide.energy;
      
      // protect against a divide by zero below
      if(totalE < 1.e-5) totalE = 1.;
      
      // loop over the entries in the map and fill the input vectors
      for (auto const& ide : ides){
        if(ide.trackID == sim::NoParticleId) continue;
        trackIDEs.emplace_back(ide.trackID,ide.energy/totalE,ide.energy);
      }
      
      
      return trackIDEs;
    }
    
    //-----------------------------------------------------------------------
    // Merge the collection of IDEs from one sim channel to another.
    // Requires an agreed upon offset for G4 trackID
    std::pair<int,int> SimChannel::MergeSimChannel(SimChannel const& channel,
                                                   int               offset)
    {
      if( this->Channel() != channel.Channel() )
        throw std::runtime_error("ERROR SimChannel Merge: Trying to merge different channels!");
      
      std::pair<int,int> range_trackID(std::numeric_limits<int>::max(),
                                       std::numeric_limits<int>::min());
      
      gar::sim::TDCIDE tdcide;
      
      for(auto const& tdc : channel.TDCIDEs()){
        tdcide.fIDEs.clear();
        tdcide.fTDC = tdc.fTDC;

        for(auto const& ide : tdc.fIDEs){
          
          tdcide.fIDEs.emplace_back(ide, offset);

          if(ide.trackID+offset  < range_trackID.first)
            range_trackID.first  = ide.trackID+offset;
          if(ide.trackID+offset  > range_trackID.second)
            range_trackID.second = ide.trackID+offset;
          
        }//end loop over IDEs
        
        
      }//end loop over TDCIDEMap
      
      return range_trackID;
      
    }
    
  } // sim
}// gar
