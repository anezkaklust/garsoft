//
//  TPCReadoutSimStandardAlg.cxx
//
//  Created by Brian Rebel on 2/14/17.
//

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar {
  namespace rosim{
    
    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::TPCReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
                                                       fhicl::ParameterSet    const& pset)
    : TPCReadoutSimAlg(engine, pset)
    {
      this->reconfigure(pset);
      
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      
      return;
    }
    
    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::~TPCReadoutSimStandardAlg()
    {
      return;
    }

    //----------------------------------------------------------------------------
    raw::RawDigit TPCReadoutSimStandardAlg::CreateRawDigit(unsigned int              channel,
                                                           std::vector<float> const& electrons)
    {
      // make a vector to store the adc information
      auto numTicks = fDetProp->NumberTimeSamples();
      std::vector<short> adcs(numTicks, 0);
      
      
      // check that the size of the electrons vector is the same as the
      // adc vector, they better be
      if(adcs.size() != electrons.size()){
        throw cet::exception("TPCReadoutSimStandardAlg")
        << "The vectors for adc and electrons have different sizes, no idea "
        << "how to handle that, bail";
      }
      
      for(size_t e = 0; e < electrons.size(); ++e){
        adcs[e] = this->ElectronsToADCs(electrons[e]);
      }
      
      // add noise to the adc vector if we want that
      if(fAddNoise) this->AddNoiseToADCs(adcs);
      
      return  raw::RawDigit(channel, adcs.size(), adcs);
    }
    
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
      fAddNoise = pset.get<bool>("AddNoise", false);
      
      return;
    }
    
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::CreateNoiseDigits(std::vector<raw::RawDigit> & digits)
    {
      if(!fAddNoise) return;
      
      // we'll need the geometry
      auto geo = gar::providerFrom<geo::Geometry>();
      
      // figure out which channels we need to make pure noise
      std::set<unsigned int> channels;
      for(auto const& dig : digits) channels.insert(dig.Channel());

      // now loop over the channels in the detector and make a noise digit
      // for any channel that does not already have signal on it
      for(unsigned int c = 0; c < geo->NChannels(); ++c){
      
        // do nothing if we already have a digit for this channel
        if(channels.count(c) > 0) continue;
        
        std::vector<short> adcs(fDetProp->NumberTimeSamples(), 0);
        this->AddNoiseToADCs(adcs);
        
        digits.emplace_back(c, adcs.size(), adcs);

      } // end make noise digits
      
      
      return;
    }
    
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::AddNoiseToADCs(std::vector<short> & adcs)
    {
      for(auto adc : adcs){
        
        // TODO: figure out how to add noise
        // It is assumed that the amount of noise added to
        // a given ADC value depends on the ADC value
        adc += 0;
        
      }
      
      return;
    }

    //----------------------------------------------------------------------------
    short TPCReadoutSimStandardAlg::ElectronsToADCs(float electrons)
    {
      return (short)(fDetProp->ElectronsToADC() * electrons);
    }
    
  } // rosim
} // gar