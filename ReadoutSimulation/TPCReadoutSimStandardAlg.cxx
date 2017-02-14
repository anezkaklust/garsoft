//
//  TPCReadoutSimStandardAlg.cxx
//
//  Created by Brian Rebel on 2/14/17.
//

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"


namespace gar {
  namespace rosim{
    
    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::TPCReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
                                                       fhicl::ParameterSet    const& pset)
    :TPCReadoutSimAlg(engine, pset)
    {
      this->reconfigure(pset);
      
      return;
    }
    
    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::~TPCReadoutSimStandardAlg()
    {
      return;
    }

    //----------------------------------------------------------------------------
    std::vector<raw::RawDigit> CreateRawDigits(std::vector<gar::sdp::IDE> const& ides)
    {
      std::vector<raw::RawDigit> digits;
      
      // loop over the IDEs to figure out where the channel boundaries are
      // then take them one channel at a time to convert the number of electrons
      // drifted as a function of TDC value to ADC values
      
      if(fAddNoise){
        // now add noise only digits to the collection
        auto this->CreateNoiseDigits();
      } // end if we are simulating noise

      return  digits;
    }
    
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
      fAddNoise = pset.get<bool>("AddNoise", false);
      
      return;
    }
    
    //----------------------------------------------------------------------------
    std::vector<raw::RawDigit> TPCReadoutSimStandardAlg::CreateNoiseDigits()
    {
      std::vector<raw::RawDigit> noiseDigits;
      
      // figure out which channels we need to make pure noise
      
      return noiseDigits;
    }

    
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::AddNoiseToADCs(raw::ADCvector_t & adcs)
    {
      for(auto adc : adcs){
        
        // TODO: figure out how to add noise
        adc += 0;
        
      }
      
      return;
    }

    
  } // rosim
} // gar