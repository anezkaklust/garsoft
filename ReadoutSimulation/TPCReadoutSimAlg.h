//
//  TPCReadoutSimAlg.h
//
//  Created by Brian Rebel on 2/14/17.
//

#ifndef GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp

#include "CLHEP/Random/RandGauss.h"

#include "RawDataProducts/RawDigit.h"

namespace fhicl {
  class ParameterSet;
}

namespace gar {

  namespace sdp {
    class IDE;
  }
  
  namespace rosim {
    
    class TPCReadoutSimAlg{
      
    public:
      
      TPCReadoutSimAlg(CLHEP::HepRandomEngine      & engine,
                       fhicl::ParameterSet    const& pset);
      
      virtual ~TPCReadoutSimAlg();
      
      // Method to take IDEs and turn them into RawDigits
      virtual std::vector<raw::RawDigit> CreateRawDigits(std::vector<gar::sdp::IDE> const& ides) = 0;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
      
    protected:
      
      CLHEP::HepRandomEngine & fEngine;   ///< random number engine
      bool                     fAddNoise; ///< flag to add noise or not
      
    private:
      
      // Methods to handle the noise simulation - CreateNoiseDigits is for
      // channels without any signal on them, AddNoiseToSignalDigits is for
      // channels where signal is also recorded.  Assume that the noise level
      // is dependent on the amount of signal recorded
      virtual std::vector<raw::RawDigit> CreateNoiseDigits()                     = 0;
      virtual void                       AddNoiseToADCs(raw::ADCvector_t & adcs) = 0;

    };
    
  } // end rosim
  
} // end gar


#endif /* GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp */
