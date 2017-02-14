//
//  TPCReadoutSimStandardAlg.h
//
//  Created by Brian Rebel on 2/14/17.
//

#ifndef GAR_READOUTSIM_TPCReadoutSimStandardAlg_h
#define GAR_READOUTSIM_TPCReadoutSimStandardAlg_h

#include "ReadoutSimulation/TPCReadoutSimAlg.h"

namespace gar{
  namespace rosim{
    
    class TPCReadoutSimStanardAlg : public TPCReadoutSimAlg {
      
      
      TPCReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
                               fhicl::ParameterSet    const& pset);
      virtual ~TPCReadoutSimStandardAlg();
      
      std::vector<raw::RawDigit> CreateRawDigits(std::vector<gar::sdp::IDE> const& ides);
      
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      // Methods to handle the noise simulation - CreateNoiseDigits is for
      // channels without any signal on them, AddNoiseToSignalDigits is for
      // channels where signal is also recorded.  Assume that the noise level
      // is dependent on the amount of signal recorded
      std::vector<raw::RawDigit> CreateNoiseDigits();
      void                       AddNoiseToADCs(raw::ADCvector_t & adcs);

      
    };
    
  }
}

#endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
