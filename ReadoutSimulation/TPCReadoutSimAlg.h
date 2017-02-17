//
//  TPCReadoutSimAlg.h
//
//  Created by Brian Rebel on 2/14/17.
//

#ifndef GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp

#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGauss.h"

#include "RawDataProducts/RawDigit.h"
#include "DetectorInfo/DetectorProperties.h"

namespace fhicl {
  class ParameterSet;
}

namespace gar {

  namespace rosim {
    
    // helper struct for tracking the information needed to make raw digits
    struct edepIDE {
      
      edepIDE();
      edepIDE(float          numE,
              unsigned int   chan,
              unsigned short tdc,
              size_t         edepIdx)
      : NumElect(numE)
      , Channel (chan)
      , TDC     (tdc)
      {
        edepLocs.push_back(edepIdx);
      }
      
      //------------------------------------------------------------------------
      void AddEDep(size_t edepLoc) { edepLocs.push_back(edepLoc); }
      
      //------------------------------------------------------------------------
      bool operator <(edepIDE const& b) const
      {
        if(Channel < b.Channel) return true;
        else if(Channel == b.Channel){
          if(TDC < b.TDC) return true;
        }
        
        return false;
      }
      
      //------------------------------------------------------------------------
      void operator +=(edepIDE const& b)
      {
        if(TDC != b.TDC ){
          LOG_WARNING("IonizationReadout")
          << "Attempting to add edepIDE with different "
          << "TDCs: "
          << TDC
          << " / "
          << b.TDC
          << " bail";
          return;
        }
        
        NumElect += b.NumElect;
        
        for(auto const e : b.edepLocs)
          edepLocs.push_back(e);
        
        return;
      }
      
      float               NumElect;
      unsigned int        Channel;
      unsigned short      TDC;
      std::vector<size_t> edepLocs;
    };

    class TPCReadoutSimAlg{
      
    public:
      
      TPCReadoutSimAlg(CLHEP::HepRandomEngine      & engine,
                       fhicl::ParameterSet    const& pset);
      
      virtual ~TPCReadoutSimAlg();
      
      // Method to take IDEs and turn them into RawDigits
      virtual raw::RawDigit CreateRawDigit(unsigned int              channel,
                                           std::vector<float> const& electrons) = 0;
      virtual void  CreateNoiseDigits(std::vector<raw::RawDigit> & digits)      = 0;


      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
      
    protected:
      
      CLHEP::HepRandomEngine &           fEngine;   ///< random number engine
      bool                               fAddNoise; ///< flag to add noise or not
      const detinfo::DetectorProperties* fDetProp;  ///< detector properties
      
      // AddNoiseToADCs is foradding noise to recorded signal.  Assume that the
      // noise level is dependent on the amount of signal recorded
      virtual void  AddNoiseToADCs(std::vector<short> & adcs) = 0;
      virtual short ElectronsToADCs(float electrons)          = 0;

    };
    
  } // end rosim
  
} // end gar


#endif /* GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp */
