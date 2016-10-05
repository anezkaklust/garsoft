/** ****************************************************************************
 * @file RawDigit.cxx
 * @brief Definition of basic raw digits
 * @author brebel@fnal.gov
 * @see  RawDigit.h raw.h
 * 
 * Compression/uncompression utilities are declared in lardata/RawData/raw.h .
 * 
 * ****************************************************************************/

#include "RawDataProducts/RawDigit.h"

// C/C++ standard libraries

namespace gar {
  namespace raw{
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit()
    : fADC(0)
    , fChannel(InvalidChannelID)
    , fSamples(0)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
      //----------------------------------------------------------------------
    RawDigit::RawDigit(ChannelID_t                  channel,
                       unsigned short               samples,
                       RawDigit::ADCvector_t const& adclist)
    : fADC(adclist)
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit(ChannelID_t             channel,
                       unsigned short          samples,
                       RawDigit::ADCvector_t&& adclist)
    : fADC(std::move(adclist))
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
    //----------------------------------------------------------------------
    void RawDigit::SetPedestal(float ped, float sigma /* = 1. */ )
    {
      
      fPedestal = ped;
      fSigma = sigma;
      
    } // RawDigit::SetPedestal()
    
    
  } // namespace raw
} // gar
////////////////////////////////////////////////////////////////////////

