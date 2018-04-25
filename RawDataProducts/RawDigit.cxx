/** ****************************************************************************
 * @file RawDigit.cxx
 * @brief Definition of basic raw digits
 * @author brebel@fnal.gov
 * @see  RawDigit.h raw.h
 * 
 * 
 * ****************************************************************************/

#include <limits>

#include "RawDataProducts/RawDigit.h"

// C/C++ standard libraries

namespace gar {
  namespace raw{
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit()
    : fADC(0)
    , fChannel(std::numeric_limits<unsigned int>::max())
    , fSamples(0)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit(Channel_t                    channel,
                       unsigned short               samples,
                       RawDigit::ADCvector_t const& adclist)
    : fADC(adclist)
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit(Channel_t               channel,
                       unsigned short          samples,
                       RawDigit::ADCvector_t&& adclist)
    : fADC(std::move(adclist))
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    {}
    
    
  } // namespace raw
} // gar
////////////////////////////////////////////////////////////////////////

