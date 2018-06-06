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
    , fCompression(gar::raw::kNone)
    {}
    
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit(Channel_t                    channel,
                       unsigned short               samples,
                       RawDigit::ADCvector_t const& adclist,
		       gar::raw::Compress_t compress)
    : fADC(adclist)
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    , fCompression(compress)
    {}
    
    
    //----------------------------------------------------------------------
    RawDigit::RawDigit(Channel_t               channel,
                       unsigned short          samples,
                       RawDigit::ADCvector_t&& adclist,
		       gar::raw::Compress_t compress)
    : fADC(std::move(adclist))
    , fChannel(channel)
    , fSamples(samples)
    , fPedestal(0.)
    , fSigma(0.)
    , fCompression(compress)
    {}
    
    
  } // namespace raw
} // gar
////////////////////////////////////////////////////////////////////////

