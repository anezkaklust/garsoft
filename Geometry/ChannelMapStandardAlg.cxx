////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapStandardAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace gar {
  namespace geo{
    
    //----------------------------------------------------------------------------
    ChannelMapStandardAlg::ChannelMapStandardAlg(fhicl::ParameterSet const& p)
    : fNumPixelsZ     (std::numeric_limits<unsigned int>::max())
    , fNumPixelsY     (std::numeric_limits<unsigned int>::max())
    , fPixelActiveSize(p.get<float>("PixelActiveSize",  0.3 )  )
    , fPixelPitch     (p.get<float>("PixelPitch",       0.35)  )
    , fPixelSize      (fPixelActiveSize + fPixelPitch          )
    {
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Initialize(GeometryCore &geo)
    {
      // start over:
      Uninitialize();
      
      // get the extent of the detector
      auto height   = geo.DetHalfHeight() * 2.;
      auto length   = geo.DetLength();
      
      // count up how many channels
      // for this channel map, we have square pixels of a size given by the
      // configuration as is the pitch between them
      fNumPixelsZ = (unsigned int)std::floor(length / fPixelSize);
      fNumPixelsY = (unsigned int)std::floor(height / fPixelSize);
      
      return;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Uninitialize()
    {
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::Nchannels() const
    {
      return (fNumPixelsZ + fNumPixelsY);
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::NearestChannel(float const* xyz) const
    {
      // place the pixels such that pixel 0 is at the bottom of the upstream end
      // and the pixels increase along the z direction.  The second row starts
      // with pixel number 0 + numPixelsZ and so on
      
      float maxZ = fNumPixelsZ * fPixelSize;
      float maxY = 0.5 * (fNumPixelsY * fPixelSize);
      
      // make sure the given position has a channel even associated with it
      if(std::abs(xyz[1]) > maxY ||
         xyz[2]           < 0.   ||
         xyz[2]           > maxZ)
        throw cet::exception("ChannelMapStandard")
        << "Closest channel requested for position outside of active volume: ("
        << xyz[0]
        << ", "
        << xyz[1]
        << ", "
        << xyz[2]
        << ")";
      
      // now figure out which pixel we are dealing with
      unsigned int zPix = (unsigned int)std::floor(xyz[2] / fPixelSize);
      
      // for the y pixel number, we have to account for the fact that
      // the y position goes between -maxY < 0 < +maxY
      unsigned int yPix = (unsigned int)std::floor((xyz[1] + maxY) / fPixelSize);
      
      return (yPix * fNumPixelsZ + zPix);      
    }
    
  } // namespace
} // namespace gar