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
    {
    }
    
    //----------------------------------------------------------------------------
    // place the pixels such that pixel 0 is at the bottom of the upstream end
    // and the pixels increase along the z direction.  The second row starts
    // with pixel number 0 + numPixelsZ and so on
    void ChannelMapStandardAlg::Initialize(GeometryCore &geo)
    {
      // start over:
      Uninitialize();
      
      // get the extent of the detector
      fHalfHeight = geo.DetHalfHeight();
      auto length = geo.DetLength();
      
      // count up how many channels
      // for this channel map, we have square pixels of a size given by the
      // configuration as is the pitch between them
      fNumPixelsZ = (unsigned int)std::floor(length / fPixelPitch);
      fNumPixelsY = (unsigned int)std::floor(2. * fHalfHeight / fPixelPitch);
      
      // get the xyz center for each pixel
      for(size_t y = 0; y < fNumPixelsY; ++y){
        fPixRowBounds.insert((y + 1) * fPixelPitch - fHalfHeight);
        for(size_t z = 0; z < fNumPixelsZ; ++z){
          if(y == 0) fPixColumnBounds.insert((z + 1) * fPixelPitch);
          fPixelCenters.emplace_back(std::numeric_limits<float>::max(),
                                     (y + 0.5) * fPixelPitch - geo.DetHalfHeight(),
                                     (z + 0.5) * fPixelPitch);
        }
      }
      
      LOG_DEBUG("ChannelMapStandardAlg")
      << "There are "
      << fPixRowBounds.size()
      << " pixel rows and "
      << fPixColumnBounds.size()
      << " pixel columns";
      
//      int ctr = 0;
//      for(auto y : fPixRowBounds){
//        LOG_VERBATIM("ChannelMapStandardAlg")
//        << "pixel "
//        << ctr
//        << " y bound: "
//        << y;
//
//        ++ctr;
//      }
//      
//      ctr = 0;
//      for(auto z : fPixColumnBounds){
//        LOG_VERBATIM("ChannelMapStandardAlg")
//        << "pixel "
//        << ctr
//        << " z bound: "
//        << z;
//        
//        ++ctr;
//      }

      return;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Uninitialize()
    {
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::Nchannels() const
    {
      return (fNumPixelsZ * fNumPixelsY);
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::NearestChannel(float const* xyz) const
    {
      // use the bounds to determine the channel
      auto yItr = fPixRowBounds.lower_bound(xyz[1]);
      auto zItr = fPixColumnBounds.lower_bound(xyz[2]);
      
      if(yItr == fPixRowBounds.end()    ||
         zItr == fPixColumnBounds.end() ){
        
        throw cet::exception("NearestChannel")
        << "y position: "
        << xyz[1]
        << " or z position: "
        << xyz[2]
        << " is out of bounds.  Max y: "
        << *(fPixRowBounds.rbegin())
        << " Max z: "
        << *(fPixColumnBounds.rbegin());
        
      }
      
      int row = std::floor((*yItr + fHalfHeight) / fPixelPitch) - 1;
      int col = std::floor((*zItr) / fPixelPitch) - 1;
      
      // we can't have negative pixel rows
      if(row < 0) row = 0;
      
      unsigned int chan = (unsigned int)((row * fNumPixelsZ) + col);
      
      return chan;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::ChannelToPosition(unsigned int chan,
                                                  float*       xyz)  const
    {
      xyz[0] = fPixelCenters[chan].x;
      xyz[1] = fPixelCenters[chan].y;
      xyz[2] = fPixelCenters[chan].z;
      
      return;
    }
    
  } // namespace
} // namespace gar