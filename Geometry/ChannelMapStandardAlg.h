////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELSTANDARDMAPALG_H
#define GEO_CHANNELSTANDARDMAPALG_H

#include <vector>
#include <set>
#include <iostream>

#include "Geometry/ChannelMapAlg.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar{
  namespace geo{
    
    struct XYZPos
    {
      XYZPos(float xPos,
             float yPos,
             float zPos)
      : x(xPos)
      , y(yPos)
      , z(zPos)
      {}
      
      float x;
      float y;
      float z;
      
    };
    
    class ChannelMapStandardAlg : public ChannelMapAlg{
      
    public:
      
      ChannelMapStandardAlg(fhicl::ParameterSet const& p);
      
      void         Initialize(GeometryCore & geo)             override;
      void         Uninitialize()                             override;
      unsigned int Nchannels()                          const override;
      unsigned int NearestChannel(float const* xyz)     const override;
      void         ChannelToPosition(unsigned int chan,
                                     float*       xyz)  const override;
      
    private:
      
      unsigned int        fNumPixelsZ;      ///< number of pixels in the z direction
      unsigned int        fNumPixelsY;      ///< number of pixels in the y direction
      float               fPixelActiveSize; ///< active size of a pixel on an edge
      float               fPixelPitch;      ///< distance between pixel centers
      std::vector<XYZPos> fPixelCenters;    ///< x,y,z center of each pixel
    };
    
    
  } // end namespace geo
} // end namespace gar
#endif // GEO_CHANNELMAPSTANDARDALG_H

