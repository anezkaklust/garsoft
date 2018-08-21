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
      {};

      XYZPos() {x=0; y=0; z=0; };  

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

      size_t              fNumSectors;                 ///<   Number of sectors -- should be 18
      float               fSectorOffsetAngleDeg;       ///<   Angle to rotate to the middle of the first sector -- should be 10 degrees
      float               fPhiSectorWidth;             ///<   width of a sector in phi (in radians)

      size_t              fNumPadRowsIROC;             ///<   Number of pad rows in the inner ROC -- 64 (TDR) or 63 (ALICE code)
      size_t              fNumPadRowsOROCI;            ///<   Number of small-pitch pad rows in the outer ROC 
      size_t              fNumPadRowsOROCO;            ///<   Number of large-pitch pad rows in the outer ROC
      float               fPadHeightIROC;              ///<   Pad height in the inner ROC (cm)
      float               fPadWidthIROC;               ///<   Pad width in the inner ROC (cm)
      float               fPadHeightOROCI;             ///<   Pad height in the outer ROC inner part (cm)
      float               fPadHeightOROCO;             ///<   Pad height in the outer ROC outer part (cm)
      float               fPadWidthOROC;               ///<   Pad width in the OROC (assumed same for both sections)

      float               fIROCInnerRadius;            ///<   Radius from the beam in cm along the midline of the sector to the inner IROC row inner edge
      float               fIROCOuterRadius;            ///<   Radius from the beam in cm along the midline of the sector to the outer IROC row outer edge
      float               fOROCInnerRadius;            ///<   Radius from the beam in cm along the midline of the sector to the inner OROC row inner edge
      float               fOROCPadHeightChangeRadius;  ///<   Radius from the beam in cm along the midline of the sector to the OROC pad height change
      float               fOROCOuterRadius;            ///<   Radius from the beam in cm along the midline of the sector to the outer OROC row outer edge

      float               fFrameWidth;                 ///< wire-fixation sides between sectors  in cm (1.2 cm)
      float               fSectorGap;                  ///< Gap between readout sectors  (0.3 cm)

      float               fXPlaneLoc;                  ///< Location of pixel plane in X (only positive.  Assume other one is at -X)

      std::vector<size_t> fNumPadsPerRow;              ///< indexed by "global" pad row number for a single sector
      std::vector<size_t> fFirstPadInRow;              ///< indexed by "global" pad row number for a single sector

      std::vector<XYZPos> fPixelCenters;               ///< pixel centers (in cm) -- for the entire detector

      size_t              fNumChansPerSector;          ///< Number of TPC pad channels per sector

      // variables for the hole filler grid

      std::vector<size_t> fCenterNumPadsPerRow;        ///< pads per row for the center hole filler
      std::vector<size_t> fCenterFirstPadInRow;        ///< first pad in row for center hole filler

      float               fCenterPadWidth;             ///< Width of square pads in center hole filler
      size_t              fNumChansCenter;             ///< Number of channels in center hole filler

      void                CheckPositions();            ///< Method to check consistency of NearestChannel and ChannelToPosition

      XYZPos              fTPCCenter;                  ///< Location of the center of the TPC
    };
    
    
  } // end namespace geo
} // end namespace gar
#endif // GEO_CHANNELMAPSTANDARDALG_H

