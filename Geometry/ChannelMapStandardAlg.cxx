////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapStandardAlg.h"
#include "TMath.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace gar {
  namespace geo{
    
    //----------------------------------------------------------------------------
    ChannelMapStandardAlg::ChannelMapStandardAlg(fhicl::ParameterSet const& p)
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
      
      // get these from the geometry?
      // all dimensions are in cm

      fNumSectors = 18;
      fSectorOffsetAngleDeg = 10.0;
      fPhiSectorWidth = 2.0*TMath::Pi()/fNumSectors;
      fNumPadRowsIROC = 63;
      fNumPadRowsOROCI = 64;
      fNumPadRowsOROCO = 32;
      fPadHeightIROC = 0.75;
      fPadWidthIROC = 0.4;
      fPadWidthOROC = 0.6;
      fPadHeightOROCI = 1.0;
      fPadHeightOROCO = 1.5;
      
      fIROCInnerRadius = 84.1 + fPadHeightIROC;  //  The 84.1 comes from the TDR, which has one extra pad row in it.  Assume it comes off the inside
      fIROCOuterRadius = 132.1;
      fOROCInnerRadius = 134.6;
      fOROCPadHeightChangeRadius = 198.6;
      fOROCOuterRadius = 246.6;

      fFrameWidth = 1.2;
      fSectorGap = 0.3;

      float TsectorH = TMath::Tan(TMath::DegToRad()*(360/(fNumSectors*2)));

      // count up channels
      // get the xyz center for each pixel and fill pad row information -- nominal geometry, negative x side

      for (size_t isector = 0; isector < fNumSectors; ++isector)
	{

	  float rotang = TMath::DegToRad()*( isector*360/fNumSectors + fSectorOffsetAngleDeg );
	  //std::cout << "trot: " << isector << " " << rotang << std::endl;
	  float crot = TMath::Cos(rotang);
	  float srot = TMath::Sin(rotang);

	  size_t ipadacc=0;

	  for (size_t irow = 0; irow < fNumPadRowsIROC; ++irow)
	    {
	      // local coordinates of a sector pointing up.  y is up, z is horizontal.  Rotate it before putting it in to

	      float yloc = fIROCInnerRadius + fPadHeightIROC * (irow + 0.5);
	      float zhalf = yloc * TsectorH - fSectorGap/2.0;
	      size_t numpads = 2*( TMath::Floor(TMath::Abs(zhalf)/fPadWidthIROC) );
	      if (isector == 0) 
		{
		  fNumPadsPerRow.push_back(numpads);
		  fFirstPadInRow.push_back(ipadacc);
		}
	      for (size_t ipad = 0; ipad < numpads; ++ipad)
		{
	          float zloc = ( (float) ipad - (float) numpads/2 + 0.5)*fPadWidthIROC;
	          XYZPos pixpos(-fXPlaneLoc, zloc*crot+yloc*srot, yloc*crot - zloc*srot);
	          fPixelCenters.push_back(pixpos);
		  ipadacc++;
		}
	    }

	  for (size_t irow = 0; irow < fNumPadRowsOROCI; ++irow)
	    {
	      // local coordinates of a sector pointing up.  y is up, z is horizontal.  Rotate it before putting it in to

	      float yloc = fOROCInnerRadius + fPadHeightOROCI * (irow + 0.5);
	      float zhalf = yloc * TsectorH - fSectorGap/2.0;
	      size_t numpads = 2*( TMath::Floor(TMath::Abs(zhalf)/fPadWidthOROC) );
	      if (isector == 0) 
		{
		  fNumPadsPerRow.push_back(numpads);
		  fFirstPadInRow.push_back(ipadacc);
		}
	      for (size_t ipad = 0; ipad < numpads; ++ipad)
		{
	          float zloc = ( (float) ipad - (float) numpads/2 + 0.5)*fPadWidthOROC;
	          XYZPos pixpos(-fXPlaneLoc, zloc*crot+yloc*srot, yloc*crot - zloc*srot);
	          fPixelCenters.push_back(pixpos);
		  ipadacc++;
		}
	    }

	  for (size_t irow = 0; irow < fNumPadRowsOROCO; ++irow)
	    {
	      // local coordinates of a sector pointing up.  y is up, z is horizontal.  Rotate it before putting it in to

	      float yloc = fOROCPadHeightChangeRadius + fPadHeightOROCO * (irow + 0.5);
	      float zhalf = yloc * TsectorH - fSectorGap/2.0;
	      size_t numpads = 2*( TMath::Floor(TMath::Abs(zhalf)/fPadWidthOROC) );
	      if (isector == 0) 
		{
		  fNumPadsPerRow.push_back(numpads);
		  fFirstPadInRow.push_back(ipadacc);
		}
	      for (size_t ipad = 0; ipad < numpads; ++ipad)
		{
	          float zloc = ( (float) ipad - (float) numpads/2 + 0.5)*fPadWidthOROC;
	          XYZPos pixpos(-fXPlaneLoc, zloc*crot+yloc*srot, yloc*crot - zloc*srot);
	          fPixelCenters.push_back(pixpos);
		  ipadacc++;
		}
	    }
	  if (isector == 0) fNumChansPerSector = ipadacc;
	}

      size_t numpixside = fPixelCenters.size();

      // duplicate the yz geometry for the opposite side -- positive x

      for (size_t ipix = 0; ipix < numpixside; ++ipix)
	{
	  XYZPos pixpos(fXPlaneLoc,fPixelCenters[ipix].y,fPixelCenters[ipix].z);
	  fPixelCenters.push_back(pixpos);
	  //std::cout << "trjpix " << fPixelCenters[ipix].z << " " << fPixelCenters[ipix].y << std::endl;
	}

      
      LOG_DEBUG("ChannelMapStandardAlg")
      << "There are "
      << fPixelCenters.size()
      << " pixels and each sector has "
      << fNumPadsPerRow.size()
      << " pad rows";
      

      return;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Uninitialize()
    {
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::Nchannels() const
    {
      return (fPixelCenters.size());
    }
    
    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::NearestChannel(float const* xyz) const
    {

      
      float phi = TMath::ATan2(xyz[1],xyz[2]);
      if (phi<0) phi += 2.0*TMath::Pi();
      float phisc = phi/( fPhiSectorWidth );
      size_t isector = TMath::Floor(phisc); // assumes the sector boundary is at phi=0

      // rotate this back down to a single sector -- to do -- test this!

      float rotang = TMath::DegToRad()*( isector*360/fNumSectors + fSectorOffsetAngleDeg );
      float crot = TMath::Cos(rotang);
      float srot = TMath::Sin(rotang);
      float zrot =    xyz[2]*crot + xyz[1]*srot;
      float yrot =  - xyz[2]*srot + xyz[1]*crot;

      // zrot is r, and yrot=0 is centered on the midplane

      float rtmp=zrot;
      if (zrot < fIROCInnerRadius) rtmp = fIROCInnerRadius;  // To FIX -- fill in the plug -- this puts everything on the inner pad row
      if (zrot > fOROCOuterRadius) rtmp = fOROCOuterRadius;
      if (zrot > fIROCOuterRadius && zrot < (fIROCOuterRadius + fOROCInnerRadius)/2.0) rtmp = fIROCOuterRadius;
      if (zrot < fOROCInnerRadius && zrot >= (fIROCOuterRadius + fOROCInnerRadius)/2.0) rtmp = fOROCInnerRadius;

      float padwidthloc = fPadWidthOROC;
      size_t irow = 0;
      if (rtmp <= fIROCOuterRadius)
	{
	  irow = TMath::Floor( (rtmp-fIROCInnerRadius)/fPadHeightIROC );
	  padwidthloc = fPadWidthIROC;
	}
      else if (rtmp < fOROCPadHeightChangeRadius)
	{
	  irow =  TMath::Floor( (rtmp-fOROCInnerRadius)/fPadHeightOROCI ) + fNumPadRowsIROC;
	  padwidthloc = fPadWidthOROC;
	}
      else
	{
	  irow = TMath::Floor( (rtmp-fOROCPadHeightChangeRadius)/fPadHeightOROCO ) + fNumPadRowsIROC + fNumPadRowsOROCI;
	  padwidthloc = fPadWidthOROC;
	}
      size_t totpadrows = fNumPadRowsIROC + fNumPadRowsOROCI + fNumPadRowsOROCO;
      if (irow >= totpadrows) irow = totpadrows-1; 

      size_t ichansector = fFirstPadInRow[irow] + TMath::Floor(yrot/padwidthloc) + fNumPadsPerRow[irow]/2;

      size_t ichan = ichansector + fNumChansPerSector * isector;

      if (xyz[0] > 0) ichan += fNumSectors*fNumChansPerSector;  // the opposite side of the TPC.

      //   throw cet::exception("NearestChannel")
      //<< "y position: "
      //<< xyz[1]
      //<< " or z position: "
      //<< xyz[2]
      //<< " is out of bounds.";

      //if (ichan>1000000) 
      //	{
      //  std::cout << "Problem Channel ID: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << ichan << std::endl;
      //}
      return ichan;
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
