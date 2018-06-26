//
//  Track.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  original Copyright © 2016 Brian Rebel. All rights reserved.
//  Additions by Tom Junk, 2018

#ifndef Track_hpp
#define Track_hpp

#include <stdio.h>


namespace gar {
  namespace rec {
    
    class Track {
      
    public:
      Track();
      
        // let the compiler provide the dtor
      
    private:
      
      float fLength;    ///< length of the track in cm
      float fMomentum_beg;  ///< momentum of the track at the vertex in GeV/c
      float fMomentum_end;  ///< momentum of the track at the end in GeV/c
      float fVertex[3]; ///< track vertex position in cm -- == "beginning" of track.  Arbitrary choice made by patrec
      float fEnd[3];    ///< track end    position in cm
      float fVtxDir[3]; ///< track vertex direction in cm
      float fEndDir[3]; ///< track end    direction in cm
      float fChisqForward; ///< chisquared forward fit
      float fChisqBackward; ///< chisquared backward fit
      size_t fNHits;        ///< number of hits

      // use the x from fVertex and fEnd to specify the independent variable -- no need to store them twice

      float fTrackParBeg[5]; ///< Track parameters at beginning of track y, z, curvature, phi, slope  -- 5-parameter track  (cm, cm, cm-1, radians, dy,z/dx)
      float fTrackParEnd[5]; ///< Track parameters at end of track y, z, curvature, phi, slope  -- 5-parameter track  (cm, cm, cm-1, radians, dy,z/dx)

      float fCovMatBeg[15]; ///< covariance matrix at beginning of track -- packed in a 1D array, assuming symmetry
      float fCovMatEnd[15]; ///< covariance matrix at end of track

#ifndef __GCCXML__
      
    public:

      // Cartesian coordinate constructor -- to be used with the MC Cheater, where
      // the uncertainties are not specified and are set to zero in the constructor.
      
      Track(float  length,
            float  momentum_beg,
            float  momentum_end,
            float *vtx,
            float *end,
            float *vtxDir,
            float *endDir,
	    size_t nhits,
	    int    charge); // units of e.  Need this to convert momentum into  curvature

      // constructor to fill after the fits -- including covariance matrix

      Track(float length,
	    size_t nhits,
	    float xbeg,           // x location at beginning of track in cm
	    float *trackparbeg,   // y, z, curvature, phi, slope  -- 5-parameter track  (cm, cm, cm-1, radians, dy,z/dx)
	    float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
	    float chisqforward,   // chisquared of forwards fit
	    float xend,           // x location at end of track
	    float *trackparend,   // y, z, curvature, phi, slope  -- 5-parameter track (cm, cm, cm-1, radians, dy,z/dx)
	    float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
	    float chisqbackward); // chisquared of forwards fit
	    
      
      const float* Vertex()   const;
      const float* End()      const;
      const float* VtxDir()   const;
      const float* EndDir()   const;
      float const& Length()   const;
      float const& Momentum_beg() const;
      float const& Momentum_end() const;
      float const& ChisqForward() const;
      float const& ChisqBackward() const;
      size_t const& NHits() const;
      const float* TrackParBeg() const;
      const float* TrackParEnd() const;
      void CovMatBegSymmetric(float *cmb) const; // fill a 5x5 array of floats owned by the caller with the values of the covariance matrix at the beginning of the track
      void CovMatEndSymmetric(float *cme) const; // fill a 5x5 array of floats owned by the caller with the values of the covariance matrix at the beginning of the track
      const float* CovMatBegPacked() const;     
      const float* CovMatEndPacked() const;
      int ChargeBeg();   // just returns +1 or -1 depending on the sign of the curvature at the track beginning point
      int ChargeEnd();   // just returns +1 or -1 depending on the sign of the curvature at the track ending point
      //   the charge assumes the track started at the beginning or the end above, and is expected to change sign
      // depending on which way the track is hypothesized to go
      
#endif
      
    };
    
    inline const float* Track::Vertex()   const { return fVertex;   }
    inline const float* Track::End()      const { return fEnd;      }
    inline const float* Track::VtxDir()   const { return fVtxDir;   }
    inline const float* Track::EndDir()   const { return fEndDir;   }
    inline float const& Track::Length()   const { return fLength;   }
    inline float const& Track::Momentum_beg() const { return fMomentum_beg; }
    inline float const& Track::Momentum_end() const { return fMomentum_end; }
    inline float const& Track::ChisqForward() const { return fChisqForward; }
    inline float const& Track::ChisqBackward() const { return fChisqBackward; }
    inline size_t const& Track::NHits() const { return fNHits; }
    inline const float* Track::TrackParBeg() const { return fTrackParBeg; }
    inline const float* Track::TrackParEnd() const { return fTrackParEnd; }
    inline const float* Track::CovMatBegPacked() const { return fCovMatBeg; }
    inline const float* Track::CovMatEndPacked() const { return fCovMatEnd; }

  } // rec
} // gar

#endif /* Track_hpp */
