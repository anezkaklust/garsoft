//
//  Track.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

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
      
      float fLength;    ///< length of the track
      float fMomentum_beg;  ///< momentum of the track at the vertex
      float fMomentum_end;  ///< momentum of the track at the end
      float fVertex[3]; ///< track vertex position
      float fEnd[3];    ///< track end    position
      float fVtxDir[3]; ///< track vertex direction
      float fEndDir[3]; ///< track end    direction
      float fChisqForward; ///< chisquared forward fit
      float fChisqBackward; ///< chisquared backward fit
      size_t fNHits;        ///< number of hits

#ifndef __GCCXML__
      
    public:
      
      Track(float  length,
            float  momentum_beg,
            float  momentum_end,
            float *vtx,
            float *end,
            float *vtxDir,
            float *endDir,
	    float chisqfw,
	    float chisqbck,
	    size_t nhits);
      
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
    
  } // rec
} // gar

#endif /* Track_hpp */
