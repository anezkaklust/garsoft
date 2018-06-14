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
      float fMomentum;  ///< momentum of the track at the vertex
      float fVertex[3]; ///< track vertex position
      float fEnd[3];    ///< track end    position
      float fVtxDir[3]; ///< track vertex direction
      float fEndDir[3]; ///< track end    direction
      
#ifndef __GCCXML__
      
    public:
      
      Track(float  length,
            float  momentum,
            float *vtx,
            float *end,
            float *vtxDir,
            float *endDir);
      
      const float* Vertex()   const;
      const float* End()      const;
      const float* VtxDir()   const;
      const float* EndDir()   const;
      float const& Length()   const;
      float const& Momentum() const;
      
#endif
      
    };
    
    inline const float* Track::Vertex()   const { return fVertex;   }
    inline const float* Track::End()      const { return fEnd;      }
    inline const float* Track::VtxDir()   const { return fVtxDir;   }
    inline const float* Track::EndDir()   const { return fEndDir;   }
    inline float const& Track::Length()   const { return fLength;   }
    inline float const& Track::Momentum() const { return fMomentum; }
    
  } // rec
} // gar

#endif /* Track_hpp */
