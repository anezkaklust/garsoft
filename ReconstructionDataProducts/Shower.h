//
//  Shower.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef Shower_hpp
#define Shower_hpp

#include <stdio.h>


namespace gar {
  namespace rec {
    
    class Shower {
      
    public:
      Shower();
      
        // let the compiler provide the dtor
      
    private:
      
      float fEnergy;    ///< energy of the shower
      float fVertex[3]; ///< Shower vertex position
      float fVtxDir[3]; ///< Shower vertex direction
      
#ifndef __GCCXML__
      
    public:
      
      Shower(float  energy,
             float *vtx,
             float *vtxDir);
      
      const float* Vertex() const;
      const float* VtxDir() const;
      float const& Energy() const;
      
#endif
      
    };
    
    inline const float* Shower::Vertex() const { return fVertex; }
    inline const float* Shower::VtxDir() const { return fVtxDir; }
    inline float float& Shower::Energy() const { return fEnergy; }
    
  } // rec
} // gar


#endif /* Shower_hpp */

