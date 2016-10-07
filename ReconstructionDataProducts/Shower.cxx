//
//  Shower.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "ReconstructionDataProducts/Shower.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Shower::Shower()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    Shower::Shower(float  energy,
                   float *vtx,
                   float *vtxDir)
    : fEnergy  (energy  )
    , fVertex  (vtx     )
    , fVtxDir  (vtxDir  )
    {
      return;
    }
    
  } // rec
} // gar
