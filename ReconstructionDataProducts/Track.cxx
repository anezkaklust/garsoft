//
//  Track.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "ReconstructionDataProducts/Track.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Track::Track()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    Track::Track(float  length,
                 float  momentum,
                 float *vtx,
                 float *end,
                 float *vtxDir,
                 float *endDir)
    : fLength  (length  )
    , fMomentum(momentum)
    , fVertex  (vtx     )
    , fEnd     (end     )
    , fVtxDir  (vtxDir  )
    , fEndDir  (endDir  )
    {
      return;
    }
    
  } // rec
} // gar

