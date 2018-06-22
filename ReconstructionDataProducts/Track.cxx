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
                 float  momentum_beg,
		 float  momentum_end,
                 float *vtx,
                 float *end,
                 float *vtxDir,
                 float *endDir,
		 float chisqfw,
		 float chisqbck,
		 size_t nhits)
    : fLength  (length  )
    , fMomentum_beg(momentum_beg)
    , fMomentum_end(momentum_end)
    , fChisqForward(chisqfw)
    , fChisqBackward(chisqbck)
    {
      fVertex[0] = vtx[0];
      fVertex[1] = vtx[1];
      fVertex[2] = vtx[2];
      
      fVtxDir[0] = vtxDir[0];
      fVtxDir[1] = vtxDir[1];
      fVtxDir[2] = vtxDir[2];

      fEnd[0] = end[0];
      fEnd[1] = end[1];
      fEnd[2] = end[2];
      
      fEndDir[0] = endDir[0];
      fEndDir[1] = endDir[1];
      fEndDir[2] = endDir[2];

      return;
    }
    
  } // rec
} // gar

