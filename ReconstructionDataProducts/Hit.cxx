//
//  Hit.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "ReconstructionDataProducts/Hit.h"

namespace gar {
  namespace rec {
   
    //--------------------------------------------------------------------------
    Hit::Hit()
    {
      return;
    }

    //--------------------------------------------------------------------------
    Hit::Hit(unsigned int chan,
             float        sig,
             float       *pos,
             float        startT,
             float        endT)
    : fChannel  (chan  )
    , fSignal   (sig   )
    , fStartTime(startT)
    , fEndTime  (endT  )
    {
      
      fPosition[0] = pos[0];
      fPosition[1] = pos[1];
      fPosition[2] = pos[2];
      
      return;
    }

    //--------------------------------------------------------------------------
    std::ostream& operator<< (std::ostream& o, gar::rec::Hit const& h)
    {
      
      o << "Hit on channel "
      << h.Channel()
      << "\n\tposition = ("
      << h.Position()[0]
      << ", "
      << h.Position()[1]
      << ", "
      << h.Position()[2]
      << ")"
      << "\n\tsignal = "
      << h.Signal()
      << "\n\t start time: "
      << h.StartTime()
      << " end time: "
      << h.EndTime();
      
      return o;
    }

  } // rec
} // gar
