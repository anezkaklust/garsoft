#ifndef GAR_GTRUTH_H
#define GAR_GTRUTH_H

#include "nusimdata/SimulationBase/GTruth.h"
#include "garana/DataProducts/GTruth.h"

using garana::GTruth;

namespace gar{
 namespace adp{
    class GarGTruth : public garana::GTruth{

      public:
        GarGTruth();
        GarGTruth(const simb::GTruth& gt);
        
        GTruth GetGaranaObj();

      private:

        void FillGT(const simb::GTruth& gt);

        GTruth fGT;
    };
 }
}

#endif
