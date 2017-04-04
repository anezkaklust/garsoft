//
//  Hit.h
//
//  Created by Brian Rebel on 10/6/16.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h

#include <iostream>

namespace gar {
  namespace rec {

    class Hit {
    
    public:
      Hit();
      
      // let the compiler provide the dtor
      
    private:
      
      unsigned int fChannel;     ///< channel recording this hit
      float        fSignal;      ///< size of the signal for this hit
      float        fPosition[3]; ///< position of the hit
      float        fStartTime;   ///< start time of the hit
      float        fEndTime;     ///< end time of the hit

#ifndef __GCCXML__
      
    public:
      
      Hit(unsigned int chan,
          float        sig,
          float       *pos,
          float        startT,
          float        endT);
      
      const float*        Position()  const;
      float        const& Signal()    const;
      unsigned int        Channel()   const;
      float               StartTime() const;
      float               EndTime()   const;
      
      friend std::ostream& operator << (std::ostream & o, gar::rec::Hit const& h);
      
#endif
      
    };
    
    inline unsigned int  Hit::Channel()   const { return fChannel;      }
    inline const  float* Hit::Position()  const { return &fPosition[0]; }
    inline float  const& Hit::Signal()    const { return fSignal;       }
    inline float         Hit::StartTime() const { return fStartTime;    }
    inline float         Hit::EndTime()   const { return fEndTime;      }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h */
