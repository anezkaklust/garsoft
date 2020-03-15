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
      float        fSignal;      ///< size of the signal for this hit  (integral of ADC values)
      float        fPosition[3]; ///< position of the hit
      float        fTime;        ///< time of hit charge arrival at the readout plane (ticks)
      float        fStartTime;   ///< start time of the hit (ticks)
      float        fEndTime;     ///< end time of the hit (ticks)
      float        fRMS;         ///< Hit width calculated with RMS (in ticks)

#ifndef __GCCXML__
      
    public:
      
      Hit(unsigned int chan,
          float        sig,
          float       *pos,
          float        startT,
          float        endT,
      float        Time,
      float        RMS);
      
      const float*        Position()  const;
      float        const& Signal()    const;
      unsigned int        Channel()   const;
      float               StartTime() const;
      float               EndTime()   const;
      float               Time()      const;
      float               RMS()       const;
      
      void operator += (gar::rec::Hit const& h);

      bool operator <  (gar::rec::Hit const& h);
      
      friend std::ostream& operator << (std::ostream & o, gar::rec::Hit const& h);
      
#endif
      
    };
    
    inline unsigned int  Hit::Channel()   const { return fChannel;      }
    inline const  float* Hit::Position()  const { return &fPosition[0]; }
    inline float  const& Hit::Signal()    const { return fSignal;       }
    inline float         Hit::StartTime() const { return fStartTime;    }
    inline float         Hit::EndTime()   const { return fEndTime;      }
    inline float         Hit::RMS()       const { return fRMS;          }
    inline float         Hit::Time()      const { return fTime;         }    
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h */
