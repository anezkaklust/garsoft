//
//  Hit.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef Hit_hpp
#define Hit_hpp


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
      double       fStartTime;   ///< start time of the hit
      double       fEndTime;     ///< end time of the hit
      
#ifndef __GCCXML__
      
    public:
      
      Hit(unsigned int chan,
          float        sig,
          float       *pos);
      
      const float*        Position()  const;
      float        const& Signal()    const;
      unsigned int        Channel()   const;
      double       const& StartTime() const;
      double       const& EndTime()   const;
      
#endif
      
    };
    
    inline unsigned int  Hit::Channel()   const { return fChannel;      }
    inline const  float* Hit::Position()  const { return &fPosition[0]; }
    inline float  const& Hit::Signal()    const { return fSignal;       }
    inline double const& Hit::StartTime() const { return fStartTime;    }
    inline double const& Hit::EndTime()   const { return fEndTime;      }
    
  } // rec
} // gar


#endif /* Hit_hpp */
