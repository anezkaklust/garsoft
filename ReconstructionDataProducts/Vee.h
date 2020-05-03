//
//  Vee.h
//
//  Created by Tom Junk on May 3, 2020.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Vee_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Vee_h

#include <iostream>
#include <TLorentzVector.h>
#include "IDNumberGen.h"

namespace gar {
  namespace rec {

    class Vee {
    
    public:
      Vee();
      
      // let the compiler provide the dtor

      // three hypotheses for mass and momentum

      typedef enum { Kshort=0, Lambda1=1, Lambda2=2 } hypothesis_t;

    private:
      static gar::rec::IDNumber const FirstNumber = 800000;
      gar::rec::IDNumber fIDnumero;
      
      float        fVertexPos[3];     ///< position of the vertex, in cm
      float        fVertexCov[3][3];  ///< uncertainties on vertex position, in cm
      double       fTime;
      float        fChisq;
      TLorentzVector fFourMomentum[3]; // four-momentum of vee, one for each hypothesis, in GeV


#ifndef __GCCXML__
      
    public:
      
      Vee(const float          *pos,
          const float          *vertexcov,
	  const float          chisq,
	  const double         time,
	  const TLorentzVector *fourmomentum);
      
      const float*           Position()   const;
      const float*           VertexCov()  const;
      double                 Time()       const;
      const TLorentzVector&  FourMomentum(const size_t i) const;
      float                  Chisq()      const;

      bool operator==(const Vee& rhs) const;
      bool operator!=(const Vee& rhs) const;
      gar::rec::IDNumber getIDNumber() const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::Vee const& h);
      
#endif
      
    };
    
    inline const  float* Vee::Position()    const { return &fVertexPos[0]; }
    inline const  float* Vee::VertexCov()   const { return &(fVertexCov[0][0]); }
    inline double        Vee::Time()        const { return fTime;        }
    inline float         Vee::Chisq()       const { return fChisq;        }
    inline const TLorentzVector& Vee::FourMomentum(const size_t i)  const { return fFourMomentum[i]; }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Vee_h */
