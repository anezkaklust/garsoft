//
//  Vertex.h
//
//  Created by Tom Junk on 7/11/16.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Vertex_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Vertex_h

#include "RtypesCore.h"
#include <stdint.h>
#include <iostream>

namespace gar {
  namespace rec {

    class Vertex {
    
    public:
      Vertex();
      
      // let the compiler provide the dtor
      
    private:
      
      float                     fPosition[3];    ///< position of vertex
      float                     fCovMat[3][3];   ///< uncertianties on the position
      ULong64_t                 fTime;           ///< timestamp (ROOT's portable version of uint64_t)

#ifndef __GCCXML__
      
    public:
      
      Vertex(float       *pos,
             float       *covmat,
	     ULong64_t   time);
      
      const float*        Position()  const;
      const float*        CovMat()    const;
      ULong64_t           Time()      const;
      
#endif
      
    };
    
    inline const  float* Vertex::Position()  const { return &fPosition[0]; }
    inline const  float* Vertex::CovMat()    const { return &(fCovMat[0][0]); }
    inline ULong64_t     Vertex::Time()      const { return fTime; }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h */
