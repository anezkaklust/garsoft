//
//  Vee.cxx
//  garsoft-mrb
//
//  Created by Tom Junk May 3, 2020

#include "ReconstructionDataProducts/Vee.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Vee::Vee(const float *pos,
             const float *covmat,
	     const double  time,
	     const TLorentzVector *FourMomentum)
      : fTime(time)
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();

      size_t icounter = 0;
      for (size_t i=0; i<3; ++i)
	{
	  fVertexPos[i] = pos[i];
	  fFourMomentum[i] = FourMomentum[i];
	  for (size_t j=0; j<3; ++j)
	    {
	      fVertexCov[i][j] = covmat[icounter++];
	    }
	}
    }

    //--------------------------------------------------------------------------
    // empty constructor -- zero all entries -- not sure if we need this -- could be misleading if we ever have one.

    Vee::Vee() 
      : fTime(0)
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();

      for (size_t i=0;i<3;++i)
	{
	  fVertexPos[i] = 0;
	  fFourMomentum[i] *= 0;
	  for (size_t j=0; j<3; ++j)
	    {
	      fVertexCov[i][j] = 0;
	    }
	}
    }



    //--------------------------------------------------------------------------
    // ID number methods
    bool Vee::operator==(const Vee& rhs) const {
        return (this->fIDnumero == rhs.fIDnumero);
    }

    bool Vee::operator!=(const Vee& rhs) const {
        return (this->fIDnumero != rhs.fIDnumero);
    }

    gar::rec::IDNumber Vee::getIDNumber() const {return fIDnumero;}

  } // rec
} // gar
