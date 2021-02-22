#ifndef GAR_FSPARTICLE_H
#define GAR_FSPARTICLE_H

#include "nusimdata/SimulationBase/MCParticle.h"
#include "garana/DataProducts/FSParticle.h"

namespace gar{

 class GarFSParticle {

  public:

      GarFSParticle() {}
      GarFSParticle(simb::MCParticle mcp);

      garana::FSParticle GetGaranaObj() { return *fFSP; };

  private:

      garana::FSParticle* fFSP = nullptr;
 
 };//class
}//namespace

#endif
