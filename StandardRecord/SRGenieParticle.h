#ifndef CAFSRGENIEPARTICLE_H
#define CAFSRGENIEPARTICLE_H

#include "garsoft/StandardRecord/SRVector3D.h"

#include <string>

namespace caf
{
  class SRGenieParticle
  {
  public:
    SRGenieParticle();

    int intidx;
    int idx;
    int pdg;
    int status;
    int firstmom;
    int lastmom;
    int firstdaugh;
    int lastdaugh;
    std::string name;
    SRVector3D p;
    float E;
    float mass;
  };
}

#endif
