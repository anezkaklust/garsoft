#ifndef CAFSRPARTICLE_H
#define CAFSRPARTICLE_H

#include "StandardRecord/SRVector3D.h"
#include "StandardRecord/SRParticleTrajectoryPoint.h"

#include <string>
#include <vector>

namespace caf
{
  class SRParticle
  {
  public:
    SRParticle();

    int trkid;
    int pdg;
    int momidx;
    int momtrkid;
    int mompdg;
    SRVector3D start;
    float time;
    SRVector3D startp;
    SRVector3D end;
    SRVector3D endp;
    std::string proc;
    std::string endproc;
    int vtxidx;

    unsigned int ntraj;
    std::vector<SRParticleTrajectoryPoint> traj;
  };
}

#endif
