#ifndef CAFSRPARTICLETRAJECTORYPOINT_H
#define CAFSRPARTICLETRAJECTORYPOINT_H

namespace caf
{
  class SRParticleTrajectoryPoint
  {
  public:
    SRParticleTrajectoryPoint();

    float x, y, z;
    float t;
    float E;
    float px, py, pz;

  };
}

#endif
