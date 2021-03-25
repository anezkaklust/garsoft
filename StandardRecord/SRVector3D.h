#ifndef CAFSRVECTOR3D_H
#define CAFSRVECTOR3D_H

namespace caf
{
  class SRVector3D
  {
  public:
    SRVector3D();
    SRVector3D(float x, float y, float z);

    float x, y, z;
  };
}

#endif
