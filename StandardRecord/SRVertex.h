#ifndef CAFSRVERTEX_H
#define CAFSRVERTEX_H

#include <cstddef>

namespace caf
{
  class SRVertex
  {
  public:
    SRVertex();

    // vertex branches
    size_t  id;
    float   x, y, z;
    size_t  t;
    int     ntrks;
    int     Q;
  };
}

#endif
