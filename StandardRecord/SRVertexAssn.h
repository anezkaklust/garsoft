#ifndef CAFSRVERTEXASSN_H
#define CAFSRVERTEXASSN_H

#include <cstddef>

namespace caf
{
  class SRVertexAssn
  {
  public:
    SRVertexAssn();

    size_t vtxid;     // Being the vertex which this Assn belongs to
    size_t trkid;
    int    trkend;
  };
}

#endif
