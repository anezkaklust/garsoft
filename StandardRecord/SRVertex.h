#ifndef CAFSRVERTEX_H
#define CAFSRVERTEX_H

#include <vector>

namespace caf
{
  class SRVertex
  {
  public:
    SRVertex();

    // vertex branches
    std::vector<size_t>  id;
    std::vector<float>   x;
    std::vector<float>   y;
    std::vector<float>   z;
    std::vector<size_t>  t;
    std::vector<int>     ntrks;
    std::vector<int>     Q;

    std::vector<size_t>   assn_vtxid;     // Being the vertex which this Assn belongs to
    std::vector<size_t>   assn_trkid;
    std::vector<int>      assn_trkend;
  };
}

#endif
