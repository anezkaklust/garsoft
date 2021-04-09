#ifndef CAFSRCLUSTERASSN_H
#define CAFSRCLUSTERASSN_H

#include <cstddef>

namespace caf
{
  class SRClusterAssn
  {
  public:
    SRClusterAssn();

    size_t clustid;
    size_t trkid;
    int    trkend;
  };
}

#endif
