#ifndef CAFSRVEEASSN_H
#define CAFSRVEEASSN_H

#include <cstddef>

namespace caf
{
  class SRVeeAssn
  {
  public:
    SRVeeAssn();

    size_t veeid;
    size_t trkid;
    int    trkend;
  };
}

#endif
