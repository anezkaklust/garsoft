#ifndef CAFSRVEE_H
#define CAFSRVEE_H

#include <vector>

namespace caf
{
  class SRVee
  {
  public:
    SRVee();

    // Vee branches
    std::vector<size_t>           id;
    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<size_t>           t;
    // TODO introduce SRVector3D
    std::vector<float>            PXKpipi;
    std::vector<float>            PYKpipi;
    std::vector<float>            PZKpipi;
    std::vector<float>            EKpipi;
    std::vector<float>            MKpipi;
    std::vector<float>            PXLppi;
    std::vector<float>            PYLppi;
    std::vector<float>            PZLppi;
    std::vector<float>            ELppi;
    std::vector<float>            MLppi;
    std::vector<float>            PXLpip;
    std::vector<float>            PYLpip;
    std::vector<float>            PZLpip;
    std::vector<float>            ELpip;
    std::vector<float>            MLpip;
    std::vector<size_t>           assn_veeid;    // Being the Vee which this Assn belongs to
    std::vector<size_t>           assn_trkid;
    std::vector<int>              assn_trkend;
  };
}

#endif
