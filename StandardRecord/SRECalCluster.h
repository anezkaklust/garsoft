#ifndef CAFSRECALCLUSTER_H
#define CAFSRECALCLUSTER_H

#include "garsoft/StandardRecord/SRVector3D.h"

#include <vector>

namespace caf
{
  class SRECalCluster
  {
  public:
    SRECalCluster();

    unsigned int                  nCluster;
    std::vector<size_t>           id;
    std::vector<unsigned int>     nhits;
    std::vector<float>            E;
    std::vector<float>            t;
    std::vector<float>            TimeDiffFirstLast;
    std::vector<float>            x;
    std::vector<float>            y;
    std::vector<float>            z;
    std::vector<float>            theta;
    std::vector<float>            phi;
    std::vector<float>            pid;
    // std::vector<float> shape;
    std::vector<SRVector3D>       mainAxis;
    std::vector<int>              MCindex;          // Branch index (NOT the GEANT track ID) of MCPartice
    std::vector<float>            MCfrac;           // that best matches & fraction of ionization therefrom

    // ECAL cluster to track association info
    std::vector<size_t>           assn_clustid; // Being the cluster which this Assn belongs to
    std::vector<size_t>           assn_trkid;  // The rec::TrackEnd (see Track.h) that extrapolated to cluster
    std::vector<int>              assn_trkend;
  };
}

#endif
