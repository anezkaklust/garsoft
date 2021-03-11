#ifndef CAFSRTRACK_H
#define CAFSRTRACK_H

#include <vector>

namespace caf
{
  class SRTrack
  {
  public:
    SRTrack();

    // track data
    std::vector<size_t>           id;
    std::vector<float>            StartX;
    std::vector<float>            StartY;
    std::vector<float>            StartZ;
    std::vector<float>            StartPX;
    std::vector<float>            StartPY;
    std::vector<float>            StartPZ;
    std::vector<int>              StartQ;

    std::vector<float>            EndX;
    std::vector<float>            EndY;
    std::vector<float>            EndZ;
    std::vector<float>            EndPX;
    std::vector<float>            EndPY;
    std::vector<float>            EndPZ;
    std::vector<int>              EndQ;

    std::vector<float>            LenF;         // from foward fit, from the Beg end to the End end
    std::vector<float>            LenB;
    std::vector<float>            Chi2F;        // from foward fit, from the Beg end to the End end
    std::vector<float>            Chi2B;
    std::vector<int>              NTPCClustersOnTrack;
    std::vector<float>            AvgIonF;      // from foward fit, from the Beg end to the End end
    std::vector<float>            AvgIonB;

    std::vector<int>              PIDF;
    std::vector<float>            PIDProbF;
    std::vector<int>              PIDB;
    std::vector<float>            PIDProbB;
    std::vector<int>              MCindex;      // Branch index (NOT the GEANT track ID) of MCPartice
    std::vector<float>            MCfrac;       // that best matchs & fraction of ionization therefrom

    //Trajectory
    std::vector<float>            TrajectoryFWDX; //forward
    std::vector<float>            TrajectoryFWDY; //forward
    std::vector<float>            TrajectoryFWDZ; //forward
    std::vector<int>              TrajectoryFWDID;

    std::vector<float>            TrajectoryBWDX; //backward
    std::vector<float>            TrajectoryBWDY; //backward
    std::vector<float>            TrajectoryBWDZ; //backward
    std::vector<int>              TrajectoryBWDID;
  };
}

#endif
