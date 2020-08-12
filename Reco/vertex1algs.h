#ifndef GARSOFT_RECO_VERTEX2ALGS_H
#define GARSOFT_RECO_VERTEX2ALGS_H

#include <vector>
#include <iostream>
#include <cstddef>
#include <TMath.h>
#include <TVector3.h>
#include <TVectorF.h>
#include <ReconstructionDataProducts/Track.h>
#include <Reco/TrackPar.h>

namespace gar {
  namespace rec {

    // fit an arbitrary number of tracks to a single vertex
    // returns 0 if success, 1 if failure
    // covmat is 3x3

    int fitVertex(std::vector<TrackPar> &tracks, 
                  std::vector<float> &xyz, 
                  float &chisquared, 
                  std::vector< std::vector<float> > &covmat, 
                  double &time,
                  std::vector<TrackEnd> usebeg);

      }
}
#endif
