#ifndef CAFSRTRUTH_H
#define CAFSRTRUTH_H

#include "garsoft/StandardRecord/SRNeutrino.h"
#include "garsoft/StandardRecord/SRGenieParticle.h"
#include "garsoft/StandardRecord/SRParticle.h"
#include "garsoft/StandardRecord/SRSimHitBranch.h"

#include <vector>

namespace caf
{
  class SRTruth
  {
  public:
    SRTruth();

    unsigned int nnu;
    std::vector<SRNeutrino> nu;

    // GENIE particles and MCParticles should both be part of `nu` (ie
    // associated with a particular neutrino). But anatree_module.cc doesn't
    // currently make that association, so place them at this level for now.

    // GENIE particle list from event record
    unsigned int ngpart;
    std::vector<SRGenieParticle> gpart;

    unsigned int npart;
    std::vector<SRParticle> part;

    // Likewise, SimHits have indices and should belong to their particles
    SRSimHitBranch simhit;
  };
}

#endif
