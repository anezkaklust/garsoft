#ifndef CAFSRNEUTRINO_H
#define CAFSRNEUTRINO_H

#include "StandardRecord/SRVector3D.h"

namespace caf
{
  class SRNeutrino
  {
  public:
    SRNeutrino();

    // MCTruth data.
    // GENIE kinematics computed ignoring Fermi momentum and the off-shellness
    // of the bound nucleon.  Well, that's what the documentation says.  Might
    // not be true!  But to get the right kinematics here, t has to be computed.
    // Mode and InteractionType given in simb::MCNeutrino.h of the nusimdata
    // product.

    int   pdg;
    int   ccnc;
    int   mode;
    int   interactionType;
    float Q2;
    float W;
    float X;
    float Y;
    float theta;
    float T;
    SRVector3D vtx;
    SRVector3D p;

    // GTruth data
    int    gint;
    int    tgtPDG;
    float  weight;
    float  gT;
  };
}

#endif
