#ifndef CAFSRTRUTH_H
#define CAFSRTRUTH_H

#include <string>
#include <vector>

namespace caf
{
  class SRTruth
  {
  public:
    SRTruth();

    // MCTruth data.
    // GENIE kinematics computed ignoring Fermi momentum and the off-shellness
    // of the bound nucleon.  Well, that's what the documentation says.  Might
    // not be true!  But to get the right kinematics here, t has to be computed.
    // Mode and InteractionType given in simb::MCNeutrino.h of the nusimdata
    // product.

    std::vector<int>              neutrinoType;
    std::vector<int>              ccnc;
    std::vector<int>              mode;
    std::vector<int>              interactionType;
    std::vector<float>            Q2;
    std::vector<float>            W;
    std::vector<float>            X;
    std::vector<float>            Y;
    std::vector<float>            theta;
    std::vector<float>            T;
    std::vector<float>            MCVertexX;
    std::vector<float>            MCVertexY;
    std::vector<float>            MCVertexZ;
    std::vector<float>            MCnuPx;
    std::vector<float>            MCnuPy;
    std::vector<float>            MCnuPz;

    // GTruth data
    std::vector<int>              gint;
    std::vector<int>              tgtPDG;
    std::vector<float>            weight;
    std::vector<float>            gT;

    //GENIE particle list from event record
    std::vector<int>              nGPart;
    std::vector<int>              GPartIntIdx;
    std::vector<int>              GPartIdx;
    std::vector<int>              GPartPdg;
    std::vector<int>              GPartStatus;
    std::vector<int>              GPartFirstMom;
    std::vector<int>              GPartLastMom;
    std::vector<int>              GPartFirstDaugh;
    std::vector<int>              GPartLastDaugh;
    std::vector<std::string>      GPartName;
    std::vector<float>            GPartPx;
    std::vector<float>            GPartPy;
    std::vector<float>            GPartPz;
    std::vector<float>            GPartE;
    std::vector<float>            GPartMass;

    // MCParticle data
    std::vector<int>              MCTrkID;
    std::vector<int>              MCPDG;
    std::vector<int>              MCMotherIndex;
    std::vector<int>              MCMotherTrkID;
    std::vector<int>              MCPDGMother;
    std::vector<float>            MCPStartX;
    std::vector<float>            MCPStartY;
    std::vector<float>            MCPStartZ;
    std::vector<float>            MCPTime;
    std::vector<float>            MCPStartPX;
    std::vector<float>            MCPStartPY;
    std::vector<float>            MCPStartPZ;
    std::vector<float>            MCPEndX;
    std::vector<float>            MCPEndY;
    std::vector<float>            MCPEndZ;
    std::vector<float>            MCPEndPX;
    std::vector<float>            MCPEndPY;
    std::vector<float>            MCPEndPZ;
    std::vector<std::string>      MCPProc;
    std::vector<std::string>      MCPEndProc;
    std::vector<int>              MCPVertIndex;

    // track trajectory of MCP
    std::vector<float>            TrajMCPX;
    std::vector<float>            TrajMCPY;
    std::vector<float>            TrajMCPZ;
    std::vector<float>            TrajMCPT;
    std::vector<float>            TrajMCPE;
    std::vector<int>              TrajMCPIndex;
    std::vector<int>              TrajMCPTrackID;
    std::vector<float>            TrajMCPPX;
    std::vector<float>            TrajMCPPY;
    std::vector<float>            TrajMCPPZ;

    // sim calo hit data
    unsigned int                  SimnHits;
    std::vector<float>            SimHitX;
    std::vector<float>            SimHitY;
    std::vector<float>            SimHitZ;
    std::vector<float>            SimHitTime;
    std::vector<float>            SimHitEnergy;
    std::vector<int>              SimHitTrackID;
    std::vector<size_t>           SimHitCellID;     // ULong64_t is size_t on 64 bit machines
    float                         SimEnergySum;

    //Muon system sim hits
    unsigned int                  SimnHits_MuID;
    std::vector<float>            SimHitX_MuID;
    std::vector<float>            SimHitY_MuID;
    std::vector<float>            SimHitZ_MuID;
    std::vector<float>            SimHitTime_MuID;
    std::vector<float>            SimHitEnergy_MuID;
    std::vector<int>              SimHitTrackID_MuID;
    std::vector<size_t>           SimHitCellID_MuID;     // ULong64_t is size_t on 64 bit machines
    float                         SimEnergySum_MuID;

  };
}

#endif
