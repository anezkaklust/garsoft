#ifndef CAFSTANDARDRECORD_H
#define CAFSTANDARDRECORD_H

#include <string>
#include <vector>

namespace caf
{
  class StandardRecord
  {
  public:
    StandardRecord();

    // global event info
    int   event;      ///< number of the event being processed
    int   run;        ///< number of the run being processed
    int   subrun;     ///< number of the sub-run being processed

    float TPC_X;      ///< center of TPC stored as per-event & compressed by root
    float TPC_Y;
    float TPC_Z;

    // MCTruth data.
    // GENIE kinematics computed ignoring Fermi momentum and the off-shellness
    // of the bound nucleon.  Well, that's what the documentation says.  Might
    // not be true!  But to get the right kinematics here, t has to be computed.
    // Mode and InteractionType given in simb::MCNeutrino.h of the nusimdata
    // product.
    // Use Rtypes.h here, as these data get used by root

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

    // Hit data
    std::vector<float>            HitX;
    std::vector<float>            HitY;
    std::vector<float>            HitZ;
    std::vector<float>            HitSig;
    std::vector<float>            HitRMS;
    std::vector<unsigned int>     HitChan;

    // TPCCluster data
    std::vector<float>            TPCClusterX;
    std::vector<float>            TPCClusterY;
    std::vector<float>            TPCClusterZ;
    std::vector<float>            TPCClusterSig;
    std::vector<float>            TPCClusterRMS;
    std::vector<size_t>           TPCClusterTrkIDNumber;
    std::vector<float>            TPCClusterCovXX;
    std::vector<float>            TPCClusterCovXY;
    std::vector<float>            TPCClusterCovXZ;
    std::vector<float>            TPCClusterCovYY;
    std::vector<float>            TPCClusterCovYZ;
    std::vector<float>            TPCClusterCovZZ;

    // track data
    std::vector<size_t>           TrackIDNumber;
    std::vector<float>            TrackStartX;
    std::vector<float>            TrackStartY;
    std::vector<float>            TrackStartZ;
    std::vector<float>            TrackStartPX;
    std::vector<float>            TrackStartPY;
    std::vector<float>            TrackStartPZ;
    std::vector<int>              TrackStartQ;

    std::vector<float>            TrackEndX;
    std::vector<float>            TrackEndY;
    std::vector<float>            TrackEndZ;
    std::vector<float>            TrackEndPX;
    std::vector<float>            TrackEndPY;
    std::vector<float>            TrackEndPZ;
    std::vector<int>              TrackEndQ;

    std::vector<float>            TrackLenF;         // from foward fit, from the Beg end to the End end
    std::vector<float>            TrackLenB;
    std::vector<float>            TrackChi2F;        // from foward fit, from the Beg end to the End end
    std::vector<float>            TrackChi2B;
    std::vector<int>              NTPCClustersOnTrack;
    std::vector<float>            TrackAvgIonF;      // from foward fit, from the Beg end to the End end
    std::vector<float>            TrackAvgIonB;

    std::vector<int>              TrackPIDF;
    std::vector<float>            TrackPIDProbF;
    std::vector<int>              TrackPIDB;
    std::vector<float>            TrackPIDProbB;
    std::vector<int>              TrackMCindex;      // Branch index (NOT the GEANT track ID) of MCPartice
    std::vector<float>            TrackMCfrac;       // that best matchs & fraction of ionization therefrom

    //TrackTrajectory
    std::vector<float>            TrackTrajectoryFWDX; //forward
    std::vector<float>            TrackTrajectoryFWDY; //forward
    std::vector<float>            TrackTrajectoryFWDZ; //forward
    std::vector<int>              TrackTrajectoryFWDID;

    std::vector<float>            TrackTrajectoryBWDX; //backward
    std::vector<float>            TrackTrajectoryBWDY; //backward
    std::vector<float>            TrackTrajectoryBWDZ; //backward
    std::vector<int>              TrackTrajectoryBWDID;

    // vertex branches
    std::vector<size_t>           VertexIDNumber;
    std::vector<float>            VertexX;
    std::vector<float>            VertexY;
    std::vector<float>            VertexZ;
    std::vector<size_t>           VertexT;
    std::vector<int>              VertexN;
    std::vector<int>              VertexQ;

    std::vector<size_t>           VTAssn_VertIDNumber;     // Being the vertex which this Assn belongs to
    std::vector<size_t>           VTAssn_TrackIDNumber;
    std::vector<int>              VTAssn_TrackEnd;

    // Vee branches
    std::vector<size_t>           VeeIDNumber;
    std::vector<float>            VeeX;
    std::vector<float>            VeeY;
    std::vector<float>            VeeZ;
    std::vector<size_t>           VeeT;
    std::vector<float>            VeePXKpipi;
    std::vector<float>            VeePYKpipi;
    std::vector<float>            VeePZKpipi;
    std::vector<float>            VeeEKpipi;
    std::vector<float>            VeeMKpipi;
    std::vector<float>            VeePXLppi;
    std::vector<float>            VeePYLppi;
    std::vector<float>            VeePZLppi;
    std::vector<float>            VeeELppi;
    std::vector<float>            VeeMLppi;
    std::vector<float>            VeePXLpip;
    std::vector<float>            VeePYLpip;
    std::vector<float>            VeePZLpip;
    std::vector<float>            VeeELpip;
    std::vector<float>            VeeMLpip;
    std::vector<size_t>           VeeTAssn_VeeIDNumber;    // Being the Vee which this Assn belongs to
    std::vector<size_t>           VeeTAssn_TrackIDNumber;
    std::vector<int>              VeeTAssn_TrackEnd;

    // raw calo digits data
    unsigned int                  DiginHits;
    std::vector<float>            DigiHitX;
    std::vector<float>            DigiHitY;
    std::vector<float>            DigiHitZ;
    std::vector<float>            DigiHitTime;
    std::vector<unsigned int>     DigiHitADC;              // Uint is unsigned 32 bit integer
    std::vector<size_t>           DigiHitCellID;

    //Muon system raw hits
    unsigned int                  DiginHits_MuID;
    std::vector<float>            DigiHitX_MuID;
    std::vector<float>            DigiHitY_MuID;
    std::vector<float>            DigiHitZ_MuID;
    std::vector<float>            DigiHitTime_MuID;
    std::vector<unsigned int>     DigiHitADC_MuID;         // Uint is unsigned 32 bit integer
    std::vector<size_t>           DigiHitCellID_MuID;

    // reco calo hit data
    unsigned int                  ReconHits;
    std::vector<size_t>           ReconHitIDNumber;
    std::vector<float>            RecoHitX;
    std::vector<float>            RecoHitY;
    std::vector<float>            RecoHitZ;
    std::vector<float>            RecoHitTime;
    std::vector<float>            RecoHitEnergy;
    std::vector<size_t>           RecoHitCellID;
    float                         RecoEnergySum;

    //Muon system reco hits
    unsigned int                  ReconHits_MuID;
    std::vector<size_t>           ReconHitIDNumber_MuID;
    std::vector<float>            RecoHitX_MuID;
    std::vector<float>            RecoHitY_MuID;
    std::vector<float>            RecoHitZ_MuID;
    std::vector<float>            RecoHitTime_MuID;
    std::vector<float>            RecoHitEnergy_MuID;
    std::vector<size_t>           RecoHitCellID_MuID;
    float                         RecoEnergySum_MuID;

    // calo cluster data
    unsigned int                  nCluster;
    std::vector<size_t>           ClusterIDNumber;
    std::vector<unsigned int>     ClusterNhits;
    std::vector<float>            ClusterEnergy;
    std::vector<float>            ClusterTime;
    std::vector<float>            ClusterTimeDiffFirstLast;
    std::vector<float>            ClusterX;
    std::vector<float>            ClusterY;
    std::vector<float>            ClusterZ;
    std::vector<float>            ClusterTheta;
    std::vector<float>            ClusterPhi;
    std::vector<float>            ClusterPID;
    // std::vector<float> fClusterShape;
    std::vector<float>            ClusterMainAxisX;
    std::vector<float>            ClusterMainAxisY;
    std::vector<float>            ClusterMainAxisZ;
    std::vector<int>              ClusterMCindex;          // Branch index (NOT the GEANT track ID) of MCPartice
    std::vector<float>            ClusterMCfrac;           // that best matches & fraction of ionization therefrom

    // ECAL cluster to track association info
    std::vector<size_t>           CALAssn_ClusIDNumber;   // Being the cluster which this Assn belongs to
    std::vector<size_t>           CALAssn_TrackIDNumber;  // The rec::TrackEnd (see Track.h) that extrapolated to cluster
    std::vector<int>               CALAssn_TrackEnd;
  };
}

#endif
