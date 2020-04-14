//
//  BackTrackerCore.h
//
//  Created by Brian Rebel on 3/13/17.
//  Modified by Leo Bellantoni, Feb-Mar 2020
//
//
/// Code to link reconstructed objects back to the MC truth information



#ifndef GAR_CHEAT_BackTrackerCore_h
#define GAR_CHEAT_BackTrackerCore_h



#include "canvas/Persistency/Common/Ptr.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "DetectorInfo/DetectorClocksService.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "Geometry/Geometry.h"

#include <unordered_map>



namespace sim {
    class ParticleList;
    class EveIdCalculator;
}



namespace gar{
    namespace cheat{



        // Two helper classes structs for doing backtracking analyses
        struct HitIDE {
            int   trackID;      ///< Geant4 trackID, as modified by ParticleListAction::PreTrackingAction
            float energyFrac;   ///< fraction of gar::rec::Hit energy from the particle with this trackID
            float energyTot;    ///< total energy for this trackID.  In units of probably-GeV.

            HitIDE() {}

            HitIDE(int id, float ef, float e)
                : trackID(id), energyFrac(ef), energyTot(e) {}
        };
        struct CalIDE {
            int   trackID;      ///< Geant4 trackID, as modified by ParticleListAction::PreTrackingAction
            float energyFrac;   ///< fraction of gar::rec::CaloHit energy from particle with this trackID
            float energyTot;    ///< total energy for this trackID.  In units of probably-GeV.

            CalIDE() {}

            CalIDE(int id, float ef, float e)
                : trackID(id), energyFrac(ef), energyTot(e) {}
        };



        class BackTrackerCore {
        public:
      
            explicit BackTrackerCore(fhicl::ParameterSet const& pset);
            ~BackTrackerCore();

            // Set the EveIdCalculator for the owned ParticleList
            void AdoptEveIdCalculator(sim::EveIdCalculator* ec) {
                fParticleList.AdoptEveIdCalculator(ec);
            }



            // Following should work as long as there are MCParticles and MCTruth
            // data products in the event, and suitable Assns between them

            // ParticleList has no copy constructor.  It does have iterators etc.
            sim::ParticleList* GetParticleList() {
                return &fParticleList;
            }

            // Returns a bare pointer to the MCParticle corresponding to a given TrackID
            // Negative trackID values, corresponding to EM shower particles, return
            // the parent particle.
            simb::MCParticle* const TrackIDToParticle(int const& id) const;

            simb::MCParticle* const FindMother(simb::MCParticle* const p) const {
                // Always walk up the ParticleList in case someday somebody changes it
                return TrackIDToParticle(p->Mother());
            }

            simb::MCParticle* const FindEve(simb::MCParticle* const p) const;

            // FindTPCEve is not a const method because it saves result of previous 
            // calls in `this' for efficiency.  2nd signature mostly for internal use.
            simb::MCParticle* const FindTPCEve(simb::MCParticle* const p) const;
            simb::MCParticle* const FindTPCEve(int const trackId) const;

            bool IsDescendedFrom(simb::MCParticle* const forebear,
                                 simb::MCParticle* const afterbear) const;
      
            // Returns an ::art::Ptr<> to simb::MCTruth
            art::Ptr<simb::MCTruth> const ParticleToMCTruth(simb::MCParticle* const p) const;



            // Following should work as long as there are also RawDigit and EnergyDeposit
            // data products in the event, and suitable Assns between them

            // These methods will return HitIDE structures containing the Geant4 track IDs of
            // the particles contributing ionization electrons to the input hit
            std::vector<HitIDE> HitToHitIDEs(art::Ptr<rec::Hit> const& hit) const;
            std::vector<HitIDE> HitToHitIDEs(         rec::Hit  const& hit) const;

            std::vector<art::Ptr<rec::Hit>>
            ParticleToHits(simb::MCParticle* const p,
                           std::vector<art::Ptr<rec::Hit>> const& allhits,
                           bool checkNeutrals=false) const;

            // Find the fraction of hits in a collection that come from the specified MCParticle
            // with statistical uncertainty.  Binomial uncertainty in weighted events is from
            // the usual primitive algorithm, which is not obviously right to me.
            std::pair<double,double>
            HitPurity(simb::MCParticle* const p,
                      std::vector<art::Ptr<rec::Hit>> const& hits,
                      bool weightByCharge=false) const;

            // method to return the fraction of all hits in an event from a specific set of 
            // Geant4 track IDs that arerepresented in a collection of hits
            std::pair<double,double>
            HitEfficiency(simb::MCParticle* const p,
                          std::vector<art::Ptr<rec::Hit> > const& hits,
                          std::vector<art::Ptr<rec::Hit> > const& allhits,
                          bool weightByCharge=false) const;



            // Following should work as long as there are CaloRawDigit and CaloDeposit
            // data products in the event, and suitable Assns between them

            // These methods will return HitIDE structures containing the Geant4 track IDs of
            // the TPC-originating particles contributing ionization electrons to the input hit.
            std::vector<CalIDE> CaloHitToCalIDEs(art::Ptr<rec::CaloHit> const& hit) const;
            std::vector<CalIDE> CaloHitToCalIDEs(         rec::CaloHit  const& hit) const;

            std::vector<art::Ptr<rec::CaloHit>>
            ParticleToCaloHits(simb::MCParticle* const p,
                               std::vector<art::Ptr<rec::CaloHit>> const& allhits) const;

            // Find the fraction of CaloHits in a collection that come from the specified 
            // MCParticle with statistical uncertainty.  Binomial uncertainty in weighted 
            // events is from the usual primitive algorithm, which is not obviously right to me.
            std::pair<double,double>
            CaloHitPurity(simb::MCParticle* const p,
                          std::vector<art::Ptr<rec::CaloHit>> const& hits,
                          bool weightByCharge=false) const;

            // method to return the fraction of all hits in an event from a specific set of 
            // Geant4 track IDs that arerepresented in a collection of hits
            std::pair<double,double>
            CaloHitEfficiency(simb::MCParticle* const p,
                              std::vector<art::Ptr<rec::CaloHit>> const& hits,
                              std::vector<art::Ptr<rec::CaloHit>> const& allhits,
                              bool weightByCharge=false) const;



            // Following should work as long as there are reconstructed Track and Hit
            // data products in the event, and suitable Assns between them
            std::vector<art::Ptr<rec::Hit>> const TrackToHits(rec::Track* const t);

            std::vector<std::pair<simb::MCParticle*,float>> TrackToMCParticles(rec::Track* const t);

            std::vector<art::Ptr<rec::Track>>
            MCParticleToTracks(simb::MCParticle* const p,
                              std::vector<art::Ptr<rec::Track>> const& tracks);



            // Following should work as long as there are reconstructed Cluster and CaloHit
            // data products in the event, and suitable Assns between them
            std::vector<art::Ptr<rec::CaloHit>> const ClusterToCaloHits(rec::Cluster* const c);

            std::vector<std::pair<simb::MCParticle*,float>> ClusterToMCParticles(rec::Cluster* const c);

            std::vector<art::Ptr<rec::Cluster>>
            MCParticleToClusters(simb::MCParticle* const p,
                              std::vector<art::Ptr<rec::Cluster>> const& clusters);


      
        private:

            std::vector<HitIDE>
            ChannelToHitIDEs(raw::Channel_t const& channel,
                             double const start, double const stop) const;

            std::vector<CalIDE>
            CellIDToCalIDEs(raw::CellID_t const& cellID, float const time) const;



        protected:

            bool fHasMC, fHasHits, fHasCalHits, fHasTracks, fHasClusters;

            const detinfo::DetectorClocks*                      fClocks;                ///< Detector clock information
            const geo::GeometryCore*                            fGeo;                   ///< pointer to the geometry

            bool                                                fDisableRebuild;        ///< for switching off backtracker's rebuild of the MCParticle tables
            std::string                                         fG4ModuleLabel;         ///< label for geant4 module
            std::string                                         fRawTPCDataLabel;       ///< label for TPC readout module
            std::string                                         fRawCaloDataLabel;      ///< label for ECAL readout module
            double                                              fECALtimeResolution;    ///< time resolution for hits in ECAL, nsec.
            double                                              fMinHitEnergyFraction;  ///< min frac of ionization a track has to count in a TPC hit
            double                                              fMinCaloHitEnergyFrac;  ///< min frac of ionization a track has to count in a CaloHit
            std::string                                         fTrackLabel;            ///< label for final track producing module
            std::string                                         fTPCClusterLabel;       ///< label for TPCCluster producing module
            double                                              fTrackFracMCP;          ///< min frac of ionization in a track for matching to an MCParticle
            std::string                                         fClusterLabel;          ///< label for ECAL cluster producing module
            double                                              fClusterFracMCP;        ///< min frac of ionization in a cluster for matching to an MCParticle

            sim::ParticleList                                   fParticleList;          ///< Maps MCParticle::TrackId() to same MCParticle
            std::vector<art::Ptr<simb::MCTruth>>                fMCTruthList;           ///< all the MCTruths for the event
            std::unordered_map<int, int>                        fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry.  Track Ids from MCParticle
                                                                                        ///< table in event store (no track ID <0)
            std::unordered_map<int, int>*                       fECALTrackToTPCTrack;   ///< results of previous FindTPCEve calls

            double                                              fInverseVelocity;       ///< inverse drift velocity
            double                                              fLongDiffConst;         ///< longitudinal diffusion constant
            // vector gives a fast lookup and is only 677864*24 = 16Mbytes.
            std::vector<std::vector<const sdp::EnergyDeposit*>> fChannelToEDepCol;      ///< convenience collections of EnergyDeposits for each channel

            // ECAL max channel id is something like 2^56.  Try an unordered map, not a vector.
            std::unordered_map<raw::CellID_t,std::vector<const sdp::CaloDeposit*>>
                                                                fCellIDToEDepCol;       ///< convenience collections of EnergyDeposit for each cell

            // Mapping final reco products to their constituents
            std::unordered_map< rec::IDNumber, std::vector<art::Ptr<rec::Hit>> >
                                                                fTrackIDToHits;         ///< Reco track ID to track's hits
            std::unordered_map< rec::IDNumber, std::vector<art::Ptr<rec::CaloHit>> >
                                                                fClusterIDToCaloHits;   ///< Reco ECAL cluster ID to CaloHits

        };
    } // namespace
}


#endif /* GAR_CHEAT_BackTrackerCore_h */
