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
#include "ReconstructionDataProducts/Hit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"




namespace sim {
    class ParticleList;
    class EveIdCalculator;
}



namespace gar{
    namespace cheat{



        // Two helper classes structs for doing backtracking analyses
        struct HitIDE {
            int   trackID;      ///< Geant4 supplied trackID
            float energyFrac;   ///< fraction of gar::rec::Hit energy from the particle with this trackID
            float energyTot;    ///< total energy for this trackID.  In units of probably-GeV.

            HitIDE() {}

            HitIDE(int id, float ef, float e)
                : trackID(id), energyFrac(ef), energyTot(e) {}
        };
        struct CalIDE {
            int   trackID;      ///< Geant4 supplied trackID
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

            // Get pointer to the ParticleList - BTW there's no copy constructor
            // in that class.  Ooooo kay!
            sim::ParticleList const& GetParticleList() const {
                return fParticleList;
            }

            // Returns a bare pointer to the MCParticle corresponding to a given TrackID
            // Negative trackID values work just fine.
            simb::MCParticle* const TrackIDToParticle(int const& id) const;

            simb::MCParticle* const FindMother(simb::MCParticle* const p) const {
                return TrackIDToParticle(p->Mother());
            }

            simb::MCParticle* const FindEve(simb::MCParticle* const p) const;

            bool IsDescendedFrom(simb::MCParticle* const forebear,
                                 simb::MCParticle* const afterbear) const;
      
            // Returns an ::art::Ptr<> to simb::MCTruth
            art::Ptr<simb::MCTruth> const ParticleToMCTruth(simb::MCParticle* const p) const;



            // Following should work as long as there are also RawDigit and EnergyDeposit
            // data products in the event, and suitable Assns between them

            // These methods will return HitIDE structures containing the Geant4 track IDs of
            // the particles contributing ionization electrons to the input hit
            std::vector<HitIDE> HitToHitIDEs(art::Ptr<gar::rec::Hit> const& hit) const;
            std::vector<HitIDE> HitToHitIDEs(         gar::rec::Hit  const& hit) const;

            std::vector<art::Ptr<gar::rec::Hit>>
            ParticleToHits(simb::MCParticle* const p,
                           std::vector<art::Ptr<gar::rec::Hit>> const& allhits,
                           bool checkNeutrals=false) const;

            // Fins the fraction of hits in a collection that come from the specified MCParticle
            // with statistical uncertainty.  Binomial uncertainty in weighted events is from
            // the usual primitive algorithm, which is not obviously right to me.
            std::pair<double,double>
            HitCollectionPurity(simb::MCParticle* const p,
                                std::vector<art::Ptr<gar::rec::Hit>> const& hits,
                                bool weightByCharge=false) const;

            // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
            // represented in a collection of hits
            std::pair<double,double>
            HitCollectionEfficiency(simb::MCParticle* const p,
                                    std::vector<art::Ptr<gar::rec::Hit> > const& hits,
                                    std::vector<art::Ptr<gar::rec::Hit> > const& allhits,
                                    bool weightByCharge=false) const;





      
      // this method determines the appropriate EnergyDeposits for a given TDC
      // on a channel
      std::map<size_t, const sdp::EnergyDeposit*> ChannelTDCToEnergyDeposit(raw::Channel_t channel) const;
      
      
      // method to return a subset of allhits that are matched to a list of TrackIDs
      std::vector<std::vector<::art::Ptr<gar::rec::Hit>>> const TrackIDsToHits(std::vector<::art::Ptr<gar::rec::Hit>> const& allhits,
                                                                               std::vector<int>                       const& tkIDs) const;
      
      // method to return the EveIDs of particles contributing ionization
      // electrons to the identified hit
      std::vector<HitIDE> HitToEveID(::art::Ptr<gar::rec::Hit> const& hit) const;
      
      // method to return the XYZ position of the weighted average energy deposition for a given hit
      std::vector<float>  HitToXYZ(::art::Ptr<gar::rec::Hit> const& hit) const;
      
      
      // method to return all TrackIDs corresponding to the given list of hits
      std::set<int> GetSetOfTrackIDs() const;



        private:

            std::vector<HitIDE>
            ChannelToHitIDEs(gar::raw::Channel_t const& channel,
                             double const start, double const stop) const;



        protected:

            bool fHasMC, fHasHits, fHasCaloHits, fHasTracks, fHasClusters;

            const detinfo::DetectorClocks*                      fClocks;                ///< Detector clock information
            const geo::GeometryCore*                            fGeo;                   ///< pointer to the geometry

            std::string                                         fG4ModuleLabel;         ///< label for geant4 module
            std::string                                         fRawDataLabel;          ///< label for geant4 module
            double                                              fMinHitEnergyFraction;  ///< minimum fraction of energy a track id has to
             bool                                               fDisableRebuild;        ///< for switching off backtracker's rebuild of the MCParticle tables
                                                                                          ///< contribute to a hit to be counted in
                                                                                            ///< purity and efficiency calculations
                                                                                            ///< based on hit collections

            sim::ParticleList                                   fParticleList;          ///< ParticleList to map track ID to sim::Particle
            std::vector<art::Ptr<simb::MCTruth>>                fMCTruthList;           ///< all the MCTruths for the event
            std::map<int, int>                                  fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry

            std::vector<std::vector<const sdp::EnergyDeposit*>> fChannelToEDepCol;      ///< convenience collections of EnergyDeposits for each channel

            double                                              fInverseVelocity;       ///< inverse drift velocity
            double                                              fLongDiffConst;         ///< longitudinal diffusion constant

        };
    } // namespace
}


#endif /* GAR_CHEAT_BackTrackerCore_h */
