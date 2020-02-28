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
			float numElectrons; ///< number of ionization electrons for this trackID

			HitIDE() {}

			HitIDE(int id, float ef, float e)
				: trackID(id), energyFrac(ef), numElectrons(e)
			{}
    	};
		struct CalIDE {
			int   trackID;      ///< Geant4 supplied trackID
			float energyFrac;   ///< fraction of gar::rec::CaloHit energy from particle with this trackID

			CalIDE() {}

			CalIDE(int id, float ef)
				: trackID(id), energyFrac(ef)
			{}
    	};



		class BackTrackerCore {
		public:
      
			explicit BackTrackerCore(fhicl::ParameterSet const& pset);
			~BackTrackerCore();

			// Set the EveIdCalculator for the owned ParticleList
			void AdoptEveIdCalculator(sim::EveIdCalculator* ec) {
				fParticleList.AdoptEveIdCalculator(ec);
			}



			// Get a reference to the ParticleLists
			sim::ParticleList const& ParticleList() const {
				return fParticleList;
			}

			// Returns a bare pointer to the MCParticle corresponding to a given TrackID
			const simb::MCParticle* TrackIDToParticle(int const& id) const;

 			const simb::MCParticle* FindMother(simb::MCParticle* const p) const;

			const simb::MCParticle* FindEve(simb::MCParticle* const p) const;

			bool IsDescendedFrom(simb::MCParticle* const p, simb::MCParticle* const c) const;
      
			// Returns an ::art::Ptr<> to simb::MCTruth
			art::Ptr<simb::MCTruth> const& ParticleToMCTruth(const simb::MCParticle* p) const;



			// These methods will return HitIDE structures containing the Geant4 track IDs of
			// the particles contributing ionization electrons to the input hit
			std::vector<HitIDE> HitToHitIDEs(art::Ptr<gar::rec::Hit> const& hit) const;
			std::vector<HitIDE> HitToHitIDEs(         gar::rec::Hit  const& hit) const;

			std::vector<::art::Ptr<gar::rec::Hit>> const ParticleToHits(simb::MCParticle* const p) const;

			// method to return the fraction of hits in a collection that come from the specified Geant4 track ids
			double HitCollectionPurity(simb::MCParticle* const p,
									   std::vector<art::Ptr<gar::rec::Hit> > const& hits,
									   bool weightByCharge=false) const;

			// method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
			// represented in a collection of hits
			double HitCollectionEfficiency(simb::MCParticle* const p,
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

			const detinfo::DetectorClocks*                        fClocks;                ///< Detector clock information
			sim::ParticleList                                     fParticleList;          ///< ParticleList to map track ID to sim::Particle
			std::vector< std::vector<const sdp::EnergyDeposit*> > fChannelToEDepCol;      ///< convenience collections of EnergyDeposits for each channel
			std::vector< ::art::Ptr<simb::MCTruth> >              fMCTruthList;           ///< all the MCTruths for the event
			std::map<int, int>                                    fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry
			std::string                                           fG4ModuleLabel;         ///< label for geant4 module
			std::string                                           fRawDataLabel;          ///< label for geant4 module
			double                                                fMinHitEnergyFraction;  ///< minimum fraction of energy a track id has to
                                                                              		      ///< contribute to a hit to be counted in
                                                                              		      ///< purity and efficiency calculations
                                                                              		      ///< based on hit collections
			double                                                fInverseVelocity;       ///< inverse drift velocity
			double                                                fLongDiffConst;         ///< longitudinal diffusion constant

			const geo::GeometryCore*                              fGeo;                   ///< pointer to the geometry
			bool                                                  fDisableRebuild;        ///< for switching off backtracker's rebuild of the MCParticle tables
		};
	} // namespace
}


#endif /* GAR_CHEAT_BackTrackerCore_h */
