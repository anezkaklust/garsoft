////////////////////////////////////////////////////////////////////////
/// \file  BackTracker.h
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GAR_CHEAT_BACKTRACKER_H
#define GAR_CHEAT_BACKTRACKER_H

#include <vector>
#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "SimulationDataProducts/SimChannel.h"
#include "ReconstructionDataProducts/Hit.h"
#include "RawDataProducts/RawDigit.h"

///code to link reconstructed objects back to the MC truth information

namespace sim {
  class ParticleList;
  class EveIdCalculator;
}

namespace gar{
  namespace cheat{
    
    // make convenience structs for doing the backtracking
    struct HitIDE {
      
      int   trackID;      ///< Geant4 supplied trackID
      float energyFrac;   ///< fraction of hit energy from the particle with this trackID
      float numElectrons; ///< number of ionization electrons for this trackID
      HitIDE()
      {}
      
      HitIDE(int   id,
             float ef,
             float e)
      : trackID(id)
      , energyFrac(ef)
      , numElectrons(e)
      {}
      
    };
    
    class BackTracker {
      
    public:
      
      BackTracker(fhicl::ParameterSet   const& pset,
                  ::art::ActivityRegistry      & reg);
      ~BackTracker();
      
      void reconfigure(fhicl::ParameterSet const& pset);
      
      // The Rebuild function rebuilds the various maps we need to answer backtracking queries.
      // It is called automatically before each event is processed. For jobs involving
      // Monte Carlo generation, this is too soon. So, we'll call rebuild after those data
      // products are put into the event in LArG4.  This is the least bad way of ensuring the
      // BackTracker works in jobs that combine MC production and reconstruction analysis based
      // on MC truth.  Don't copy this design pattern without talking to brebel@fnal.gov first
      void Rebuild(::art::Event const& evt);
      
      // Get a reference to the ParticleList
      sim::ParticleList const&             ParticleList()  const { return fParticleList; }
      
      // Set the EveIdCalculator for the owned ParticleList
      void                                 AdoptEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }
      
      // Return a pointer to the simb::MCParticle object corresponding to
      // the given TrackID
      const simb::MCParticle*              TrackIDToParticle(int const& id)       const;
      const simb::MCParticle*              TrackIDToMotherParticle(int const& id) const;
      
      // Get ::art::Ptr<> to simb::MCTruth and related information
      const ::art::Ptr<simb::MCTruth>&       TrackIDToMCTruth(int const& id)                        const;
      const ::art::Ptr<simb::MCTruth>&       ParticleToMCTruth(const simb::MCParticle* p)           const;
      std::vector<const simb::MCParticle*> MCTruthToParticles(::art::Ptr<simb::MCTruth> const& mct) const;
      const std::vector< ::art::Ptr<simb::MCTruth> >& MCTruthVector() const { return fMCTruthList; }
      
      // this method will return the Geant4 track IDs of
      // the particles contributing ionization electrons to the identified hit
      std::vector<HitIDE> HitToTrackID(::art::Ptr<gar::rec::Hit> const& hit);
      
      // method to return a subset of allhits that are matched to a list of TrackIDs
      const std::vector<std::vector<::art::Ptr<gar::rec::Hit>>> TrackIDsToHits(std::vector<::art::Ptr<gar::rec::Hit>> const& allhits,
                                                                               std::vector<int>                       const& tkIDs);
      
      // method to return the EveIDs of particles contributing ionization
      // electrons to the identified hit
      std::vector<HitIDE> HitToEveID(::art::Ptr<gar::rec::Hit> const& hit);
      
      // method to return the XYZ position of the weighted average energy deposition for a given hit
      //std::vector<double>                   HitToXYZ(::art::Ptr<gar::rec::Hit> const& hit);
      
      // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
      double HitCollectionPurity(std::set<int>                                   trackIDs,
                                 std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                 bool                                            weightByCharge=false);
      
      // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
      // represented in a collection of hits
      double HitCollectionEfficiency(std::set<int>                                   trackIDs,
                                     std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                     std::vector< ::art::Ptr<gar::rec::Hit> > const& allhits,
                                     bool                                            weightByCharge=false);
      
      // method to return all EveIDs corresponding to the current sim::ParticleList
      std::set<int> GetSetOfEveIDs();
      
      // method to return all TrackIDs corresponding to the current sim::ParticleList
      std::set<int> GetSetOfTrackIDs();
      
      // method to return all EveIDs corresponding to the given list of hits
      std::set<int> GetSetOfEveIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits);
      
      // method to return all TrackIDs corresponding to the given list of hits
      std::set<int> GetSetOfTrackIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits);
      
    private:
      
      void                   ChannelToTrackID(std::vector<HitIDE>      & hitIDEs,
                                              gar::raw::Channel_t const& channel,
                                              unsigned short      const  start_tdc,
                                              unsigned short      const  end_tdc);
      
      sim::ParticleList                        fParticleList;          ///< ParticleList to map track ID to sim::Particle
      std::vector< std::vector<sdp::IDE> >     fChannelToIDEs;         ///< convenience collections of IDEs for each channel
      std::vector< ::art::Ptr<simb::MCTruth> > fMCTruthList;           ///< all the MCTruths for the event
      std::map<int, int>                       fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry
      std::string                              fG4ModuleLabel;         ///< label for geant4 module
      std::string                              fIonizationModuleLabel; ///< label for geant4 module
      double                                   fMinHitEnergyFraction;  ///< minimum fraction of energy a track id has to
                                                                       ///< contribute to a hit to be counted in
                                                                       ///< purity and efficiency calculations
                                                                       ///< based on hit collections
    };
  } // namespace
}

DECLARE_ART_SERVICE(gar::cheat::BackTracker, LEGACY)

#endif // CHEAT_BACKTRACK_H
