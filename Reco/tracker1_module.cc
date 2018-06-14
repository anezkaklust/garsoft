////////////////////////////////////////////////////////////////////////
// Class:       tracker1
// Plugin Type: producer (art v2_11_02)
// File:        tracker1_module.cc
//
// Generated at Thu Jun 14 15:47:04 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Sort hits in X and add nearby hits.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <memory>

#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"

namespace gar {
  namespace rec {

    class tracker1 : public art::EDProducer {
    public:
      explicit tracker1(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tracker1(tracker1 const &) = delete;
      tracker1(tracker1 &&) = delete;
      tracker1 & operator = (tracker1 const &) = delete;
      tracker1 & operator = (tracker1 &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      float fHitResolYZ;     ///< resolution in cm of a hit in YZ (pad size)
      float fHitResolX;      ///< resolution in cm of a hit in X (drift direction)
      float fSigmaRoad;      ///< how many sigma away from a track a hit can be and still add it
      float fXGapToEndTrack; ///< how big a gap must be before we end a track and start a new one
      int   fMinNumHits;     ///< minimum number of hits to define a track
      std::string fHitLabel; ///< label of module creating hits
    };


    tracker1::tracker1(fhicl::ParameterSet const & p)
    {

      produces< std::vector<rec::Track> >();
      produces< ::art::Assns<rec::Hit, rec::Track> >();

      fHitResolYZ     = p.get<float>("HitResolYZ",0.7);
      fHitResolX      = p.get<float>("HitResolX",0.5);  // this is probably much better
      fSigmaRoad      = p.get<float>("SigmaRoad",3.0);
      fXGapToEndTrack = p.get<float>("XGapToEndTrack",3.0) ;
      fMinNumHits     = p.get<float>("MinNumHits",10);
      fHitLabel       = p.get<std::string>("HitLabel","hit");

    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< ::art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

      auto hitHandle = evt.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      // make an array of hit indices sorted by hit X position
      std::vector<float> hitx;
      for (size_t i=0; i<hits.size(); ++i)
	{
	  hitx.push_back(hit.Position()[0]);
	}
      std::vector<int> hsi(hitx.size());
      TMath::Sort(hitx.size(),hitx,hsi);

      // to do -- start with the first hit in the sort (lowest X), and accumulate a track using
      // local linear extrapolations.  Look for nearby hits in Y,Z at small delta X.  If none are found

      e.put(std::move(trkCol));
      e.put(std::move(hitTrkAssns));
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar
