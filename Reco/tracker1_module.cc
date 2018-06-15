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
      fSigmaRoad      = p.get<float>("SigmaRoad",5.0);
      fMinNumHits     = p.get<float>("MinNumHits",10);
      fHitLabel       = p.get<std::string>("HitLabel","hit");

    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< ::art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      // make an array of hit indices sorted by hit X position
      std::vector<float> hitx;
      for (size_t i=0; i<hits.size(); ++i)
	{
	  hitx.push_back(hits[i].Position()[0]);
	}
      std::vector<int> hsi(hitx.size());
      TMath::Sort((int) hitx.size(),hitx.data(),hsi.data());

      float roadsq = fSigmaRoad*fSigmaRoad;

      // record which hits we have assigned to which tracks
      std::vector< std::vector<int> > hitlist;

      for (size_t i=0; i<hits.size(); ++i)
	{
	  const float *hpos = hits[hsi[i]].Position();
	  float bestsignifs = -1;
	  int ibest = -1;
	  for (size_t itcand = 0; itcand < hitlist.size(); ++itcand)
	    {
	      const float *cpos = hits[hsi[hitlist[itcand].back()]].Position();
	      float signifs = TMath::Sq( (hpos[0]-cpos[0])/fHitResolX ) + 
		TMath::Sq( (hpos[1]-cpos[1])/fHitResolYZ ) +
		TMath::Sq( (hpos[3]-cpos[3])/fHitResolYZ );
	      if (bestsignifs < 0 || signifs < bestsignifs)
		{
		  bestsignifs = signifs;
		  ibest = itcand;
		}
	    }
	  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
	    {
	      std::vector<int> hloc;
	      hloc.push_back(i);
	      hitlist.push_back(hloc);
	    }
	  else  // add the hit to the existing best track
	    {
	      hitlist[ibest].push_back(i);
	    }
	}

      // now that we have the hits assigned, fit the tracks and adjust the hit lists.

      e.put(std::move(trkCol));
      e.put(std::move(hitTrkAssns));
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar
