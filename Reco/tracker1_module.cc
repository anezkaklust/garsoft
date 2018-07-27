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
#include "art/Persistency/Common/PtrMaker.h"

#include <memory>

#include "TVectorF.h"
#include "TMatrix.h"

#include "Geant4/G4ThreeVector.hh"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/MagneticField/MagneticField.h"

// GArSoft Includes
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"

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

      float fHitResolYZ;           ///< resolution in cm of a hit in YZ (pad size)
      float fHitResolX;            ///< resolution in cm of a hit in X (drift direction)
      float fSigmaRoad;            ///< how many sigma away from a track a hit can be and still add it during patrec
      float fXGapToEndTrack;       ///< how big a gap must be before we end a track and start a new one
      unsigned int fMinNumHits;    ///< minimum number of hits to define a track
      std::string fHitLabel;       ///< label of module creating hits
      int fPrintLevel;             ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      int fTrackPass;              ///< which pass of the tracking to save as the tracks in the event
      int fDumpTracks;             ///< 0: do not print out tracks, 1: print out tracks

      float fHitResolYZinFit;      ///< Hit resolution parameter to use in fit
      float fRoadYZinFit;          ///< cut in cm for dropping hits from tracks in fit

      int initial_trackpar_estimate(art::ValidHandle<std::vector<Hit> > &hitHandle, 
				    std::vector<std::vector<int> >      &hitlist,
				    std::vector<int>                    &hsi,
				    int itrack, 
				    bool isForwards,
				    float &curvature_init,
				    float &slope_init,
				    float &phi_init,
				    float &xpos,
				    float &ypos,
				    float &zpos);

      int KalmanFit( art::ValidHandle<std::vector<Hit> > &hitHandle, 
		     std::vector<std::vector<int> > &hitlist, 
		     std::vector<int> &hsi, 
		     int itrack, 
		     bool isForwards, 
		     std::vector<float> &trackparatend,
		     float &chisquared,
		     float &length,
		     float *covmat,    // 5x5 covariance matrix
		     std::vector<int> &unused_hits);

      int FitHelix(art::ValidHandle<std::vector<Hit> > &hitHandle, 
		   std::vector<std::vector<int> > &hitlist, 
		   std::vector<int> &hsi, 
		   int itrack, 
		   bool isForwards, 
		   std::vector<int> &unused_hits,
		   TrackPar &trackpar
		   );

      size_t ifob(size_t ihit, size_t nhits, bool isForwards);
 
      float capprox(float x1,float y1,
		    float x2,float y2,
		    float x3,float y3);  ///< initial guess of curvature calculator -- from ALICE

      float capprox2(float y0, float z0, float y1, float z1, float y2, float z2);  // redo -- returns abs value of curvature

    };


    tracker1::tracker1(fhicl::ParameterSet const & p)
    {

      produces< std::vector<rec::Track> >();
      produces< art::Assns<rec::Hit, rec::Track> >();

      fHitResolYZ      = p.get<float>("HitResolYZ",1.0); // TODO -- think about what this value is
      fHitResolX       = p.get<float>("HitResolX",0.5);  // this is probably much better
      fSigmaRoad       = p.get<float>("SigmaRoad",5.0);
      fMinNumHits      = p.get<unsigned int>("MinNumHits",20);
      fHitLabel        = p.get<std::string>("HitLabel","hit");
      fPrintLevel      = p.get<int>("PrintLevel",0);
      fTrackPass       = p.get<int>("TrackPass",2);
      fDumpTracks      = p.get<int>("DumpTracks",2);
      fHitResolYZinFit = p.get<float>("HitResolYZinFit",4.0);
      fRoadYZinFit     = p.get<float>("RoadYZinFit",1.0);
    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e, *this);
      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, hitHandle.id());

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

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
      std::vector< std::vector<int> > hitlistf;
      std::vector<int> trackf(hits.size());
      std::vector< std::vector<int> > hitlistb;
      std::vector<int> trackb(hits.size());

      std::vector< std::vector<int> > hitlist;
      std::vector<int> whichtrack(hits.size(),-1);

      float resolSq = fHitResolYZ*fHitResolYZ;

      // idea -- find all of the track candidates a hit can be added to that have signfs <= roadsq, and
      // pick the one that is closest to the local linear extrapolation of the existing hits.

      // do this twice, once going forwards through the hits, once backwards, and then collect hits
      // in groups and split them when either the forwards or the backwards list says to split them.

      for (size_t ihit=0; ihit<hits.size(); ++ihit)
	{
	  const float *hpos = hits[hsi[ihit]].Position();
	  float bestsignifs = -1;
	  int ibest = -1;
	  for (size_t itcand = 0; itcand < hitlistf.size(); ++itcand)
	    {
	      const float *cpos = hits[hsi[hitlistf[itcand].back()]].Position();
	      float signifs = // TMath::Sq( (hpos[0]-cpos[0])/fHitResolX ) + 
		(TMath::Sq( (hpos[1]-cpos[1]) ) +
		 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
	      if (bestsignifs < 0 || signifs < bestsignifs)
		{
		  bestsignifs = signifs;
		  ibest = itcand;
		}
	    }
	  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
	    {
	      ibest = hitlistf.size();
	      std::vector<int> vtmp;
	      hitlistf.push_back(vtmp);
	    }
	  hitlistf[ibest].push_back(ihit);
	  trackf[ihit] = ibest;
	}

      for (int ihit=hits.size()-1; ihit >= 0; --ihit)
	{
	  const float *hpos = hits[hsi[ihit]].Position();
	  float bestsignifs = -1;
	  int ibest = -1;
	  for (size_t itcand = 0; itcand < hitlistb.size(); ++itcand)
	    {
	      const float *cpos = hits[hsi[hitlistb[itcand].back()]].Position();
	      float signifs = // TMath::Sq( (hpos[0]-cpos[0])/fHitResolX ) + 
		(TMath::Sq( (hpos[1]-cpos[1]) ) +
		 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
	      if (bestsignifs < 0 || signifs < bestsignifs)
		{
		  bestsignifs = signifs;
		  ibest = itcand;
		}
	    }
	  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
	    {
	      ibest = hitlistb.size();
	      std::vector<int> vtmp;
	      hitlistb.push_back(vtmp);
	    }
	  hitlistb[ibest].push_back(ihit);
	  trackb[ihit] = ibest;
	}

      // make a list of tracks that is the set of disjoint subsets

      for (size_t itrack=0; itrack<hitlistf.size(); ++itrack)
	{
	  int itrackl = 0;
	  for (size_t ihit=0; ihit<hitlistf[itrack].size(); ++ihit)
	    {
	      int ihif = hitlistf[itrack][ihit];
	      if (ihit == 0 || (trackb[ihif] != itrackl))
		{
		  std::vector<int> vtmp;
		  hitlist.push_back(vtmp);
		  itrackl = trackb[ihif];
		}
	      hitlist.back().push_back(ihif);
	    }
	}



      // now that we have the hits assigned, fit the tracks and adjust the hit lists.
      // fit them in both directions.  The Kalman filter gives the most precise measurement
      // of position and momentum at the end.

      // do a first pass of fitting the tracks

      std::vector<TrackPar> firstpass_tracks;
      std::vector<int> firstpass_tid;
      std::vector<TrackPar> secondpass_tracks;
      std::vector<int> secondpass_tid;
      float covmatbeg[25];
      float covmatend[25];
      size_t ntracks = hitlist.size();

      if (fDumpTracks > 0)
	{
	  int ntracktmp = 0;
	  for (size_t itrack=0;itrack<ntracks;++itrack)
	    {
	      if (hitlist[itrack].size() >= fMinNumHits) ntracktmp++;
	    }
	  std::cout << "Trkdump: " << ntracktmp << std::endl;
	}
      for (size_t itrack=0; itrack<ntracks; ++itrack)
	{
	  size_t nhits = hitlist[itrack].size();
	  if ( nhits >= fMinNumHits)
	    {
	      for (size_t ihit=0; ihit<nhits; ++ihit)
		{
		  whichtrack[hitlist[itrack][ihit]] = itrack;  // fill this here so we only get the ones passing the nhits cut
		}

	      if (fPrintLevel)
		{
		  std::cout << "Starting a new Pass1 track: " << itrack << " Number of hits: " << nhits << std::endl;
		}

	      // variables:  x is the independent variable
	      // 0: y
	      // 1: z
	      // 2: curvature
	      // 3: phi
	      // 4: slope = d(yz distance)/dx
	      // 5: x   /// added on to the end

	      // the "forward" fit is just in increasing x.  Track parameters are at the end of the fit

	      std::vector<float> tparend(6);
	      float chisqforwards = 0;
	      float lengthforwards = 0;
	      std::vector<int> unused_hits;  
	      int retcode = KalmanFit(hitHandle,hitlist,hsi,itrack,true,tparend,chisqforwards,lengthforwards,covmatend,unused_hits);
	      if (retcode != 0) continue;
              unused_hits.clear();  // todo -- use these hits for something

	      // the "backwards" fit is in decreasing x.  Track paramters are at the end of the fit, the other end of the track

	      std::vector<float> tparbeg(6);
	      float chisqbackwards = 0;
	      float lengthbackwards = 0;
	      retcode = KalmanFit(hitHandle,hitlist,hsi,itrack,false,tparbeg,chisqbackwards,lengthbackwards,covmatbeg,unused_hits);
	      if (retcode != 0) continue;
              unused_hits.clear();  // todo -- use these hits for something

	      firstpass_tracks.emplace_back(lengthforwards,
					    lengthbackwards,
					    nhits,
					    tparbeg[5],
					    tparbeg.data(),
					    covmatbeg,
					    chisqforwards,
					    tparend[5],
					    tparend.data(),
					    covmatend,
					    chisqbackwards,
					    0);  // zero timestamp for now
	      firstpass_tid.push_back(itrack);

	      if (fDumpTracks > 0)
		{
		  std::cout << "Trkdump: " << itrack << std::endl;
		  std::cout << "Trkdump: " << tparbeg[5] << std::endl;
		  for (int i=0; i<5;++i) std::cout << "Trkdump: " << tparbeg[i] << std::endl;
		  std::cout << "Trkdump: " << tparend[5] << std::endl;
		  for (int i=0; i<5;++i) std::cout << "Trkdump: " << tparend[i] << std::endl;
		  std::cout << "Trkdump: " << nhits << std::endl;
		  for (size_t ihit=0;ihit<nhits;++ihit)
		    {
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[0] << std::endl;
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[1] << std::endl;
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[2] << std::endl;
		    }
		}
	    }
	}

      // Rearrange the hit lists -- ask ourselves which track each hit is best assigned to.  Make a new hit list, hitlist2
      // start only with hits that were assigned to tracks in pass1 (i.e. don't try to add new hits that made short, faraway tracks)
      // todo -- change this last requirement to an absolute distance cut.  Debug the track parameter extrapolation first.

      std::vector< std::vector<int> > hitlist2(firstpass_tracks.size());

      for (size_t ihit=0; ihit< hits.size(); ++ihit)
	{
	  if (whichtrack[ihit] < 0) continue;
	  const float *hpos = hits[hsi[ihit]].Position();
          float mindist = 0;
	  size_t ibest = 0;
	  for (size_t itrack=0; itrack<firstpass_tracks.size(); ++itrack)
	    {
	      float dist = firstpass_tracks[itrack].DistXYZ(hpos);
	      if (itrack == 0 || dist < mindist) 
		{
		  mindist = dist;
		  ibest = itrack;
		}
	    }
          hitlist2[ibest].push_back(ihit);	  
	}

      size_t ntracks2 = hitlist2.size();
      for (size_t itrack=0; itrack<ntracks2; ++itrack)
	{
	  size_t nhits = hitlist2[itrack].size();
	  if ( nhits >= fMinNumHits)
	    {
	      if (fPrintLevel)
		{
		  std::cout << "Starting a new Pass2 track: " << itrack << " Number of hits: " << nhits << std::endl;
		}

	      std::vector<float> tparend(6);
	      float chisqforwards = 0;
	      float lengthforwards = 0;
	      std::vector<int> unused_hits;
	      int retcode = KalmanFit(hitHandle,hitlist2,hsi,itrack,true,tparend,chisqforwards,lengthforwards,covmatend,unused_hits);
	      if (retcode != 0) continue;
              unused_hits.clear();  // todo -- use these hits for something

	      std::vector<float> tparbeg(6);
	      float chisqbackwards = 0;
	      float lengthbackwards = 0;
	      retcode = KalmanFit(hitHandle,hitlist2,hsi,itrack,false,tparbeg,chisqbackwards,lengthbackwards,covmatbeg,unused_hits);
	      if (retcode != 0) continue;
              unused_hits.clear();  // todo -- use these hits for something

	      secondpass_tracks.emplace_back(lengthforwards,
					     lengthbackwards,
					     nhits,
					     tparbeg[5],
					     tparbeg.data(),
					     covmatbeg,
					     chisqforwards,
					     tparend[5],
					     tparend.data(),
					     covmatend,
					     chisqbackwards,
					     0);  // zero timestamp for now
	      secondpass_tid.push_back(itrack);
	    }
	}

      //  Remove stray hits.  Dig through unassociated hits and try to make extra tracks out of them.
      //  May need to wait until vertex finding is done so we know where to concentrate the effort

      // currently -- put second-pass tracks and associations with hits in the event

      if (fTrackPass == 1)
	{
	  for (size_t itrack=0; itrack<firstpass_tracks.size(); ++itrack)
	    {
	      trkCol->push_back(firstpass_tracks[itrack].CreateTrack());
	      auto const trackpointer = trackPtrMaker(itrack);
	      for (size_t ihit=0; ihit<hitlist[firstpass_tid[itrack]].size(); ++ ihit)
		{
		  auto const hitpointer = hitPtrMaker(hsi[hitlist[firstpass_tid[itrack]][ihit]]);
		  hitTrkAssns->addSingle(hitpointer,trackpointer);
		}
	    }
	}
      else if (fTrackPass == 2)
	{
	  for (size_t itrack=0; itrack<secondpass_tracks.size(); ++itrack)
	    {
	      trkCol->push_back(secondpass_tracks[itrack].CreateTrack());
	      auto const trackpointer = trackPtrMaker(itrack);
	      for (size_t ihit=0; ihit<hitlist2[secondpass_tid[itrack]].size(); ++ ihit)
		{
		  auto const hitpointer = hitPtrMaker(hsi[hitlist2[secondpass_tid[itrack]][ihit]]);
		  hitTrkAssns->addSingle(hitpointer,trackpointer);
		}
	    }
	}
      else
	{
	  throw cet::exception("Tracker1") << "Invalid track pass number requested: " << fTrackPass;
	}

      e.put(std::move(trkCol));
      e.put(std::move(hitTrkAssns));
    }

    //_____________________________________________________________________________
    float tracker1::capprox(float x1,float y1,
			    float x2,float y2,
			    float x3,float y3)
    {
      //-----------------------------------------------------------------
      // Initial approximation of the track curvature -- copied from ALICE
      // here x is y and y is z for us
      //-----------------------------------------------------------------
      x3 -=x1;
      x2 -=x1;
      y3 -=y1;
      y2 -=y1;
      //  
      float det = x3*y2-x2*y3;
      if (TMath::Abs(det)<1e-10){
	return 100;
      }
      //
      float u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
      float x0 = x3*0.5-y3*u;
      float y0 = y3*0.5+x3*u;
      float c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
      if (det<0) c2*=-1;
      return c2;
    }

    // redo this using y and z as coordinates, and also return a positive value.  The caller needs to decide the
    // screw orientation of the helix -- not enough to go on with just y and z coordinates

    float tracker1::capprox2(float y0, float z0, float y1, float z1, float y2, float z2)
    {
      float A = y1*(z2 - z0) - z1*(y2 - y0) + (y2*z0 - z2*y0);

      float B = -( (y1*y1 + z1*z1)*(z2 - z0)
		   -z1*( (y2*y2 + z2*z2) - (y0*y0 + z0*z0) )
		   + ( (y2*y2 + z2*z2)*z0 - z2*(y0*y0 + z0*z0) ) );

      float C = (y1*y1 + z1*z1)*(y2 - y0)
	-y1*( (y2*y2 + z2*z2) - (y0*y0 + z0*z0))
	+ ( (y2*y2 + z2*z2)*y0 - y2*(y0*y0 + z0*z0) );

      float D = - ( (y1*y1 + z1*z1)*(y2*z0 - z2*y0) 
		    - y1*( (y2*y2 + z2*z2)*z0 - z2*(y0*y0 + z0*z0) ) 
		    + z1*( (y2*y2 + z2*z2)*y0 - y2*(y0*y0 + z0*z0) ));

      if (TMath::Abs(A) < 1E-10) 
	{ return 0;}

      float yc = -B/(2.0*A);
      float zc = -C/(2.0*A);
      float rs = yc*yc + zc*zc -D/A;

      if (fPrintLevel > 1)
	{
	  std::cout << " In capprox2, A, B, C, D, rs: " << A << " " << B << " " << C << " " << D << " " << rs << std::endl;
	}
      if (rs <= 0) 
	{ throw cet::exception("tracker1capprox2: negative input to sqrt"); }
      float curv = 1.0/TMath::Sqrt(rs);
      return curv;
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    // KalmanFit does a forwards or backwards Kalman fit using the sorted hit list
    // variables:  x is the independent variable
    // 0: y
    // 1: z
    // 2: curvature
    // 3: phi
    // 4: slope

    int tracker1::KalmanFit( art::ValidHandle<std::vector<Hit> > &hitHandle, 
			     std::vector<std::vector<int> >      &hitlist,
			     std::vector<int>                    &hsi,
			     int itrack, 
			     bool isForwards, 
			     std::vector<float> &trackparatend,
			     float &chisquared,
			     float &length,
			     float *covmat,                     // 5x5 covariance matrix
			     std::vector<int> &unused_hits) 
    {

      // set some default values in case we return early

      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();
      chisquared = 0;
      length = 0;
      for (size_t i=0; i<5; ++i) trackparatend[i] = 0;
      for (size_t i=0; i<25; ++i) covmat[i] = 0;

      float roadsq = fRoadYZinFit*fRoadYZinFit;

      // estimate curvature, slope, phi, xpos from the initial track parameters
      float curvature_init=0.1;
      float phi_init = 0;
      float slope_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;
      if ( initial_trackpar_estimate(hitHandle, 
				     hitlist,
				     hsi,
				     itrack, 
				     isForwards,
				     curvature_init,
				     slope_init,
				     phi_init,
				     xpos_init,
				     ypos_init,
				     zpos_init) != 0)
	{
	  return 1;
	}

      // Kalman fitter variables

      float xpos = xpos_init; 

      TMatrixF P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      P[0][0] = TMath::Sq(1); // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1); // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5); // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // slope uncertainty
	  
      TMatrixF PPred(5,5);

      // per-step additions to the covariance matrix
      TMatrixF Q(5,5);
      Q.Zero();
      Q[2][2] = 0.0001;  // allow for some curvature uncertainty between points  (1%)
      Q[3][3] = 0.0001;   // phi
      Q[4][4] = 0.0001;   // slope

      // uncertainties on the measured points  (big for now)
      TMatrixF R(2,2);
      R.Zero();
      R[0][0] = TMath::Sq(fHitResolYZinFit);  // in cm^2
      R[1][1] = TMath::Sq(fHitResolYZinFit);  // in cm^2

      // add the hits and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
      TMatrixF F(5,5);
      TMatrixF FT(5,5);
      TVectorF parvec(5);
      parvec[0] = ypos_init;
      parvec[1] = zpos_init;
      parvec[2] = curvature_init;
      parvec[3] = phi_init;
      parvec[4] = slope_init;
      TVectorF predstep(5);

      TMatrixF H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
      TMatrixF HT(5,2);

      TVectorF z(2);
      TVectorF ytilde(2);
      TVectorF hx(2);
      TMatrixF S(2,2);
      TMatrixF K(5,2);

      TMatrixF I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;

      for (size_t ihit=1; ihit<nhits; ++ihit)
	{

	  size_t ihf = ifob(ihit,nhits,isForwards);
	  float xh = hits[hsi[hitlist[itrack][ihf]]].Position()[0];
	  float yh = hits[hsi[hitlist[itrack][ihf]]].Position()[1];
	  float zh = hits[hsi[hitlist[itrack][ihf]]].Position()[2];

	  if (fPrintLevel > 0)
	    {
	      std::cout << std::endl;
	      std::cout << "Adding a new hit: " << xh << " " << yh << " " << zh << std::endl;
	    }

	  // for readability

	  float curvature = parvec[2];
	  float phi = parvec[3];
	  float slope = parvec[4];
	  float dx = xh - xpos;

	  // update prediction to the plane containing x

	  F.Zero();

	  // y = yold + slope*dx*Sin(phi).   F[0][i] = df/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

	  F[0][0] = 1.; 
	  F[0][3] = dx*slope*TMath::Cos(phi);
	  F[0][4] = dx*TMath::Sin(phi);

	  // z = zold + slope*dx*Cos(phi)
	  F[1][1] = 1.;
	  F[1][3] = -dx*slope*TMath::Sin(phi);
	  F[1][4] = dx*TMath::Cos(phi);

	  // curvature = old curvature -- doesn't change but put in an uncertainty
	  F[2][2] = 1.;

	  // phi = old phi + curvature*slope*dx  
	  // need to take the derivative of a product here
	  F[3][2] = dx*slope;
	  F[3][3] = 1.;
	  F[3][4] = dx*curvature;

	  // slope -- same -- but put in an uncertainty in case it changes
	  F[4][4] = 1.;

	  // predicted step

	  if (fPrintLevel > 1)
	    {
	      std::cout << "F Matrix: " << std::endl;
	      F.Print();
	      std::cout << "P Matrix: " << std::endl;
	      P.Print();
	    }
	  if (fPrintLevel > 0)
	    {
	      std::cout << "x: " << xpos << " dx: " << dx <<  std::endl;
	      std::cout << " Parvec:   y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " slope " << parvec[4] << std::endl;
	    }

	  predstep = parvec;
	  predstep[0] += slope*dx*TMath::Sin(phi);  // update y
	  predstep[1] += slope*dx*TMath::Cos(phi);  // update z
	  predstep[3] += slope*dx*curvature;        // update phi

	  if (fPrintLevel > 1)
	    {
	      std::cout << " Predstep: y " << predstep[0] << " z " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " slope " << predstep[4] << std::endl;
	    }
	  // equations from the extended Kalman filter
	  FT.Transpose(F);
	  PPred = F*P*FT + Q;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "PPred Matrix: " << std::endl;
	      PPred.Print();
	    }

	  ytilde[0] = yh - predstep[0];
	  ytilde[1] = zh - predstep[1];
	  float ydistsq = ytilde.Norm2Sqr();
	  if (ydistsq > roadsq) 
	    {
	      unused_hits.push_back(ihf);
	      continue;
	    }
	  chisquared += ytilde.Norm2Sqr()/TMath::Sq(fHitResolYZ); 
	  if (fPrintLevel > 0)
	    {
	      std::cout << "ytilde (residuals): " << std::endl;
	      ytilde.Print();
	    }
	  if (fPrintLevel > 1)
	    {
	      std::cout << "H Matrix: " << std::endl;
	      H.Print();
	    }

	  HT.Transpose(H);
	  S = H*PPred*HT + R;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "S Matrix: " << std::endl;
	      S.Print();
	    }

	  S.Invert();
	  if (fPrintLevel > 1)
	    {
	      std::cout << "Inverted S Matrix: " << std::endl;
	      S.Print();
	    }

	  K = PPred*HT*S;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "K Matrix: " << std::endl;
	      K.Print();
	    }

	  float yprev = parvec[0];
	  float zprev = parvec[1];
	  parvec = predstep + K*ytilde;
	  P = (I-K*H)*PPred;
	  xpos = xh;

	  length += TMath::Sqrt( dx*dx + TMath::Sq(parvec[0]-yprev) + TMath::Sq(parvec[1]-zprev) );
	}

      for (size_t i=0; i<5; ++i)
	{
	  trackparatend[i] = parvec[i];
	}
      trackparatend[5] = xpos;  // tack this on so we can specify where the track endpoint is
      if (fPrintLevel > 1)
	{
	  std::cout << "Track params at end (y, z, curv, phi, slope) " << trackparatend[0] << " " << trackparatend[1] << " " <<
	    trackparatend[2] << " " << trackparatend[3] <<" " << trackparatend[4] << std::endl;
	  S.Print();
	}

      // just for visualization:  put in track parameter guess at beginning as the ending track parameter so we can
      // use the plotting tools -- remove or comment this out when the initial track parameters are validated
      //trackparatend[0] = trackbeg[1];
      //trackparatend[1] = trackbeg[2];
      //trackparatend[2] = curvature_init;
      //trackparatend[3] = phi_init;
      //trackparatend[4] = slope_init;
      //trackparatend[5] = trackbeg[0];


      size_t icov=0;
      for (size_t i=0; i<5; ++i)
	{
	  for (size_t j=0; j<5; ++j)
	    {
	      covmat[icov] = P[i][j];
	    }
	}

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------


    size_t tracker1::ifob(size_t ihit, size_t nhits, bool isForwards)
    {
      if (ihit >= nhits)
	{
	  throw cet::exception("Tracker1") << "Invalid hit index in ifob: " << ihit << " nhits: " << nhits;
	}
      if (isForwards)
	{
	  return ihit;
	}
      else
	{
	  return (nhits - ihit - 1);
	}
    }


    //--------------------------------------------------------------------------------------------------------------

    int tracker1::initial_trackpar_estimate(art::ValidHandle<std::vector<Hit> > &hitHandle, 
					    std::vector<std::vector<int> >      &hitlist,
					    std::vector<int>                    &hsi,
					    int itrack, 
					    bool isForwards,
					    float &curvature_init,
					    float &slope_init,
					    float &phi_init,
					    float &xpos,
					    float &ypos,
					    float &zpos)
    {
      // form a rough guess of track parameters

      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();

      size_t firsthit = ifob(0,nhits,isForwards);
      size_t inthit = ifob(fMinNumHits/2,nhits,isForwards);
      size_t farhit = ifob(fMinNumHits-1,nhits,isForwards);
      //size_t inthit = ifob(nhits/2,nhits,isForwards);
      //size_t farhit = ifob(nhits-1,nhits,isForwards);

      float trackbeg[3] = {hits[hsi[hitlist[itrack][firsthit]]].Position()[0],
			   hits[hsi[hitlist[itrack][firsthit]]].Position()[1],
			   hits[hsi[hitlist[itrack][firsthit]]].Position()[2]};

      float tp1[3] = {hits[hsi[hitlist[itrack][inthit]]].Position()[0],
		      hits[hsi[hitlist[itrack][inthit]]].Position()[1],
		      hits[hsi[hitlist[itrack][inthit]]].Position()[2]};

      float tp2[3] = {hits[hsi[hitlist[itrack][farhit]]].Position()[0],
		      hits[hsi[hitlist[itrack][farhit]]].Position()[1],
		      hits[hsi[hitlist[itrack][farhit]]].Position()[2]};

      if (fPrintLevel>1)
	{
	  std::cout << "Printing the first " << fMinNumHits << " hits" << std::endl;
	  for (size_t i=0;i<fMinNumHits;++i)
	    {
	      size_t ihf = ifob(i,nhits,isForwards);
	      std::cout << i << " : " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[0] << " " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[1] << " " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[2] << std::endl;
	    }
	}
      if (fPrintLevel>0)
	{
	  std::cout << "isForwards: " << isForwards << std::endl;
	  std::cout << "first hit: " << firsthit << ", inter hit: " << inthit << " " << " far hit: " << farhit << std::endl;
	  std::cout << "in the hit list: " << hsi[hitlist[itrack][firsthit]] << " " << hsi[hitlist[itrack][inthit]] << " " << hsi[hitlist[itrack][farhit]] << std::endl;
	  std::cout << "First hit x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
	  std::cout << "Inter hit x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
	  std::cout << "Far   hit x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
	}

      xpos = trackbeg[0];
      ypos = trackbeg[1];
      zpos = trackbeg[2];

      // linear guess for the initial slope.
      float dx1 = tp2[0] - xpos;
      float dyz1 = TMath::Sqrt( TMath::Sq( tp2[1] - trackbeg[1] ) +
				TMath::Sq( tp2[2] - trackbeg[2] ) );
      if (dx1 != 0) 
	{ slope_init = dyz1/dx1; }
      else
	{ 
	  slope_init = 0;
	  return 1; 
	} // got fMinNumHits all at exactly the same value of x (they were sorted).  Reject track.

      // phi in the y,z plane

      phi_init = TMath::ATan2( tp2[1]-trackbeg[1], tp2[2]-trackbeg[2] ); 

      curvature_init = TMath::Abs(capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2]));

      //float curvature_init2 = capprox2(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2]);
      //std::cout << "Compare initial curvatures: " << curvature_init << " " <<  curvature_init2 << std::endl;

      if (! isForwards) curvature_init *= -1;

      if (fPrintLevel>0)
	{
	  std::cout << "phi calc: dz, dy " << tp2[2]-trackbeg[2] << " " <<  tp2[1]-trackbeg[1] << std::endl;
	  std::cout << "initial curvature, phi, slope: " << curvature_init << " " << phi_init << " " << slope_init << std::endl;
	}
      return 0;
    }

    //--------------------------------------------------------------
    // the isForwards switch is only used to select which end to use to estimate the initial track parameters. Since this is a single helix
    // fit, the track parameters are the same, except for an evaluation of x, y, and z at the beginning and the end.

    int tracker1::FitHelix(art::ValidHandle<std::vector<Hit> > &hitHandle, 
			   std::vector<std::vector<int> > &hitlist, 
			   std::vector<int> &hsi, 
			   int itrack, 
			   bool isForwards, 
			   std::vector<int> &unused_hits,
			   TrackPar &trackpar
			   )
    {
      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();

      // estimate curvature, slope, phi, xpos from the initial track parameters
      float curvature_init=0.1;
      float phi_init = 0;
      float slope_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;
      if ( initial_trackpar_estimate(hitHandle, 
				     hitlist,
				     hsi,
				     itrack, 
				     isForwards,
				     curvature_init,
				     slope_init,
				     phi_init,
				     xpos_init,
				     ypos_init,
				     zpos_init) != 0)
	{
	  return 1;
	}

      float tpi[5] = {ypos_init, zpos_init, curvature_init, phi_init, slope_init};
      float covmat[25] = {0};

      // only need this to compute chisquared, so set the track parameters the same at the beginning and end

      TrackPar tpar(0,0,nhits,xpos_init,tpi,covmat,0,xpos_init,tpi,covmat,0,0);

      float c2sum = 0;
      for (size_t ihit=0; ihit<nhits; ++ihit)
	{
	  TVector3 hitpos(hits[hsi[hitlist[itrack][ihit]]].Position()[0],
		          hits[hsi[hitlist[itrack][ihit]]].Position()[1],
		          hits[hsi[hitlist[itrack][ihit]]].Position()[2]);
	  TVector3 helixpos = tpar.getPosAtX(hitpos.X(),isForwards);
	  c2sum += (hitpos-helixpos).Mag2();
	}

      // todo -- wrap a minimizer around this chisquared function to get the best track parameters
      return 0;
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar


