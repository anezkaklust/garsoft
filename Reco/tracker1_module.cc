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

      float fHitResolYZ;           ///< resolution in cm of a hit in YZ (pad size)
      float fHitResolX;            ///< resolution in cm of a hit in X (drift direction)
      float fSigmaRoad;            ///< how many sigma away from a track a hit can be and still add it
      float fXGapToEndTrack;       ///< how big a gap must be before we end a track and start a new one
      unsigned int fMinNumHits;    ///< minimum number of hits to define a track
      std::string fHitLabel;       ///< label of module creating hits
      int fPrintLevel;             ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all

      int KalmanFit( art::ValidHandle<std::vector<Hit> > &hitHandle, 
		     std::vector<std::vector<int> > &hitlist, 
		     std::vector<int> &hsi, 
		     int itrack, 
		     bool Forwards, 
		     std::vector<float> &trackparatend,
		     float &chisquared,
		     float &length);

      size_t ifob(size_t ihit, size_t nhits, bool Forwards);
 
      float capprox(float x1,float y1,
		    float x2,float y2,
		    float x3,float y3);  ///< initial guess of curvature calculator

    };


    tracker1::tracker1(fhicl::ParameterSet const & p)
    {

      produces< std::vector<rec::Track> >();
      produces< ::art::Assns<rec::Hit, rec::Track> >();

      fHitResolYZ     = p.get<float>("HitResolYZ",1.0); // TODO -- think about what this value is
      fHitResolX      = p.get<float>("HitResolX",0.5);  // this is probably much better
      fSigmaRoad      = p.get<float>("SigmaRoad",5.0);
      fMinNumHits     = p.get<unsigned int>("MinNumHits",20);
      fHitLabel       = p.get<std::string>("HitLabel","hit");
      fPrintLevel     = p.get<int>("PrintLevel",0);
    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< ::art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

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
      std::vector< std::vector<int> > hitlist;

      float resolSq = fHitResolYZ*fHitResolYZ;

      // idea -- find all of the track candidates a hit can be added to that have signfs <= roadsq, and
      // pick the one that is closest to the local linear extrapolation of the existing hits.

      for (size_t i=0; i<hits.size(); ++i)
	{
	  const float *hpos = hits[hsi[i]].Position();
	  float bestsignifs = -1;
	  int ibest = -1;
	  for (size_t itcand = 0; itcand < hitlist.size(); ++itcand)
	    {
	      const float *cpos = hits[hsi[hitlist[itcand].back()]].Position();
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
      // fit them in both directions.  The Kalman filter gives the most precise measurement
      // of position and momentum at the end.

      size_t ntracks = hitlist.size();
      for (size_t itrack=0; itrack<ntracks; ++itrack)
	{
	  size_t nhits = hitlist[itrack].size();
	  if ( nhits >= fMinNumHits)
	    {

	      if (fPrintLevel)
		{
		  std::cout << "Starting a new track: " << itrack << " Number of hits from patrec: " << nhits << std::endl;
		}

	      // variables:  x is the independent variable
	      // 0: y
	      // 1: z
	      // 2: curvature
	      // 3: phi
	      // 4: slope = d(yz distance)/dx
	      // 5: x   /// added on to the end

	      std::vector<float> tparend(6);
	      float chisqforwards = 0;
	      float lengthforwards = 0;
	      int retcode = KalmanFit(hitHandle,hitlist,hsi,itrack,true,tparend,chisqforwards,lengthforwards);
	      if (retcode != 0) continue;

	      std::vector<float> tparbeg(6);
	      float chisqbackwards = 0;
	      float lengthbackwards = 0;
	      retcode = KalmanFit(hitHandle,hitlist,hsi,itrack,false,tparbeg,chisqbackwards,lengthbackwards);
	      if (retcode != 0) continue;

	      float length = 0.5*(lengthforwards + lengthbackwards);

	      float momentum_beg = 0;
	      if (tparbeg[2] != 0) momentum_beg = 0.3*magfield[0]/tparbeg[2];  // TODO -- check constant?  In kGauss or T?
	      float momentum_end = 0;
	      if (tparend[2] != 0) momentum_end = 0.3*magfield[0]/tparend[2];

	      float posbeg[3] = {tparbeg[5],tparbeg[0],tparbeg[1]};
	      float posend[3] = {tparend[5],tparend[0],tparend[1]};
	      
	      float hbeg = TMath::Sqrt( 1.0 + TMath::Sq(tparbeg[4]) );
	      float sbeg = TMath::Sin(tparbeg[3]);
	      float cbeg = TMath::Cos(tparbeg[3]);
	      float dirbeg[3] = {tparbeg[4]/hbeg,
				 cbeg/hbeg,
				 sbeg/hbeg};

	      float hend = TMath::Sqrt( 1.0 + TMath::Sq(tparend[4]) );
	      float send = TMath::Sin(tparend[3]);
	      float cend = TMath::Cos(tparend[3]);
	      float dirend[3] = {tparend[4]/hend,
				 cend/hend,
				 send/hend};

	      trkCol->emplace_back(length,
				  momentum_beg,
				  momentum_end,
				  posbeg,
				  posend,
				  dirbeg,
				  dirend,
				  chisqbackwards,
				  chisqforwards,
				  nhits);

	      auto const trackpointer = trackPtrMaker(trkCol->size()-1);

	      for (size_t ihit=0; ihit<nhits; ++ ihit)
		{
		  auto const hitpointer = hitPtrMaker(hsi[hitlist[itrack][ihit]]);
		  hitTrkAssns->addSingle(hitpointer,trackpointer);
		}
	    }
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


    // does a forwards or backwards Kalman fit using the sorted hit list
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
			     bool Forwards, 
			     std::vector<float> &trackparatend,
			     float &chisquared,
			     float &length)
    {

      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();
      chisquared = 0;
      length = 0;

      // form a rough guess of track parameters

      size_t firsthit = ifob(0,nhits,Forwards);
      size_t inthit = ifob(fMinNumHits/2,nhits,Forwards);
      size_t farhit = ifob(fMinNumHits,nhits,Forwards);

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
	      size_t ihf = ifob(i,nhits,Forwards);
	      std::cout << i << " : " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[0] << " " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[1] << " " << 
		hits[hsi[hitlist[itrack][ihf]]].Position()[2] << std::endl;
	    }
	}
      if (fPrintLevel>0)
	{
	  std::cout << "Forwards: " << Forwards << std::endl;
	  std::cout << "first hit: << " << firsthit << ", inter hit: " << inthit << " " << " far hit: " << farhit << std::endl;
	  std::cout << "in the hit list: " << hsi[hitlist[itrack][firsthit]] << " " << hsi[hitlist[itrack][inthit]] << " " << hsi[hitlist[itrack][farhit]] << std::endl;
	  std::cout << "First hit x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
	  std::cout << "Inter hit x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
	  std::cout << "Far   hit x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
	}

      float curvature_init = 0.1; // initial guess -- maybe use the min and max above to refine this guess
      float phi_init = 0;        // initial guess at the track beginning
      float slope_init = 1;   // d(yz distance)/dx

      float xpos = trackbeg[0];

      // linear guess for the initial slope.
      float dx1 = tp2[0] - xpos;
      float dyz1 = TMath::Sqrt( TMath::Sq( tp2[1] - trackbeg[1] ) +
				TMath::Sq( tp2[2] - trackbeg[2] ) );
      if (dx1 != 0) 
	{ slope_init = dyz1/dx1; }
      else
	{ return 1; } // got fMinNumHits all at exactly the same value of x (they were sorted).  Reject track.
      // phi in the y,z plane
      phi_init = TMath::ATan2( tp2[2]-trackbeg[2], tp2[1]-trackbeg[1] ); 
      curvature_init = capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2]);

      if (fPrintLevel>0)
	{
	  std::cout << "initial slope, phi, curvature: " << slope_init << " " << phi_init << " " << curvature_init << std::endl;
	}

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
      R[0][0] = TMath::Sq(fHitResolYZ);  // in cm^2
      R[1][1] = TMath::Sq(fHitResolYZ);  // in cm^2

      // add the hits and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
      TMatrixF F(5,5);
      TMatrixF FT(5,5);
      TVectorF parvec(5);
      parvec[0] = trackbeg[1];
      parvec[1] = trackbeg[2];
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

	  size_t ihf = ifob(ihit,nhits,Forwards);
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

	  // y = yold + slope*dx*Sin(phi)
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
      return 0;
    }

    size_t tracker1::ifob(size_t ihit, size_t nhits, bool Forwards)
    {
      if (ihit >= nhits)
	{
	  throw cet::exception("Tracker1") << "Invalid hit index in ifob: " << ihit << " nhits: " << nhits;
	}
      if (Forwards)
	{
	  return ihit;
	}
      else
	{
	  return (nhits - ihit - 1);
	}
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar


