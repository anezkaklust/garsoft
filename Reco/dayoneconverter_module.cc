////////////////////////////////////////////////////////////////////////
// Class:       dayoneconverter
// Plugin Type: producer (art v3_00_00)
// File:        dayoneconverter_module.cc
//
// Generated at Tue Jun 23 12:28:55 2020 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
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
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include <memory>

#include "TMath.h"
#include "TVector3.h"

#include "Geant4/G4ThreeVector.hh"
#include "nutools/MagneticField/MagneticField.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandGauss.h"

// GArSoft Includes

#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "RecoAlg/TrackPropagator.h"

#include "Geometry/BitFieldCoder.h"
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/Geometry.h"

class dayoneconverter : public art::EDProducer {
public:
  explicit dayoneconverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  dayoneconverter(dayoneconverter const&) = delete;
  dayoneconverter(dayoneconverter&&) = delete;
  dayoneconverter& operator=(dayoneconverter const&) = delete;
  dayoneconverter& operator=(dayoneconverter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  std::string fInputEdepLabel;  ///<   Input label for edeps
  std::string fInputEdepInstanceTPC;  ///<   Input instance for TPC edeps
  std::string fInputEdepInstanceMuID;  ///<  Input instance for MuID edeps
  bool        fIncludeMuIDhits;       ///< Include MuID hits as TPCClusters
  float       fSmearX;         ///< amount by which to smear X, in cm
  float       fSmearY;         ///< amount by which to smear Y, in cm
  float       fSmearT;         ///< amount by which to smear T, in ns
  float       fZCut1;          ///< Cut to ensure TPC Clusters are on different planes
  float       fZCut2;          ///< Cut to ensure TPC clusters are on the same plane
  float       fRCut;           ///< Road in the YZ plane to add hits on a circle

  cet::exempt_ptr<CLHEP::HepRandomEngine> fEngine;  //< random engine

  const gar::geo::GeometryCore* fGeo;               ///< geometry information
  gar::geo::BitFieldCoder *fFieldDecoderTrk;
  gar::geo::BitFieldCoder *fFieldDecoderMuID;

  int makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar);

};


dayoneconverter::dayoneconverter(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
// More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  fInputEdepLabel = p.get<std::string>("InputLabel","edepconvert");
  fInputEdepInstanceTPC = p.get<std::string>("InputInstanceTPC","TrackerSc");
  fInputEdepInstanceMuID = p.get<std::string>("InputInstanceMuID","MuID");
  fIncludeMuIDhits = p.get<bool>("IncludeMuIDhits", false);
  fSmearX = p.get<float>("SmearX",0.3);  // in cm
  fSmearY = p.get<float>("SmearY",0.3);  // in cm
  fSmearT = p.get<float>("SmearT",1.0);  // in ns
  fZCut1   = p.get<float>("ZCut1",10.0); // in cm 
  fZCut2   = p.get<float>("ZCut2",10.0); // in cm 
  fRCut    = p.get<float>("RCut" ,5.0);  // in cm 

  art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
  art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
  consumes<std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
  consumes<std::vector<gar::sdp::CaloDeposit> >(muidedeptag);
  produces<std::vector<gar::rec::Track> >();
  produces<std::vector<gar::rec::TPCCluster> >();
  produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();

  // setup the random number service
  // obtain the random seed from NuRandomService
  ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "dayoneconverter", p, "RandomSeed");

  ::art::ServiceHandle<::art::RandomNumberGenerator> rng;
  auto& engine = rng->getEngine(art::ScheduleID::first(),
                                p.get<std::string>("module_label"),"dayoneconverter");
  fEngine = cet::make_exempt_ptr(&engine);

  fGeo = gar::providerFrom<gar::geo::Geometry>();
  std::string fEncoding = fGeo->GetMinervaCellIDEncoding();
  fFieldDecoderTrk = new gar::geo::BitFieldCoder( fEncoding );
  std::string fEncodingMuID = fGeo->GetMuIDCellIDEncoding();
  fFieldDecoderMuID = new gar::geo::BitFieldCoder( fEncodingMuID );
}

void dayoneconverter::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::unique_ptr< std::vector<gar::rec::TPCCluster> > TPCClusterCol(new std::vector<gar::rec::TPCCluster>);
  std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
  std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

  art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
  //art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
  auto tccdHandle = e.getValidHandle< std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
  auto const& tccds = *tccdHandle;

  // put these back when we have tracks and can associate them
  auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
  auto const tpcclusPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e);

  art::ServiceHandle<mag::MagneticField> magFieldService;
  G4ThreeVector zerovec(0,0,0);
  G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

  CLHEP::RandGauss GaussRand(*fEngine);

  // first stab, just make all the TPC clusters in the event, and then worry later
  // about pattern recognition

  for (auto const& cd : tccds)
    {
      float fcpos[3];
      //Need to check in which direction to smear (depends on the segmenetation)
      //Can use the cellID decoder to know looking at cellX and cellY values
      gar::raw::CellID_t cID = cd.CellID();
      int cellX = fFieldDecoderTrk->get(cID, "cellX");
      int cellY = fFieldDecoderTrk->get(cID, "cellY");

      fcpos[0] = cd.X(); // default values
      fcpos[1] = cd.Y();

      if(cellX == 1) {
        //Segmented in Y
        fcpos[0] = cd.X() + GaussRand.fire(0., fSmearX);
        fcpos[1] = cd.Y();
      }

      if(cellY == 1) {
        //Segmented in X
        fcpos[0] = cd.X();
        fcpos[1] = cd.Y() + GaussRand.fire(0., fSmearY);
      }

      fcpos[2] = cd.Z();

      float time = cd.Time() + GaussRand.fire(0., fSmearT);

      float covmat[6] = {0,0,0,0,0,0};  // TODO -- fill this in with something reasonble

      TPCClusterCol->emplace_back(cd.Energy(),   // do we need a scale factor here?
                                  fcpos,
                                  time,     // time is in ns
                                  time,
                                  time,
                                  fSmearX,
				  covmat);
    }

  if(fIncludeMuIDhits) {
    art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
    auto muIDHandle = e.getValidHandle< std::vector<gar::sdp::CaloDeposit> >(muidedeptag);
    auto const& muids = *muIDHandle;

    for (auto const& cd : muids)
      {
        float fcpos[3];
        //Need to check in which direction to smear (depends on the segmenetation)
        //Can use the cellID decoder to know looking at cellX and cellY values
        gar::raw::CellID_t cID = cd.CellID();
        int cellX = fFieldDecoderMuID->get(cID, "cellX");
        int cellY = fFieldDecoderMuID->get(cID, "cellY");

        fcpos[0] = cd.X(); // default values
        fcpos[1] = cd.Y();

        if(cellX == 1) {
          //Segmented in Y
          fcpos[0] = cd.X() + GaussRand.fire(0., fSmearX);
          fcpos[1] = cd.Y();
        }

        if(cellY == 1) {
          //Segmented in X
          fcpos[0] = cd.X();
          fcpos[1] = cd.Y() + GaussRand.fire(0., fSmearY);
        }

        fcpos[2] = cd.Z();

        float time = cd.Time() + GaussRand.fire(0., fSmearT);

        float covmat[6] = {0,0,0,0,0,0};  // TODO -- fill this in with something reasonble

        TPCClusterCol->emplace_back(cd.Energy(),   // do we need a scale factor here?
                                    fcpos,
                                    time,     // time is in ns
                                    time,
                                    time,
                                    fSmearX,
				    covmat);
      }
  }

  if (true)
    {

  // sort the TPC Clusters along Z

  std::sort(TPCClusterCol->begin(),TPCClusterCol->end(),
	    [](const gar::rec::TPCCluster &a, const gar::rec::TPCCluster &b)->bool
	    { return a.Position()[2] < b.Position()[2]; } );

  size_t ntpcclus = TPCClusterCol->size();

  // look for best-fit tracks in the list of TPC clusters
  // find the best triplet of TPC Clusters with spacing at least fZCut that fit on a circle
  // to think about -- time and energy cuts

  size_t bestnpts = 0;
  float bestsumr2 = 0;
  std::vector<size_t> besttpcclusindex;    
  gar::rec::Track besttrack;

  for (size_t i=0; i<ntpcclus; ++i)
    {
      const float *ipos = TPCClusterCol->at(i).Position();
      for (size_t j=i+1; j<ntpcclus; ++j)
        {
          const float *jpos = TPCClusterCol->at(j).Position();
          if (TMath::Abs( ipos[2] - jpos[2] ) < fZCut1) continue;
          for (size_t k=j+1; k<ntpcclus; ++k)
            {
              const float *kpos = TPCClusterCol->at(k).Position();
              if (TMath::Abs( ipos[2] - kpos[2] ) < fZCut1) continue;
              if (TMath::Abs( jpos[2] - kpos[2] ) < fZCut1) continue;
	      std::vector<gar::rec::TPCCluster> triplet;
	      triplet.push_back(TPCClusterCol->at(i));
	      triplet.push_back(TPCClusterCol->at(j));
	      triplet.push_back(TPCClusterCol->at(k));
	      gar::rec::TrackPar triplettrack;
	      makepatrectrack(triplet,triplettrack);

	      // pick the best TPC clusters at each Z position.  Save their
	      // indices in tpcclusindex and sum the squares of distances in sumr2

	      std::vector<size_t> tpcclusindex;
	      float sumr2 = 0;

	      float zcur = -2E9;
	      int tpcclusindexb = -1;
	      float dbest=1E9;
	      gar::rec::Track tpt;

	      for (size_t k2=0; k2<ntpcclus; ++k2)
		{
		  const float *k2pos = TPCClusterCol->at(k2).Position();

		  // clusters are sorted along Z.  If we found a new Z, put the best point on the list

		  if ((TMath::Abs(zcur - k2pos[2]) > fZCut2) &&
		      (tpcclusindexb > -1) &&
		      (dbest < fRCut) )
		    {
		      tpcclusindex.push_back(tpcclusindexb);
		      sumr2 += dbest*dbest;
		      dbest = 1E9;
		      tpcclusindexb = -1;
		      zcur = k2pos[2];
		    }

		  float dist=0;
		  tpt = triplettrack.CreateTrack();
		  int retcode = util::TrackPropagator::DistXYZ(tpt.TrackParBeg(),tpt.Vertex(),k2pos,dist);
		  if (retcode != 0) continue;
		  if (dist > fRCut) continue;
		  if (dist<dbest)
		    {
		      dbest = dist;
		      tpcclusindexb = k2;
		    }
		  // last point -- check to see if it gets added.
		  if (k2 == ntpcclus-1 && tpcclusindexb > -1)
		    {
		      tpcclusindex.push_back(tpcclusindexb);
		      sumr2 += dbest*dbest;
		    }		  
		}  // end loop over k2 -- assigning clusters to this track

	      if (tpcclusindex.size() > bestnpts || 
		  ((tpcclusindex.size() == bestnpts) && 
		   (sumr2 < bestsumr2)))
		{
		  bestnpts = tpcclusindex.size();
		  bestsumr2 = sumr2;
		  besttpcclusindex = tpcclusindex;
		  besttrack = tpt;
		}
	    } // end loop over k in triplet
	} // end loop over j in triplet
    } // end loop over i in triplet

  // so far can only make one track.  Look at other points in collection not yet used and make more tracks

  if (bestnpts > 0)
    {
      trkCol->push_back(besttrack);
      auto const trackpointer = trackPtrMaker(trkCol->size()-1);
      for (size_t i=0; i<besttpcclusindex.size(); ++i)
	{
	  auto const tpccluspointer = tpcclusPtrMaker(besttpcclusindex.at(i));
	  TPCClusterTrkAssns->addSingle(tpccluspointer,trackpointer);
	}
    }

    }
  e.put(std::move(trkCol));
  e.put(std::move(TPCClusterCol));
  e.put(std::move(TPCClusterTrkAssns));
}

// maybe refactor this so we don't have to duplicate it.  But need to pass in all the config parameters.
// temporary: just hardcode the parameters to see if we can get something to work

int dayoneconverter::makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar)
{
  // track parameters:  x is the independent variable
  // 0: y
  // 1: z
  // 2: curvature
  // 3: phi
  // 4: lambda = angle from the cathode plane
  // 5: x   /// added on to the end

  float fSortTransWeight = 1.0;
  int fInitialTPNTPCClusters = 100;
  float fSortDistBack = 2.0;
  int fPrintLevel = 0;

  float lengthforwards = 0;
  std::vector<int> hlf;
  float lengthbackwards = 0;
  std::vector<int> hlb;

  gar::rec::sort_TPCClusters_along_track(trackTPCClusters,hlf,hlb,fPrintLevel,lengthforwards,lengthbackwards,fSortTransWeight,fSortDistBack);

  std::vector<float> tparbeg(6,0);
  float xother = 0;
  if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlf, tparbeg[2], tparbeg[4], 
					   tparbeg[3], tparbeg[5], tparbeg[0], tparbeg[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0) 
    {
      return 1;
    }

  std::vector<float> tparend(6,0);
  if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlb, tparend[2], tparend[4], 
					   tparend[3], tparend[5], tparend[0], tparend[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0)
    {
      return 1;
    }

  // no chisquare or covariance in patrec tracks
  float covmatbeg[25];
  float covmatend[25];
  for (size_t i=0; i<25; ++i) // no covmat in patrec tracks
    {
      covmatend[i] = 0;
      covmatbeg[i] = 0;
    }

  trackpar.setNTPCClusters(trackTPCClusters.size());
  trackpar.setTime(0);
  trackpar.setChisqForwards(0);
  trackpar.setChisqBackwards(0);
  trackpar.setLengthForwards(lengthforwards);
  trackpar.setLengthBackwards(lengthbackwards);
  trackpar.setCovMatBeg(covmatbeg);
  trackpar.setCovMatEnd(covmatend);
  trackpar.setTrackParametersBegin(tparbeg.data());
  trackpar.setXBeg(tparbeg[5]);
  trackpar.setTrackParametersEnd(tparend.data());
  trackpar.setXEnd(tparend[5]);

  return 0;
}


DEFINE_ART_MODULE(dayoneconverter)
