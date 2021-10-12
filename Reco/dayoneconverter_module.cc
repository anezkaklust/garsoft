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
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include <memory>

#include "TMath.h"
#include "TVector3.h"
#include "TF1.h"

#include "Geant4/G4ThreeVector.hh"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandGauss.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft Includes

#include "DetectorInfo/GArMagneticField.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "RecoAlg/TrackPropagator.h"

#include "Geometry/BitFieldCoder.h"
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/GeometryGAr.h"


namespace gar {
  namespace rec {

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
      float       fPECm;           ///< conversion factor from cm step length to pe
      float       fSmearLY;        ///< amount by which to smear the LY
      float       fThrPE;          ///< threshold cut in pe
      float       fZCut1;          ///< Cut to ensure TPC Clusters are on different planes
      float       fZCut2;          ///< Cut to ensure TPC clusters are on the same plane
      float       fRCut;           ///< Road in the YZ plane to add hits on a circle
      float       fTimeBunch;      ///< Time bunch spread, in ns            

      bool        fUseOnlyTrueMuonHits;   ///< whether to cheat and use true muon hits
      std::string fG4ModuleLabel;     ///< Use this to get the MCParticles

      CLHEP::HepRandomEngine &fEngine;  //< random engine

      const gar::geo::GeometryCore* fGeo;               ///< geometry information
      gar::geo::BitFieldCoder *fFieldDecoderTrk;
      gar::geo::BitFieldCoder *fFieldDecoderMuID;
      TF1* fMu2e;

      int makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar);
      void digitizeCaloHitsSimple(gar::sdp::CaloDeposit cd, float *fcpos, float &energy, float &time);
      void digitizeCaloHitsMu2e(gar::sdp::CaloDeposit cd, float *fcpos, float &energy, float &time);
    };

    dayoneconverter::dayoneconverter(fhicl::ParameterSet const& p)
      : EDProducer{p}
      , fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,p,"Seed"))
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
        fPECm  = p.get<float>("PeCm", 50.);  // in pe/cm
        fSmearLY  = p.get<float>("SmearLY", 10.0);  // in pe
        fThrPE    = p.get<float>("ThrPE", 5.0);  // in pe
        fZCut1    = p.get<float>("ZCut1",10.0); // in cm
        fZCut2    = p.get<float>("ZCut2",0.5); // in cm
        fRCut     = p.get<float>("RCut" ,5.0);  // in cm
        fTimeBunch = p.get<float>("TimeBunch",150000); // in ns            
	fUseOnlyTrueMuonHits = p.get<bool>("UseOnlyTrueMuonHits",false);  // whether to cheat or not
	fG4ModuleLabel = p.get<std::string>("G4ModuleLabel","geant");  // for retrieving MCParticles

        art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
        art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
        consumes<std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
        consumes<std::vector<gar::sdp::CaloDeposit> >(muidedeptag);
	if (fUseOnlyTrueMuonHits)
	  {
	    consumes<std::vector<simb::MCParticle> >(fG4ModuleLabel);
	  }
        produces<std::vector<gar::rec::Track> >();
        produces<std::vector<gar::rec::TPCCluster> >();
        produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();

        fGeo = gar::providerFrom<gar::geo::GeometryGAr>();
        std::string fEncoding = fGeo->GetMinervaCellIDEncoding();
        fFieldDecoderTrk = new gar::geo::BitFieldCoder( fEncoding );
        std::string fEncodingMuID = fGeo->GetMuIDCellIDEncoding();
        fFieldDecoderMuID = new gar::geo::BitFieldCoder( fEncodingMuID );

        fMu2e = new TF1("mu2e_pe", "[0] + [1]*x", 0, 600);
        fMu2e->SetParameter(0, 50.);//pe
        fMu2e->SetParameter(1, -0.06);//pe/cm
      }

    void dayoneconverter::produce(art::Event& e)
    {
      // Implementation of required member function here.
      std::unique_ptr< std::vector<gar::rec::TPCCluster> > TPCClusterCol(new std::vector<gar::rec::TPCCluster>);
      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      art::InputTag tpcedeptag(fInputEdepLabel,fInputEdepInstanceTPC);
      auto tccdHandle = e.getValidHandle< std::vector<gar::sdp::CaloDeposit> >(tpcedeptag);
      auto const& tccds = *tccdHandle;

      // put these back when we have tracks and can associate them
      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const tpcclusPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e);

      auto const *magFieldService = gar::providerFrom<mag::MagneticFieldService>();
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::Handle<std::vector<simb::MCParticle > > mcphandle;
      std::unordered_map<int, const simb::MCParticle * > TrkIDMap;
      if (fUseOnlyTrueMuonHits)
	{
	  mcphandle = e.getHandle<std::vector<simb::MCParticle> >(fG4ModuleLabel);
	  if (!mcphandle)
	    {
	       throw cet::exception("DayOneConverter::Produce")
                    << "Attempting to identify muon hits without MC particles";
	    }
	  for (auto const &mcp : *mcphandle)
	    {
	      TrkIDMap[mcp.TrackId()] = &mcp;
	    }
	}

      // first stab, just make all the TPC clusters in the event, and then worry later
      // about pattern recognition

      for (auto const& cd : tccds)
        {
	  const int trackid = cd.TrackID();
          if (fUseOnlyTrueMuonHits)
	    {
	      auto mcpt = TrkIDMap.find(trackid);
	      if (mcpt == TrkIDMap.end()) continue;
	      if (std::abs(mcpt->second->PdgCode()) != 13) continue;
	    }

          float energy = 0.;
          float time = 0.;
          float fcpos[3] = {0., 0., 0.};
	  
          // digitizeCaloHitsSimple(cd, fcpos, energy, time);
          digitizeCaloHitsMu2e(cd, fcpos, energy, time);

          float covmat[6] = {0,0,0,0,0,0};  // TODO -- fill this in with something reasonble

          if(energy <= 0.) continue;

          TPCClusterCol->emplace_back(energy,
                                      fcpos,
                                      time,     // time is in ns
                                      time,
                                      time,
                                      fSmearX,
                                      covmat);
        }

      if(fIncludeMuIDhits)
        {
          art::InputTag muidedeptag(fInputEdepLabel,fInputEdepInstanceMuID);
          auto muIDHandle = e.getValidHandle< std::vector<gar::sdp::CaloDeposit> >(muidedeptag);
          auto const& muids = *muIDHandle;

          for (auto const& cd : muids)
            {

	      const int trackid = cd.TrackID();
              if (fUseOnlyTrueMuonHits)
	        {
	          auto mcpt = TrkIDMap.find(trackid);
	          if (mcpt == TrkIDMap.end()) continue;
	          if (std::abs(mcpt->second->PdgCode()) != 13) continue;
	        }

              float energy = 0.;
              float time = 0.;
              float fcpos[3] = {0., 0., 0.};

              // digitizeCaloHitsSimple(cd, fcpos, energy, time);
              digitizeCaloHitsMu2e(cd, fcpos, energy, time);

              if(energy <= 0.) continue;

              float covmat[6] = {0,0,0,0,0,0};  // TODO -- fill this in with something reasonble

              TPCClusterCol->emplace_back(energy,
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

          // sort the TPC Clusters based on time
          bool debug=false;

          std::vector<std::vector<gar::rec::TPCCluster>> TPCClusterColTime;
          std::vector<gar::rec::TPCCluster> TPCClusterColTimeBlock;
          std::sort(TPCClusterCol->begin(),TPCClusterCol->end(),
                    [](const gar::rec::TPCCluster &a, const gar::rec::TPCCluster &b)->bool
                    { return a.Time() < b.Time(); } );
          size_t ntpcclus = TPCClusterCol->size();

          if (ntpcclus!=0&&debug) std::cout<<"Time start and end: "<<TPCClusterCol->at(0).Time()<<" "<<TPCClusterCol->at(TPCClusterCol->size()-1).Time()<<std::endl;

          float startt= 0;
          if (ntpcclus!=0)
            {
              startt=TPCClusterCol->at(0).Time();
              for(size_t i=0; i<ntpcclus; i++)
                {
                        
                  if(TPCClusterCol->at(i).Time()-startt<=fTimeBunch)
                    {
                      TPCClusterColTimeBlock.push_back(TPCClusterCol->at(i));
                    }
                  else
                    {
                      TPCClusterColTime.push_back(TPCClusterColTimeBlock);
                      startt=TPCClusterCol->at(i).Time();
                      TPCClusterColTimeBlock.clear();     
                      TPCClusterColTimeBlock.push_back(TPCClusterCol->at(i));             
                    }
                        
                }

              if(TPCClusterColTimeBlock.size()!=0) TPCClusterColTime.push_back(TPCClusterColTimeBlock);
            }


          if(debug)
            {
              for(size_t c=0; c<TPCClusterColTime.size(); c++)
                {
                    
                  std::cout << "l"<<c<<" = { ";
                  for (size_t d=0; d<TPCClusterColTime.at(c).size(); d++) {
                    std::cout << TPCClusterColTime.at(c).at(d).Time() << ", ";
                  }
                  std::cout << "};\n";
                
                }
            }
                
          for(size_t t=0; t<TPCClusterColTime.size(); t++)
            {

              ////Sort the clusters along Z
              std::sort(TPCClusterColTime.at(t).begin(),TPCClusterColTime.at(t).end(),
                        [](const gar::rec::TPCCluster &a, const gar::rec::TPCCluster &b)->bool
                        { return a.Position()[2] < b.Position()[2]; } );


              size_t ntpcclust = TPCClusterColTime.at(t).size();
                    
                    
              // look for best-fit tracks in the list of TPC clusters
              // find the best triplet of TPC Clusters with spacing at least fZCut that fit on a circle
              // to think about -- time and energy cuts

              size_t bestnpts = 0;
              float bestsumr2 = 0;
              std::vector<size_t> besttpcclusindex;
              gar::rec::Track besttrack;

              for (size_t i=0; i<ntpcclust; ++i)
                {
                  const float *ipos = TPCClusterColTime.at(t).at(i).Position();
                  for (size_t j=i+1; j<ntpcclust; ++j)
                    {
                      const float *jpos = TPCClusterColTime.at(t).at(j).Position();
                      if (TMath::Abs( ipos[2] - jpos[2] ) < fZCut1) continue;
                      for (size_t k=j+1; k<ntpcclust; ++k)
                        {
                          const float *kpos = TPCClusterColTime.at(t).at(k).Position();
                          if (TMath::Abs( ipos[2] - kpos[2] ) < fZCut1) continue;
                          std::vector<gar::rec::TPCCluster> triplet;
                          triplet.push_back(TPCClusterColTime.at(t).at(i));
                          triplet.push_back(TPCClusterColTime.at(t).at(j));
                          triplet.push_back(TPCClusterColTime.at(t).at(k));
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

                          for (size_t k2=0; k2<ntpcclust; ++k2)
                            {
                              const float *k2pos = TPCClusterColTime.at(t).at(k2).Position();

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
                              if (k2 == ntpcclust-1 && tpcclusindexb > -1)
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

              // so far we can only make one track.  Look at other points in collection not yet used and make more tracks

              if (bestnpts > 0)
                {
                  // "besttrack" above only has track parameters from the triplet.  make a new track from
                  // all the TPC clusters
                  //trkCol->push_back(besttrack);
                  std::vector<gar::rec::TPCCluster> tcv;
                  for (size_t i=0;i<besttpcclusindex.size(); ++i)
                    {
                      tcv.push_back(TPCClusterColTime.at(t).at(besttpcclusindex.at(i)));
                    }

                  gar::rec::TrackPar btp;
                  makepatrectrack(tcv,btp);
                  gar::rec::Track btt = btp.CreateTrack();
                  trkCol->push_back(btt);

                  auto const trackpointer = trackPtrMaker(trkCol->size()-1);
                  for (size_t i=0; i<besttpcclusindex.size(); ++i)
                    {
                      auto const tpccluspointer = tpcclusPtrMaker(besttpcclusindex.at(i));
                      TPCClusterTrkAssns->addSingle(tpccluspointer,trackpointer);
                    }
                }

            }
        }
      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterCol));
      e.put(std::move(TPCClusterTrkAssns));
            
    }

    // digitize the plane hits based on minerva numbers

    void dayoneconverter::digitizeCaloHitsSimple(gar::sdp::CaloDeposit cd, float *fcpos, float &energy, float &time)
    {
      //Need to check in which direction to smear (depends on the segmenetation)
      //Can use the cellID decoder to know looking at cellX and cellY values
      CLHEP::RandGauss GaussRand(fEngine);

      gar::raw::CellID_t cID = cd.CellID();
      int cellX = fFieldDecoderTrk->get(cID, "cellX");
      int cellY = fFieldDecoderTrk->get(cID, "cellY");

      fcpos[0] = cd.X(); // default values
      fcpos[1] = cd.Y();
      fcpos[2] = cd.Z();
      energy = cd.Energy();
      time = cd.Time();

      if(cellX == 0) {
        //Segmented in Y
        fcpos[0] += GaussRand.fire(0., fSmearX);
      }

      if(cellY == 0) {
        //Segmented in X
        fcpos[1] += GaussRand.fire(0., fSmearY);
      }

      time += GaussRand.fire(0., fSmearT);

      //convert energy to pe
      energy = fPECm * cd.StepLength();
      energy += GaussRand.fire(0., fSmearLY);

      if(energy < fThrPE) energy = 0.;

      return;
    }

    //digitize based on Mu2e strip numbers

    void dayoneconverter::digitizeCaloHitsMu2e(gar::sdp::CaloDeposit cd, float *fcpos, float &energy, float &time)
    {
      //Need to check in which direction to smear (depends on the segmenetation)
      //Can use the cellID decoder to know looking at cellX and cellY values
      CLHEP::RandGauss GaussRand(fEngine);

      gar::raw::CellID_t cID = cd.CellID();
      int cellX = fFieldDecoderTrk->get(cID, "cellX");
      int cellY = fFieldDecoderTrk->get(cID, "cellY");

      fcpos[0] = cd.X(); // default values
      fcpos[1] = cd.Y();
      fcpos[2] = cd.Z();
      energy = cd.Energy();
      time = cd.Time();

      if(cellX == 0) {
        //Segmented in Y
        fcpos[0] += GaussRand.fire(0., fSmearX);
      }

      if(cellY == 0) {
        //Segmented in X
        fcpos[1] += GaussRand.fire(0., fSmearY);
      }

      time += GaussRand.fire(0., fSmearT);

      //convert energy to pe based on the distance in the strip
      //Curve as "linear" with parameters a = 0.06 pe/cm, b = 50 pe (80*0.6 from Mu2e)
      double local_distance = 0.;
      std::array<double, 3> point = {cd.X(), cd.Y(), cd.Z()};
      std::array<double, 3> pointLocal;
      gar::geo::LocalTransformation<TGeoHMatrix> trans;
      fGeo->WorldToLocal(point, pointLocal, trans);
      TVector3 tpoint(point[0], point[1], point[2]);
      std::string name = fGeo->VolumeName(tpoint);

      if(cellX == 0) {
        //Segmented in Y -> get the x pos
        local_distance = 300. - pointLocal[0];
      }

      if(cellY == 0) {
        //Segmented in X -> get the y pos
        local_distance = 300. - pointLocal[1];
      }

      energy = fMu2e->Eval(local_distance);
      // std::cout << "Volume " << name << std::endl;
      // std::cout << "CellX " << cellX << " CellY " << cellY << std::endl;
      // std::cout << "Global Position " << cd.X() << ", " << cd.Y() << ", " << cd.Z() << std::endl;
      // std::cout << "Local Position " << pointLocal[0] << ", " << pointLocal[1] << ", " << pointLocal[2] << std::endl;
      // std::cout << "local distance in strip " << local_distance << " cm has " << energy << " pe detected" << std::endl;
      energy += GaussRand.fire(0., fSmearLY);

      if(energy < fThrPE) energy = 0.;

      return;
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
      float chisqforwards = 0;
      float tvpos[3] = {tparbeg[5],tparbeg[0],tparbeg[1]};
      for (size_t i=0; i<trackTPCClusters.size(); ++i)
        {
          float dist = 0;
          const float *pos = trackTPCClusters.at(i).Position();
          int retcode = util::TrackPropagator::DistXYZ(tparbeg.data(),tvpos,pos,dist);
          if (retcode == 0)
            {
              chisqforwards += dist*dist/(fSmearX*fSmearX);   // crude parameterization of errors
            }
        }

      std::vector<float> tparend(6,0);
      if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlb, tparend[2], tparend[4],
                                               tparend[3], tparend[5], tparend[0], tparend[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0)
        {
          return 1;
        }
      float chisqbackwards = 0;
      float tepos[3] = {tparend[5],tparend[0],tparend[1]};
      for (size_t i=0; i<trackTPCClusters.size(); ++i)
        {
          float dist = 0;
          const float *pos = trackTPCClusters.at(i).Position();
          int retcode = util::TrackPropagator::DistXYZ(tparend.data(),tepos,pos,dist);
          if (retcode == 0)
            {
              chisqbackwards += dist*dist/(fSmearX*fSmearX);   // crude parameterization of errors
            }
        }


      // no covariance matrix in patrec tracks
      float covmatbeg[25];
      float covmatend[25];
      for (size_t i=0; i<25; ++i) // no covmat in patrec tracks
        {
          covmatend[i] = 0;
          covmatbeg[i] = 0;
        }

      trackpar.setNTPCClusters(trackTPCClusters.size());
      trackpar.setTime(0);
      trackpar.setChisqForwards(chisqforwards);
      trackpar.setChisqBackwards(chisqbackwards);
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
  }
}

namespace gar {
  namespace rec {
    DEFINE_ART_MODULE(dayoneconverter)
  }
}
