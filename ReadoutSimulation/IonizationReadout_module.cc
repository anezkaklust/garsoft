////////////////////////////////////////////////////////////////////////
/// \file  IonizationReadout_module.cc
/// \brief Create the readout from the ionization signal
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_IONIZATIONREADOUT
#define GAR_READOUTSIMULATION_IONIZATIONREADOUT

// C++ Includes
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <sys/stat.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/GeometryGAr.h"
#include "CoreUtils/ServiceUtil.h"

// ROOT Includes
#include "TVector3.h"
#include "TFile.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {

    /**
     * @brief Runs Readout simulation including propagation of electrons and photons to readout
     *
     *
     * Randomness
     * -----------
     *
     * The random number generators used by this process are:
     * - 'GEANT' instance: used by Geant4
     * - 'propagation' instance: used in electron propagation
     *
     */
    class IonizationReadout : public ::art::EDProducer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit IonizationReadout(fhicl::ParameterSet const& pset);
      virtual ~IonizationReadout();

      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      void reconfigure(fhicl::ParameterSet const& pset);

    private:

      void DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                   std::vector<edepIDE>                 & edepIDEs);
      void CombineIDEs(std::vector<edepIDE>                 & edepIDEs,
                       std::vector<sdp::EnergyDeposit> const& edepCol);
      void CreateSignalDigit(unsigned int                                      const& channel,
                             std::vector<float>                                     & electrons,
                             std::set<size_t>                                       & eDepLocs,
                             std::deque<float>                                      & eDepWeights,
                             std::vector<raw::RawDigit>                             & digCol,
                             ::art::ValidHandle< std::vector<sdp::EnergyDeposit> >  & eDepCol,
                             ::art::Assns<sdp::EnergyDeposit, raw::RawDigit, float> & erassn,
                             ::art::Event                                           & evt);
      void CheckChannelToEnergyDepositMapping(unsigned int       const& channel,
                                              sdp::EnergyDeposit const& edep,
                                              std::string        const& id);

      std::string                         fG4Label;     ///< label of G4 module
      std::unique_ptr<ElectronDriftAlg>   fDriftAlg;    ///< algorithm to drift ionization electrons
      const gar::detinfo::DetectorClocks* fTime;        ///< electronics clock
      std::unique_ptr<TPCReadoutSimAlg>   fROSimAlg;    ///< algorithm to simulate the electronics
      fhicl::ParameterSet                 fISCalcPars;  ///< parameter set for the IS calculator
      size_t                              fNumTicks;    ///< number of TDC samples
      const gar::geo::GeometryCore*       fGeo;         ///< geometry information
      bool                                fCheckChan;   ///< flag to check mapping of energy deposits to channels
      CLHEP::HepRandomEngine              &fEngine;  ///< random engine
      std::string                         fPRFFileName; ///< where to find the pad response function histograms 

      TH2F                               *fHFILLPRF;   ///< pad response function for hole-filler chamber
      TH2F                               *fIROCPRF;    ///< pad response function for IROC
      TH2F                               *fIOROCPRF;   ///< pad response function for IOROC
      TH2F                               *fOOROCPRF;   ///< pad response function for OOROC

      bool                               fUsePRF;      ///< switch to turn on PRF modeling, otherwise just use the arrival pad
    };

  } // namespace rosim





  namespace rosim {

    //----------------------------------------------------------------------
    // Constructor
    IonizationReadout::IonizationReadout(fhicl::ParameterSet const& pset) : art::EDProducer{pset},
      fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,pset,"Seed")) {

      fTime  = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();

      fNumTicks = gar::providerFrom<detinfo::DetectorPropertiesService>()->NumberTimeSamples();

      fGeo = gar::providerFrom<geo::GeometryGAr>();

      this->reconfigure(pset);

      produces< std::vector<raw::RawDigit>                      >();
      produces< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit, float> >();


      // read in the pad response function histograms

      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(fPRFFileName, fullname);
      struct stat sb;
      if (fullname.empty() || stat(fullname.c_str(), &sb)!=0)
        throw cet::exception("IonizationReadout") << "Input pad response function file "
                          << fPRFFileName
                          << " not found in FW_SEARCH_PATH\n";

      TFile infile(fullname.c_str(),"READ");  // file will close when infile goes out of scope
      fHFILLPRF = (TH2F*) infile.Get("respHFILL");
      fIROCPRF  = (TH2F*) infile.Get("respIROC");
      fIOROCPRF = (TH2F*) infile.Get("respIOROC");
      fOOROCPRF = (TH2F*) infile.Get("respOOROC");
      fHFILLPRF->SetDirectory(0);
      fIROCPRF->SetDirectory(0);
      fIOROCPRF->SetDirectory(0);
      fOOROCPRF->SetDirectory(0);

      return;
    }



    //----------------------------------------------------------------------
    // Destructor
    IonizationReadout::~IonizationReadout() {}



    //----------------------------------------------------------------------
    void IonizationReadout::reconfigure(fhicl::ParameterSet const& pset) {

      MF_LOG_DEBUG("IonizationReadout") << "Debug: IonizationReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;

      fISCalcPars  = pset.get<fhicl::ParameterSet>("ISCalcPars"                 );
      fG4Label     = pset.get<std::string        >("G4ModuleLabel",      "geant");
      fCheckChan   = pset.get<bool               >("CheckChannelMapping", false );
      fPRFFileName = pset.get<std::string        >("PRFFileName",        "MPD/TPCPRF/mpdtpcprf_v1.root");
      fUsePRF      = pset.get<bool               >("UsePRF"             , true  );

      auto driftAlgPars = pset.get<fhicl::ParameterSet>("ElectronDriftAlgPars");
      auto driftAlgName = driftAlgPars.get<std::string>("DriftAlgType");


      if (driftAlgName.compare("Standard") == 0) {
        fDriftAlg = std::make_unique<gar::rosim::ElectronDriftStandardAlg>(fEngine,
                                                                           driftAlgPars);
      } else {
        throw cet::exception("IonizationReadout")
          << "Unable to determine which electron drift algorithm to use, bail";
      }

      auto tpcROAlgPars = pset.get<fhicl::ParameterSet>("TPCReadoutSimAlgPars");
      auto tpcROAlgName = tpcROAlgPars.get<std::string>("TPCReadoutSimType");

      if (tpcROAlgName.compare("Standard") == 0){
        fROSimAlg = std::make_unique<gar::rosim::TPCReadoutSimStandardAlg>(fEngine,
                                                                           tpcROAlgPars);
      } else {
        throw cet::exception("IonizationReadout")
          << "Unable to determine which TPC readout simulation algorithm to use, bail";
      }

      return;
    }



    //----------------------------------------------------------------------
    void IonizationReadout::beginJob() {

      //auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());

      // create the ionization and scintillation calculator;
      // this is a singleton (!) so we just need to make the instance in one
      // location
      IonizationAndScintillation::CreateInstance(fEngine,
                                                 fISCalcPars);

      return;
    }



    //--------------------------------------------------------------------------
    void IonizationReadout::beginRun(::art::Run& /* run */) {
      return;
    }



    //--------------------------------------------------------------------------
    void IonizationReadout::produce(::art::Event& evt) {

      MF_LOG_DEBUG("IonizationReadout") << "produce()";
      //std::cout << "Ionization readout produce" << std::endl;
      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<raw::RawDigit>                      > rdCol (new std::vector<raw::RawDigit>                     );
      std::unique_ptr< ::art::Assns<sdp::EnergyDeposit, raw::RawDigit, float> > erassn(new ::art::Assns<sdp::EnergyDeposit, raw::RawDigit, float>);

      // first get the energy deposits from the event record
      auto eDepCol = evt.getValidHandle< std::vector<sdp::EnergyDeposit> >(fG4Label);

      if(eDepCol->size() > 0){

        std::vector<edepIDE> eDepIDEs;

	// drift the ionization electrons to the readout and create edepIDE objects
  //std::cout << "Drift electrons" << std::endl;
	this->DriftElectronsToReadout(*eDepCol, eDepIDEs);
  //std::cout << "Drifted electrons" << std::endl;
  //std::cout << "eDepIDEs.size() " << eDepIDEs.size() <<std::endl; 
	if (eDepIDEs.size()>0)
	  {
     //std::cout << "eDepIDEs.size()>0" << std::endl;
            // IDEs have been combined already; there are no repeat TDC values for any channel
	    unsigned int       prevChan = eDepIDEs.front().Channel;
	    std::set<size_t>   digitEDepLocs;
            std::deque<float>  digitEDepWeights;
	    std::vector<float> electrons(fNumTicks, 0.);

	    // make the signal raw digits and set their associations to the energy deposits
	    for(auto edide : eDepIDEs){

	      MF_LOG_DEBUG("IonizationReadout")
		<< "Current eDepIDE channel is "
		<< edide.Channel
		<< " previous channel is "
		<< prevChan;

	      if(edide.Channel != prevChan){
		MF_LOG_DEBUG("IonizationReadout")
		  << "There are  "
		  << digitEDepLocs.size()
		  << " locations for "
		  << edide.Channel
		  << " rdCol size is currently "
		  << rdCol->size();

		// this method clears the electrons and digitEDepLocs collections
		// after creating the RawDigit
    //std::cout << "Create signal digit" << std::endl;
		this->CreateSignalDigit(prevChan,
					electrons,
					digitEDepLocs,
                                        digitEDepWeights,
					*rdCol,
					eDepCol,
					*erassn,
					evt);

		// reset the previous channel info
		prevChan = edide.Channel;

	      }

  	      // put overflow times in the last bin.  Is this okay?  TODO
              size_t esize = electrons.size();
              if (esize>0) {
                if (edide.TDC >= esize) {
                  electrons[esize - 1] = edide.NumElect;
                } else {
                  electrons[edide.TDC] = edide.NumElect;
                }
              }

	      for(auto loc : edide.edepLocs) digitEDepLocs.insert(loc);
              for(auto w : edide.edepWeights) digitEDepWeights.push_front(w);

	    } // end loop to fill signal raw digit vector and make EnergyDeposit associations

	    // still one more digit to make because we ran out of channels to compare against
	    this->CreateSignalDigit(eDepIDEs.back().Channel,
				    electrons,
				    digitEDepLocs,
                                    digitEDepWeights,
				    *rdCol,
				    eDepCol,
				    *erassn,
				    evt);

	    MF_LOG_DEBUG("IonizationReadout")
	      << "Created "
	      << rdCol->size()
	      << " raw digits from signal";

	    // now make the noise digits
	    // to do -- only make noise digits on channels we haven't
	    // yet considered for noise digits, but which may have
	    // been entirely zero-suppressed -- may need to keep a
	    // list of channels and pass it in
	    fROSimAlg->CreateNoiseDigits(*rdCol);

        } // end if the EdepIDEs have any size
      } // end if there were energy deposits to use

      evt.put(std::move(rdCol));
      evt.put(std::move(erassn));

      return;
    } // IonizationReadout::produce()

    //--------------------------------------------------------------------------
    void IonizationReadout::DriftElectronsToReadout(std::vector<sdp::EnergyDeposit> const& edepCol,
                                                    std::vector<edepIDE>                 & edepIDEs)
    {
      auto geo = gar::providerFrom<geo::GeometryGAr>();
      //std::cout << "geometry: "  << geo << std::endl;

      float        xyz[3] = {0.};
      unsigned int chan   = 0;

      std::vector<edepIDE> edepIDEaccumulator;

      // now instantiate an ElectronDriftInfo object to keep track of the
      // drifted locations of each electron cluster from each energy deposit
      rosim::ElectronDriftInfo driftInfo;
      //std::cout << "edepCol.size() " << edepCol.size() << std::endl;
      // loop over the energy deposits
      for (size_t e = 0; e < edepCol.size(); ++e) {
        //std::cout << "edepCol e " << e << std::endl;

        // get the positions, arrival times, and electron cluster sizes
        // for this energy deposition
        fDriftAlg->DriftElectronsToReadout(edepCol[e], driftInfo);

        // auto is vector<double> or vector<int> in the clusterSize case.  Size of vector is
        // number of clusters produced by DriftElectronsToReadout.
        auto clusterXPos = driftInfo.ClusterXPos();
        auto clusterYPos = driftInfo.ClusterYPos();
        auto clusterZPos = driftInfo.ClusterZPos();
        auto clusterTime = driftInfo.ClusterTime();
        auto clusterSize = driftInfo.ClusterSize();

        // the vectors should all have the same size by the time we get them
        // here (verified by the ElectronDriftInfo object when they are filled)
        for (size_t c = 0; c < clusterXPos.size(); ++c) {

          xyz[0] = clusterXPos[c];
          xyz[1] = clusterYPos[c];
          xyz[2] = clusterZPos[c];

          // map the cluster's drift point to channels.  Use the method that also gives
          // us a list of neighboring channels.

          gar::geo::ChanWithNeighbors cwn;
          geo->NearestChannelInfo(xyz, cwn);

          // the first channel in the list is the nearest one to the xyz point
          chan = cwn.at(0).id;
          //std::cout << "the first channel in the list is the nearest one to the xyz point " << chan << std::endl; 

          // if charge is deposited on the cover electrodes or is otherwise in a gap, skip it

          if (chan == geo->GapChannelNumber()) continue;

          // incorporate pad response function. 
          // figure out what fraction of the charge of this cluster is to be deposted in each of the channels
          // in chanset.  Need to know if the channel is iroc, ioroc, ooroc, or hole-filler, and use the pad response functions

          TVector3 pos = cwn.at(0).pos;
          TVector3 pproj(0,xyz[1],xyz[2]);  // assume the pad planes are in the YZ plane

          std::vector<float> chanweight;
          size_t ncdistrib = cwn.size();
          if (!fUsePRF) ncdistrib = 1;    //  localize to just one channel if we aren't using the PRF
          float sumw = 0;

          //std::cout << "edepCol e " << e << std::endl;
          //std::cout << "cluster c (x,y,z) " << c << "(" << clusterXPos[c] << ", " << clusterYPos[c] << ", " << clusterZPos[c] << ")" << std::endl;
          //std::cout << "channel (x, y, z) " << chan << "(" << cwn.at(0).pos.X() << ", " << cwn.at(0).pos.Y() << ", " << cwn.at(0).pos.Z() << ")" << std::endl;

          MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << "pproj: " << pproj.Y() << " " << pproj.Z() << std::endl;
          //std::cout<< "ChanWithNeighbors size " << ncdistrib << std::endl;
          for (size_t icd=0; icd<ncdistrib; ++icd) {
            //std::cout<< "ChanWithNeighbors size icd " << icd << std::endl;
            TH2F *prfhist=0;
            if (cwn.at(icd).roctype == gar::geo::HFILLER) {
              prfhist = fHFILLPRF;
              MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << " hole filler" << std::endl;
            } else if (cwn.at(icd).roctype == gar::geo::IROC) {
              prfhist = fIROCPRF;
              MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << " iroc " << std::endl;
            } else if (cwn.at(icd).roctype == gar::geo::IOROC) {
              prfhist = fIOROCPRF;
              MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << " ioroc " << std::endl;
            } else if (cwn.at(icd).roctype == gar::geo::OOROC) {
              prfhist = fOOROCPRF;
              MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << " ooroc " << std::endl;
            } else {
              throw cet::exception("IonizationReadout::DriftElectronsToReadout") << "Ununderstood readout chamber type "
                                << cwn.at(icd).roctype << "\n";
            }
            
            //std::cout<< "Which ROC " << prfhist << std::endl;
            TVector3 dproj = pproj - cwn.at(icd).pos;
            dproj.SetX(0);
            MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << " Pad loc: " << cwn.at(icd).pos.Y() << 
              " " << cwn.at(icd).pos.Z() << std::endl;

            float dist_along_padrow = dproj.Dot(cwn.at(icd).padrowdir);  
            float dist_perp_padrow = (dproj - dist_along_padrow*cwn.at(icd).padrowdir).Mag(); // assume symmetric PRF
            MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout") << "along, perp: " << dist_along_padrow << 
              " " << dist_perp_padrow << std::endl;

            dist_along_padrow = TMath::Abs(dist_along_padrow);
            //std::cout << "dist_along_padrow " << dist_along_padrow << std::endl;
            if (dist_along_padrow < prfhist->GetXaxis()->GetBinUpEdge(prfhist->GetNbinsX()) &&
                dist_perp_padrow < prfhist->GetYaxis()->GetBinUpEdge(prfhist->GetNbinsY())) {
              chanweight.push_back(prfhist->GetBinContent(prfhist->FindBin(dist_along_padrow,dist_perp_padrow)));
            } else {
              chanweight.push_back(0);
    
            }
            sumw += chanweight.back();
          }

          //std::cout << "sumw " << sumw << std::endl; 

          if (sumw == 0) {
            //std::cout<< "Skipped" << std::endl;
            //continue;
            throw cet::exception("IonizationReadout::DriftElectronsToReadout") << 
            "Weight sum is zero, even when including the closest channel " << std::endl;
          }
          float rsumw = 1.0/sumw;
          //std::cout << "rsumw " << rsumw << std::endl; 
          //std::cout << "chan weight size " << chanweight.size() << std::endl; 
          for (size_t i=0; i<chanweight.size(); ++i) {
            //std::cout<< "chan weight i " << i << std::endl;
            if (chanweight.at(i)>0 && clusterSize.at(c) > 0) {
              edepIDEaccumulator.emplace_back(clusterSize.at(c)*chanweight.at(i)*rsumw,
                                              cwn.at(i).id,
                                              fTime->TPCG4Time2TDC(clusterTime.at(c)),
                                              e, chanweight.at(i));
              
              this->CheckChannelToEnergyDepositMapping(edepIDEaccumulator.back().Channel,
                                                       edepCol[e],
                                                       "DriftElectronsToReadout");
	            // compress as we go to save memory
	            if (edepIDEaccumulator.size() > 10000){
                //std::cout << "IonizationReadout::DriftElectronsToReadout 4" << std::endl;
                this->CombineIDEs(edepIDEaccumulator, edepCol);
		            edepIDEs.insert(edepIDEs.end(),edepIDEaccumulator.begin(),edepIDEaccumulator.end());
		            edepIDEaccumulator.clear();
		          }
            }
          }
          
          MF_LOG_DEBUG("IonizationReadout::DriftElectronsToReadout")
            << "cluster time: "
            << clusterTime[c]
            << " TDC "
            << fTime->TPCG4Time2TDC(clusterTime[c])
            << " "
            << fTime->G4ToElecTime(clusterTime[c])
            << " "
            << fTime->TPCClock().TickPeriod();

        }
        //std::cout<< "Hello" << std::endl;
      } // end loop over deposit collections

      //std::cout << "IonizationReadout::DriftElectronsToReadout 5" << std::endl;

      // one last collection of accumulated edepIDEs 
      edepIDEs.insert(edepIDEs.end(),edepIDEaccumulator.begin(),edepIDEaccumulator.end());
      //std::cout<< "edepIDEs size before Combine IDEs " << edepIDEs.size() << std::endl;
      this->CombineIDEs(edepIDEs, edepCol);

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CombineIDEs(std::vector<edepIDE>                 & edepIDEs,
                                        std::vector<sdp::EnergyDeposit> const& edepCol) {

      MF_LOG_DEBUG("IonizationReadout::CombineIDEs")
        << "starting with "
        << edepIDEs.size()
        << " energy deposits";

      if (edepIDEs.size()==0) return;

      std::vector<edepIDE> temp;

      // sort the edepIDE objects.  This is sorting by channel 1st & TDC
      // 2nd via the < operator of edepIDE
      std::sort(edepIDEs.begin(), edepIDEs.end());

      for(auto itr : edepIDEs){
        for(auto edloc : itr.edepLocs) {
          this->CheckChannelToEnergyDepositMapping(itr.Channel,
        										   edepCol[edloc],
        										   "CombineIDEsAfterSort");
        }
      }

      edepIDE prev = edepIDEs.front();
      edepIDE sum(edepIDEs.front());
      edepIDE cur(edepIDEs.front());

      // now loop over the sorted vector and combine any edepIDEs
      // with the same channel and tdc values for the IDEs
      // start with entry 1 as sum is already holding the information
      // from entry 0
      for (size_t e = 1; e < edepIDEs.size(); ++e) {
        cur = edepIDEs[e];

        MF_LOG_DEBUG("IonizationReadout::CombineIDEs")
          << "current edepIDE: "
          << cur.NumElect
          << " "
          << cur.Channel
          << " "
          << cur.TDC
          << " "
          << cur.edepLocs.size();

        if (cur != prev) {
          MF_LOG_DEBUG("IonizationReadout::CombineIDEs")
            << "storing edepIDE sum: "
            << sum.NumElect
            << " "
            << sum.Channel
            << " "
            << sum.TDC
            << " "
            << sum.edepLocs.size();

          if (fCheckChan) {
        	for (auto edloc : sum.edepLocs) {
        	  this->CheckChannelToEnergyDepositMapping(sum.Channel,
        											   edepCol[edloc],
        											   "CombineIDEsStore");
        	}
          }

          // put the summed edepIDE into the temp vector
          temp.push_back(sum);

          // start over with a fresh sum
          sum  = cur;
          prev = cur;

        } else {
          MF_LOG_DEBUG("IonizationReadout::CombineIDEs")
        	<< "summing current edepIDE";
          sum  += cur;
          prev  = cur;
        }

      } // end loop to sum edepIDEs

      // now swap the input vector with the temp vector
      temp.swap(edepIDEs);

      MF_LOG_DEBUG("IonizationReadout::CombineIDEs")
        << "ending with "
        << edepIDEs.size()
        << " energy deposits";

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CreateSignalDigit(unsigned int                                      const& channel,
                                              std::vector<float>                                     & electrons,
                                              std::set<size_t>                                       & eDepLocs,
                                              std::deque<float>                                      & eDepWeights,
                                              std::vector<raw::RawDigit>                             & digCol,
                                              ::art::ValidHandle< std::vector<sdp::EnergyDeposit> >  & eDepCol,
                                              ::art::Assns<sdp::EnergyDeposit, raw::RawDigit, float> & erassn,
                                              ::art::Event                                           & evt)
    {

      // could be that all the adc's fell below threshold, so test if we got a raw digit at all.

      bool todrop=false;
      raw::RawDigit tmpdigit = fROSimAlg->CreateRawDigit(channel, electrons, todrop);

      if (!todrop) {
        digCol.emplace_back(tmpdigit);

        MF_LOG_DEBUG("IonizationReadout::CreateSignalDigit")
          << "Associating "
          << eDepLocs.size()
          << " energy deposits to digit for channel "
          << channel;

        // loop over the locations in the eDepCol to make the associations
        size_t index = 0;
        for (auto ed : eDepLocs) {
          auto const ptr = art::Ptr<sdp::EnergyDeposit>(eDepCol, ed);

          this->CheckChannelToEnergyDepositMapping(channel, *ptr, "CreateSignalDigit");
          MF_LOG_DEBUG("IonizationReadout::CreateSignalDigit")
            << "Making association: " << digCol.size() << " " << ptr << std::endl;
	  util::CreateAssnD(*this, evt, digCol, ptr, eDepWeights.at(index), erassn);
          index++;
        }
      }

      eDepLocs.clear();
      eDepWeights.clear();
      electrons.clear();
      electrons.resize(fNumTicks, 0);

      return;
    }

    //--------------------------------------------------------------------------
    void IonizationReadout::CheckChannelToEnergyDepositMapping(unsigned int       const& channel,
                                                               sdp::EnergyDeposit const& edep,
                                                               std::string        const& id)
    {
      if (!fCheckChan) return;

      // check that the channel for this cluster is close to what we expect
      // for the energy deposit

      float xyz[3] = {0.};
      fGeo->ChannelToPosition(channel, xyz);
      //std::cout<< "Channel to position " << std::endl;

      if(std::abs(edep.Y() - xyz[1]) > 1 ||
         std::abs(edep.Z() - xyz[2]) > 1){

        MF_LOG_DEBUG("IonizationReadout::CheckChannelToEnergyDepositMapping")
          << "In function "
          << id
          << ": Channel "
          << channel
          << " is off from the energy deposit: ("
          << xyz[1]
          << ", "
          << xyz[2]
          << ") vs ("
          << edep.Y()
          << ", "
          << edep.Z()
          << ") " << std::endl;
      }

      return;
    }

  } // namespace rosim

  namespace rosim {

    DEFINE_ART_MODULE(IonizationReadout)

  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_IONIZATIONREADOUT
