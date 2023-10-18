////////////////////////////////////////////////////////////////////////
// $Id: GENIEGen_plugin.cc,v 1.4 2010/04/27 19:48:46 brebel Exp $
//
//
// GENIE neutrino event generator
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_GENIEGEN_H
#define EVGEN_GENIEGEN_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unistd.h>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TStopwatch.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Core/EDProducer.h"

#include "GENIE/Framework/Ntuple/NtpMCEventRecord.h"
#include "GENIE/Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nugen/EventGeneratorBase/GENIE/GENIEHelper.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// GArSoft includes
#include "Geometry/GeometryGAr.h"
#include "SummaryDataProducts/RunData.h"
#include "SummaryDataProducts/POTSummary.h"
#include "SimulationDataProducts/GenieParticle.h"
#include "SimulationDataProducts/BeamTypes.h"
#include "SimulationDataProducts/BeamGateInfo.h"
#include "SimulationDataProducts/sim.h"
#include "Utilities/AssociationUtil.h"
#include "CoreUtils/ServiceUtil.h"

///Event Generation using GENIE, cosmics or single particles
namespace gar {
    namespace evgen {
        /**
        * @brief A module to check the results from the Monte Carlo generator
        *
        * Note on random number generator
        * --------------------------------
        *
        * GENIE uses a TRandom generator for its purposes.
        * Since art's RandomNumberGenerator service only provides
        * `CLHEP::HepRandomEngine`, the standard GArSoft/art mechanism for handling
        * the random stream can't be used.
        * GENIEHelper, interface to GENIE provided by nutools, creates a TRandom
        * that GENIE can use. It initializes it with a random seed read from
        * *RandomSeed* configuration parameter. This and all the other parameters
        * are inherited from the art module (that is, `GENIEGen`) configuration.
        * GArSoft meddles with this mechanism to provide support for the standard
        * "Seed" parameter and NuRandomService service.
        *
        * Configuration parameters
        * -------------------------
        *
        * - *RandomSeed* (integer, optional): if specified, this value is used as
        *   seed for GENIE random number generator engine
        * - *Seed* (unsigned integer, optional): if specified, this value is used as
        *   seed for GENIE random number generator engine; if *RandomSeed* is also
        *   specified, this value is ignored (but in the future this could turn into
        *   a configuration error)
        *
        * As custom, if the random seed is not provided by the configuration, one is
        * fetched from `NuRandomService` (if available), with the behaviour in
        * gar::rndm::FetchRandomSeed().
        */
        class GENIEGen : public ::art::EDProducer {
        public:
            explicit GENIEGen(fhicl::ParameterSet const &pset);
            virtual ~GENIEGen();

            void produce(::art::Event& evt);
            void beginJob();
            void beginRun(::art::Run& run);
            void endRun(::art::Run& run);

        private:

            std::string ParticleStatus(int StatusCode);
            std::string ReactionChannel(int ccnc,int mode);

            void FillHistograms(simb::MCTruth mc);

            evgb::GENIEHelper  *fGENIEHelp;       ///< GENIEHelper object
            int                 fPassEmptySpills; ///< whether or not to kill evnets with no interactions
            TStopwatch          fStopwatch;       ///keep track of how long it takes to run the job

            double fGlobalTimeOffset;             /// The start of a simulated "beam gate".
            double fRandomTimeOffset;             /// The width of a simulated "beam gate".
            gar::sdp::BeamType_t fBeamType;       /// The type of beam

            TH1F* fGenerated[6];  ///< Spectra as generated

            TH1F* fVertexX;    ///< vertex location of generated events in x
            TH1F* fVertexY;    ///< vertex location of generated events in y
            TH1F* fVertexZ;    ///< vertex location of generated events in z

            TH2F* fVertexXY;   ///< vertex location in xy
            TH2F* fVertexXZ;   ///< vertex location in xz
            TH2F* fVertexYZ;   ///< vertex location in yz

            TH1F* fDCosX;      ///< direction cosine in x
            TH1F* fDCosY;      ///< direction cosine in y
            TH1F* fDCosZ;      ///< direction cosine in z

            TH1F* fMuMomentum; ///< momentum of outgoing muons
            TH1F* fMuDCosX;    ///< direction cosine of outgoing mu in x
            TH1F* fMuDCosY;    ///< direction cosine of outgoing mu in y
            TH1F* fMuDCosZ;    ///< direction cosine of outgoing mu in z

            TH1F* fEMomentum;  ///< momentum of outgoing electrons
            TH1F* fEDCosX;     ///< direction cosine of outgoing e in x
            TH1F* fEDCosY;     ///< direction cosine of outgoing e in y
            TH1F* fEDCosZ;     ///< direction cosine of outgoing e in z

            TH1F* fCCMode;      ///< CC interaction mode
            TH1F* fNCMode;      ///< CC interaction mode

            TH1F* fDeltaE;     ///< difference in neutrino energy from MCTruth::Enu() vs TParticle
            TH1F* fECons;      ///< histogram to determine if energy is conserved in the event

        };
    }

    namespace evgen{

        //____________________________________________________________________________
        GENIEGen::GENIEGen(fhicl::ParameterSet const& pset)
        : EDProducer{pset},
        fGENIEHelp(0),
        fPassEmptySpills (pset.get< bool   >("PassEmptySpills")),
        fGlobalTimeOffset(pset.get< double >("GlobalTimeOffset",0)),
        fRandomTimeOffset(pset.get< double >("RandomTimeOffset",1600.)), // BNB default value
        fBeamType(gar::sdp::kBNB)
        {
            fStopwatch.Start();

            produces< std::vector<simb::MCTruth> >();
            produces< std::vector<simb::MCFlux>  >();
            produces< std::vector<simb::GTruth>  >();
            produces< std::vector<sdp::GenieParticle> >();
            produces< sumdata::RunData, ::art::InRun >();
            produces< sumdata::POTSummary, ::art::InRun >();
            produces< ::art::Assns<simb::MCTruth, simb::MCFlux> >();
            produces< ::art::Assns<simb::MCTruth, simb::GTruth> >();

            std::string beam_type_name = pset.get<std::string>("BeamName");

            if     (beam_type_name == "numi"   ) fBeamType = gar::sdp::kNuMI;
            else if(beam_type_name == "booster") fBeamType = gar::sdp::kBNB;
            else                                 fBeamType = gar::sdp::kUnknown;

            auto geo = gar::providerFrom<geo::GeometryGAr>();

            fhicl::ParameterSet GENIEconfig(pset);
            if (!GENIEconfig.has_key("RandomSeed")) {
                // no RandomSeed specified; check for the GArSoft-style "Seed" instead:
                // obtain the random seed from a service,
                // unless overridden in configuration with key "Seed"

                unsigned int seed;
                if (!GENIEconfig.has_key("Seed"))
                    seed = ::art::ServiceHandle<rndm::NuRandomService>()->getSeed();

                // The seed is not passed to RandomNumberGenerator,
                // since GENIE uses a TRandom generator that is owned by the GENIEHelper.
                // Instead, we explicitly configure the random seed for GENIEHelper:
                GENIEconfig.put("RandomSeed", seed);
            } // if no RandomSeed present


            fGENIEHelp = new evgb::GENIEHelper(GENIEconfig,
                                               geo->ROOTGeoManager(),
                                               geo->ROOTFile(),
                                               geo->TotalMass(pset.get< std::string>("TopVolume").c_str()));


        }

        //____________________________________________________________________________
        GENIEGen::~GENIEGen()
        {
            // if(fGENIEHelp) delete fGENIEHelp;
            fStopwatch.Stop();
            mf::LogInfo("GENIEProductionTime") << "real time to produce file: " << fStopwatch.RealTime();
        }

        //____________________________________________________________________________
        void GENIEGen::beginJob(){
            fGENIEHelp->Initialize();

            // Get access to the TFile service.
            ::art::ServiceHandle<::art::TFileService> tfs;

            fGenerated[0] = tfs->make<TH1F>("fGenerated_necc",  "", 100, 0.0, 20.0);
            fGenerated[1] = tfs->make<TH1F>("fGenerated_nebcc", "", 100, 0.0, 20.0);
            fGenerated[2] = tfs->make<TH1F>("fGenerated_nmcc",  "", 100, 0.0, 20.0);
            fGenerated[3] = tfs->make<TH1F>("fGenerated_nmbcc", "", 100, 0.0, 20.0);
            fGenerated[4] = tfs->make<TH1F>("fGenerated_nnc",   "", 100, 0.0, 20.0);
            fGenerated[5] = tfs->make<TH1F>("fGenerated_nbnc",  "", 100, 0.0, 20.0);

            fDCosX = tfs->make<TH1F>("fDCosX", ";dx/ds", 200, -1., 1.);
            fDCosY = tfs->make<TH1F>("fDCosY", ";dy/ds", 200, -1., 1.);
            fDCosZ = tfs->make<TH1F>("fDCosZ", ";dz/ds", 200, -1., 1.);

            fMuMomentum = tfs->make<TH1F>("fMuMomentum", ";p_{#mu} (GeV/c)", 500, 0., 50.);
            fMuDCosX    = tfs->make<TH1F>("fMuDCosX", ";dx/ds;", 200, -1., 1.);
            fMuDCosY    = tfs->make<TH1F>("fMuDCosY", ";dy/ds;", 200, -1., 1.);
            fMuDCosZ    = tfs->make<TH1F>("fMuDCosZ", ";dz/ds;", 200, -1., 1.);

            fEMomentum  = tfs->make<TH1F>("fEMomentum", ";p_{e} (GeV/c)", 500, 0., 50.);
            fEDCosX     = tfs->make<TH1F>("fEDCosX", ";dx/ds;", 200, -1., 1.);
            fEDCosY     = tfs->make<TH1F>("fEDCosY", ";dy/ds;", 200, -1., 1.);
            fEDCosZ     = tfs->make<TH1F>("fEDCosZ", ";dz/ds;", 200, -1., 1.);

            fCCMode = tfs->make<TH1F>("fCCMode", ";CC Interaction Mode;", 4, 0., 4.);
            fCCMode->GetXaxis()->SetBinLabel(1, "QE");
            fCCMode->GetXaxis()->SetBinLabel(2, "Res");
            fCCMode->GetXaxis()->SetBinLabel(3, "DIS");
            fCCMode->GetXaxis()->SetBinLabel(4, "Coh");
            fCCMode->GetXaxis()->CenterLabels();

            fNCMode = tfs->make<TH1F>("fNCMode", ";NC Interaction Mode;", 4, 0., 4.);
            fNCMode->GetXaxis()->SetBinLabel(1, "QE");
            fNCMode->GetXaxis()->SetBinLabel(2, "Res");
            fNCMode->GetXaxis()->SetBinLabel(3, "DIS");
            fNCMode->GetXaxis()->SetBinLabel(4, "Coh");
            fNCMode->GetXaxis()->CenterLabels();

            fDeltaE = tfs->make<TH1F>("fDeltaE", ";#Delta E_{#nu} (GeV);", 200, -1., 1.);
            fECons  = tfs->make<TH1F>("fECons", ";#Delta E(#nu,lepton);", 500, -5., 5.);

            auto   geo  = gar::providerFrom<geo::GeometryGAr>();
            double x    = 2.1 * geo->TPCRadius();
            double y    = 2.1 * geo->TPCRadius();
            double z    = 2.  * geo->TPCLength();
            int    xdiv = TMath::Nint(2 * x / 5.);
            int    ydiv = TMath::Nint(2 * y / 5.);
            int    zdiv = TMath::Nint(2 * z / 5.);

            fVertexX = tfs->make<TH1F>("fVertexX", ";x (cm)", xdiv, -x, x);
            fVertexY = tfs->make<TH1F>("fVertexY", ";y (cm)", ydiv, -y, y);
            fVertexZ = tfs->make<TH1F>("fVertexZ", ";z (cm)", zdiv, -z, z);

            fVertexXY = tfs->make<TH2F>("fVertexXY", ";x (cm);y (cm)", xdiv, -x, x, ydiv, -y, y);
            fVertexXZ = tfs->make<TH2F>("fVertexXZ", ";z (cm);x (cm)", zdiv, -z, z, xdiv, -x, x);
            fVertexYZ = tfs->make<TH2F>("fVertexYZ", ";z (cm);y (cm)", zdiv, -z, z, ydiv, -y, y);

        }

        //____________________________________________________________________________
        void GENIEGen::beginRun(::art::Run& run)
        {
            // grab the geometry object to see what geometry we are using
            auto geo = gar::providerFrom<geo::GeometryGAr>();
            std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
            run.put(std::move(runcol),art::fullRun());
            return;
        }

        //____________________________________________________________________________
        void GENIEGen::endRun(::art::Run& run)
        {
            std::unique_ptr<sumdata::POTSummary> p(new sumdata::POTSummary());
            p->totpot     = fGENIEHelp->TotalExposure();
            p->totgoodpot = fGENIEHelp->TotalExposure();
            run.put(std::move(p), art::fullRun());
            return;
        }

        //____________________________________________________________________________
        void GENIEGen::produce(::art::Event& evt)
        {
            std::unique_ptr< std::vector<simb::MCTruth> >                truthcol (new std::vector<simb::MCTruth>);
            std::unique_ptr< std::vector<simb::MCFlux>  >                fluxcol  (new std::vector<simb::MCFlux >);
            std::unique_ptr< std::vector<simb::GTruth>  >                gtruthcol(new std::vector<simb::GTruth >);
            std::unique_ptr< ::art::Assns<simb::MCTruth, simb::MCFlux> > tfassn   (new ::art::Assns<simb::MCTruth, simb::MCFlux>);
            std::unique_ptr< ::art::Assns<simb::MCTruth, simb::GTruth> > tgtassn  (new ::art::Assns<simb::MCTruth, simb::GTruth>);
            std::unique_ptr< std::vector<gar::sdp::GenieParticle> > geniepartcol( new std::vector<gar::sdp::GenieParticle> );

            while(truthcol->size() < 1)
            {
                unsigned int_idx = 0;
                while(!fGENIEHelp->Stop())
                {
                    simb::MCTruth truth;        // These variables - well, truth in particular - need to have their constructor
                    simb::MCFlux  flux;         // called directly prior to the invocation of GENIEHelper::Sample.  It's the
                    simb::GTruth  gTruth;       // fNeutrinoSet field of MCTruth that's a problem, if not others as well

                    // GENIEHelper returns a false in the sample method if
                    // either no neutrino was generated, or the interaction
                    // occurred beyond the detector's z extent - ie something we
                    // would never see anyway.
                    if(fGENIEHelp->Sample(truth, flux, gTruth)) {

                        //Store the list of GENIE particles
                        genie::EventRecord *event = fGENIEHelp->GetGenieEventRecord();
                        genie::GHepParticle* p = nullptr;
                        TObjArrayIter piter(event);
                        unsigned int idx = 0;

                        while( (p = (genie::GHepParticle*) piter.Next()) ) {
                            geniepartcol->emplace_back(int_idx, idx, p->Pdg(), p->Status(), p->Name(), p->FirstMother(), p->LastMother(), p->FirstDaughter(), p->LastDaughter(), p->Px(), p->Py(), p->Pz(), p->E(), p->Mass(), p->RescatterCode());
                            idx++;
                        }

                        truthcol ->push_back(truth);
                        fluxcol  ->push_back(flux);
                        gtruthcol->push_back(gTruth);
                        util::CreateAssn(*this, evt, *truthcol, *fluxcol,   *tfassn,  fluxcol->size()   - 1, fluxcol->size()  );
                        util::CreateAssn(*this, evt, *truthcol, *gtruthcol, *tgtassn, gtruthcol->size() - 1, gtruthcol->size());

                        this->FillHistograms(truth);
                    }// end if genie was able to make an event
                    int_idx++;
                }// end event generation loop

                // check to see if we are to pass empty spills
                if(truthcol->size() < 1 && fPassEmptySpills) {
                    MF_LOG_DEBUG("GENIEGen") << "no events made for this spill but "
                    << "passing it on and ending the event anyway";
                    break;
                }

            }// end loop while no interactions are made

            // put the collections in the event
            evt.put(std::move(truthcol));
            evt.put(std::move(fluxcol));
            evt.put(std::move(gtruthcol));
            evt.put(std::move(tfassn));
            evt.put(std::move(tgtassn));
            evt.put(std::move(geniepartcol));

            return;
        }

        //......................................................................
        std::string GENIEGen::ParticleStatus(int StatusCode)
        {
            int code = StatusCode;
            std::string ParticleStatusName;

            switch(code)
            {
            case -1:
                ParticleStatusName = "kIStUndefined";
                break;
            case 0:
                ParticleStatusName = "kIStInitialState";
                break;
            case 1:
                ParticleStatusName = "kIStStableFinalState";
                break;
            case 2:
                ParticleStatusName = "kIStIntermediateState";
                break;
            case 3:
                ParticleStatusName = "kIStDecayedState";
                break;
            case 11:
                ParticleStatusName = "kIStNucleonTarget";
                break;
            case 12:
                ParticleStatusName = "kIStDISPreFragmHadronicState";
                break;
            case 13:
                ParticleStatusName = "kIStPreDecayResonantState";
                break;
            case 14:
                ParticleStatusName = "kIStHadronInTheNucleus";
                break;
            case 15:
                ParticleStatusName = "kIStFinalStateNuclearRemnant";
                break;
            case 16:
                ParticleStatusName = "kIStNucleonClusterTarget";
                break;
            default:
                ParticleStatusName = "Status Unknown";
            }
            return ParticleStatusName;
        }

        //......................................................................
        std::string GENIEGen::ReactionChannel(int ccnc,int mode)
        {
            std::string ReactionChannelName=" ";

            if(ccnc==0)
            ReactionChannelName = "kCC";
            else if(ccnc==1)
            ReactionChannelName = "kNC";
            else std::cout<<"Current mode unknown!! "<<std::endl;

            if(mode==0)
            ReactionChannelName += "_kQE";
            else if(mode==1)
            ReactionChannelName += "_kRes";
            else if(mode==2)
            ReactionChannelName += "_kDIS";
            else if(mode==3)
            ReactionChannelName += "_kCoh";
            else std::cout<<"interaction mode unknown!! "<<std::endl;

            return ReactionChannelName;
        }

        //......................................................................
        void GENIEGen::FillHistograms(simb::MCTruth mc)
        {
            // Decide which histograms to put the spectrum in
            int id = -1;
            if (mc.GetNeutrino().CCNC()==simb::kCC) {
                fCCMode->Fill(mc.GetNeutrino().Mode());
                if      (mc.GetNeutrino().Nu().PdgCode() ==  12) id = 0;
                else if (mc.GetNeutrino().Nu().PdgCode() == -12) id = 1;
                else if (mc.GetNeutrino().Nu().PdgCode() ==  14) id = 2;
                else if (mc.GetNeutrino().Nu().PdgCode() == -14) id = 3;
                else return;
            }
            else {
                fNCMode->Fill(mc.GetNeutrino().Mode());
                if (mc.GetNeutrino().Nu().PdgCode() > 0) id = 4;
                else                                     id = 5;
            }
            if (id==-1) abort();

            // Fill the specta histograms
            fGenerated[id]->Fill(mc.GetNeutrino().Nu().E() );

            ///< fill the vertex histograms from the neutrino - that is always
            ///< particle 0 in the list
            simb::MCNeutrino       mcnu = mc.GetNeutrino();
            const simb::MCParticle nu   = mcnu.Nu();

            fVertexX->Fill(nu.Vx());
            fVertexY->Fill(nu.Vy());
            fVertexZ->Fill(nu.Vz());

            fVertexXY->Fill(nu.Vx(), nu.Vy());
            fVertexXZ->Fill(nu.Vz(), nu.Vx());
            fVertexYZ->Fill(nu.Vz(), nu.Vy());

            double mom = nu.P();
            if(std::abs(mom) > 0.){
                fDCosX->Fill(nu.Px()/mom);
                fDCosY->Fill(nu.Py()/mom);
                fDCosZ->Fill(nu.Pz()/mom);
            }


            MF_LOG_DEBUG("GENIEInteractionInformation")
            << std::endl
            << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode())
            << std::endl
            << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl
            << std::setiosflags(std::ios::left)
            << std::setw(20) << "PARTICLE"
            << std::setiosflags(std::ios::left)
            << std::setw(32) << "STATUS"
            << std::setw(18) << "E (GeV)"
            << std::setw(18) << "m (GeV/c2)"
            << std::setw(18) << "Ek (GeV)"
            << std::endl << std::endl;

            const TDatabasePDG* databasePDG = TDatabasePDG::Instance();

            // Loop over the particle stack for this event
            for(int i = 0; i < mc.NParticles(); ++i){
                simb::MCParticle part(mc.GetParticle(i));
                std::string name = databasePDG->GetParticle(part.PdgCode())->GetName();
                int code = part.StatusCode();
                std::string status = ParticleStatus(code);
                double mass = part.Mass();
                double energy = part.E();
                double Ek = (energy-mass); // Kinetic Energy (GeV)
                if(status == "kIStStableFinalState" ||
                status == "kIStHadronInTheNucleus")
                MF_LOG_DEBUG("GENIEFinalState")
                << std::setiosflags(std::ios::left) << std::setw(20) << name
                << std::setiosflags(std::ios::left) << std::setw(32) <<status
                << std::setw(18)<< energy
                << std::setw(18)<< mass
                << std::setw(18)<< Ek <<std::endl;
                else
                MF_LOG_DEBUG("GENIEFinalState")
                << std::setiosflags(std::ios::left) << std::setw(20) << name
                << std::setiosflags(std::ios::left) << std::setw(32) << status
                << std::setw(18) << energy
                << std::setw(18) << mass <<std::endl;
            }


            if(mc.GetNeutrino().CCNC() == simb::kCC){

                ///look for the outgoing lepton in the particle stack
                ///just interested in the first one
                for(int i = 0; i < mc.NParticles(); ++i){
                    simb::MCParticle part(mc.GetParticle(i));
                    if(std::abs(part.PdgCode()) == 11){
                        fEMomentum->Fill(part.P());
                        fEDCosX->Fill(part.Px()/part.P());
                        fEDCosY->Fill(part.Py()/part.P());
                        fEDCosZ->Fill(part.Pz()/part.P());
                        fECons->Fill(nu.E() - part.E());
                        break;
                    }
                    else if(std::abs(part.PdgCode()) == 13){
                        fMuMomentum->Fill(part.P());
                        fMuDCosX->Fill(part.Px()/part.P());
                        fMuDCosY->Fill(part.Py()/part.P());
                        fMuDCosZ->Fill(part.Pz()/part.P());
                        fECons->Fill(nu.E() - part.E());
                        break;
                    }
                }// end loop over particles
            }//end if CC interaction

            return;
        }

    }

    namespace evgen{

        DEFINE_ART_MODULE(GENIEGen)

    }
} // end gar
#endif
