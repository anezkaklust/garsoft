#ifndef CONVEDEP2ART_H
#define CONVEDEP2ART_H

// C++ Includes
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <iostream>
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
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/G4Base/DetectorConstruction.h"
#include "nutools/G4Base/UserActionManager.h"
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// ROOT Includes
#include "TGeoManager.h"

///Geant4 interface
namespace gar {
    namespace garg4 {

        class ConvertEdep2Art : public ::art::EDProducer{
        public:

            /// Standard constructor and destructor for an FMWK module.
            explicit ConvertEdep2Art(fhicl::ParameterSet const& pset);
            virtual ~ConvertEdep2Art();

            /// The main routine of this module: Fetch the primary particles
            /// from the event, simulate their evolution in the detctor, and
            /// produce the detector response.
            void produce (::art::Event& evt);
            void beginJob();
            void beginRun(::art::Run& run);

        private:


        };

    } // namespace garg4

    namespace garg4 {

        //----------------------------------------------------------------------
        // Constructor
        ConvertEdep2Art::ConvertEdep2Art(fhicl::ParameterSet const& pset)
        : art::EDProducer{pset}
        {


        }

        //----------------------------------------------------------------------
        // Destructor
        ConvertEdep2Art::~ConvertEdep2Art()
        {

        }

        //----------------------------------------------------------------------
        void ConvertEdep2Art::beginJob()
        {

            return;
        }

        //--------------------------------------------------------------------------
        void ConvertEdep2Art::beginRun(::art::Run& run)
        {

            return;
        }

        //--------------------------------------------------------------------------
        void ConvertEdep2Art::produce(::art::Event& evt)
        {

            return;
        }

    } // namespace garg4

    namespace garg4 {

        DEFINE_ART_MODULE(ConvertEdep2Art)

    } // namespace garg4
} // gar

#endif
