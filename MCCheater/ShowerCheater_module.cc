////////////////////////////////////////////////////////////////////////
/// \file  ShowerCheater_module.cc
/// \brief Create rec::Shower products using MC truth information
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_MCCHEATER_SHOWERCHEATER
#define GAR_MCCHEATER_SHOWERCHEATER

  // C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

  // Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "cetlib/exception.h"

  // GArSoft Includes
#include "MCCheater/BackShowerer.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "RawDataProducts/RawDigit.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Shower.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

  // Forward declarations

namespace gar {
  namespace cheat {
    
    class ShowerCheater : public ::art::EDProducer{
    public:
      
        // Standard constructor and destructor for an FMWK module.
      explicit ShowerCheater(fhicl::ParameterSet const& pset);
      virtual ~ShowerCheater();
      
      void produce (::art::Event& evt);
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      std::string fHitLabel; ///< label of module creating raw digits
      std::string fG4Label;  ///< label of module creating mc particles
      
    };
    
  } // cheat
  
  namespace cheat {
    
    //--------------------------------------------------------------------------
    ShowerCheater::ShowerCheater(fhicl::ParameterSet const& pset)
    {
      this->reconfigure(pset);
      
      produces< std::vector<rec::Shower>                    >();
      produces< ::art::Assns<rec::Hit,         rec::Shower> >();
      produces< ::art::Assns<simb::MCParticle, rec::Shower> >();
      
      return;
    }
    
    //--------------------------------------------------------------------------
    ShowerCheater::~ShowerCheater()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    void ShowerCheater::reconfigure(fhicl::ParameterSet const& pset)
    {
      fHitLabel = pset.get<std::string >("HitLabel",      "hit"  );
      fG4Label  = pset.get<std::string >("G4ModuleLabel", "geant");
      return;
    }
    
    //--------------------------------------------------------------------------
    void ShowerCheater::produce(::art::Event & evt)
    {
      if( evt.isRealData() ){
        throw cet::exception("ShowerCheater")
        << "Attempting to cheat the shower reconstruction for real data, "
        << "that will never work";
        return;
      }
      
      std::unique_ptr< std::vector<rec::Shower>                    > shwCol     (new std::vector<rec::Shower>                   );
      std::unique_ptr< ::art::Assns<rec::Hit,         rec::Shower> > hitShwAssns(new ::art::Assns<rec::Hit,         rec::Shower>);
      std::unique_ptr< ::art::Assns<simb::MCParticle, rec::Shower> > pShwAssns  (new ::art::Assns<simb::MCParticle, rec::Shower>);
      
      // get the MCParticles from the event and make a vector of Ptrs of them
      auto partCol = evt.getValidHandle< std::vector<simb::MCParticle> >(fG4Label);
      std::vector<art::Ptr<simb::MCParticle> > partVec;
      art::fill_ptr_vector(partVec, partCol);
      
      // get the FindMany for the digit to EnergyDeposit association
      ::art::FindManyP<rec::Hit> fmphit(partVec, evt, fHitLabel);
      
      // test that the FindMany is valid - don't worry about the digCol, the
      // getValidHandle throws if it can't find a valid handle
      if(!fmphit.isValid() ){
        throw cet::exception("ShowerCheater")
        << "Unable to find valid FindMany<Hit> "
        << fmhit.isValid()
        << " this is a problem for cheating";
      }
      
      std::vector< art::Ptr<rec::Hit> > hits;
      simb::MCParticle const& part;
      float                   energy    = 0.;
      float                   vtx[3]    = {0.};
      float                   vtxDir[3] = {0.};
      
      for(size_t p = 0; p < partVec.size(); ++p){
        fmphit.get(p, hits);
        if(hits.size() < 1) continue;
        
        part = *partVec[p];
        
        // ignore if we don't have an EM shower
        if(std::abs(part->PdgCode()) != 11 ||
           std::abs(part->PdgCode()) != 111) continue;
        
        momentum  = part->E();
        vtx[0]    = part->Vx();
        vtx[1]    = part->Vx();
        vtx[2]    = part->Vx();
        vtxDir[0] = part->Px() / momentum;
        vtxDir[1] = part->Py() / momentum;
        vtxDir[2] = part->Pz() / momentum;
        
        shwCol->emplace_back(energy,
                             vtx,
                             vtxDir);
        
        // make the hit to MCParticle association
        util::CreateAssn(*this, evt, *shwCol, partVec[p], *pShwAssns);
        
        // make the hit to track associations
        for(auto hit : hits) util::CreateAssn(*this, evt, *shwCol, hit, *hitShwAssns);
        
      } // end loop over particles
      
      evt.put(std::move(shwCol     ));
      evt.put(std::move(hitShwAssns));
      evt.put(std::move(pShwAssns  ));
      
      return;
    }
    
    
  } // cheat
  
  namespace cheat {
    
    DEFINE_ART_MODULE(ShowerCheater)
    
  } // cheat
} // gar
#endif // GAR_MCCHEATER_SHOWERCHEATER