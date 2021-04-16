//
//  G4Particle.h
//
//  Created by C. Hilgenberg
//
//  Inital and final 4-vectors with parantage info
//
//  Extension of garana::G4Particle to facilitate
//   construction using nusimdata::MCParticle

#ifndef GAR_ANALYSISDATAPRODUCTS_G4Particle_h
#define GAR_ANALYSISDATAPRODUCTS_G4Particle_h

#include "nusimdata/SimulationBase/MCParticle.h"
#include "AnalysisDataProducts/AnaUtils.h"
#include "garana/DataProducts/G4Particle.h"

#include <TLorentzVector.h>

#include <vector>
#include <utility>

using std::vector;
using std::pair;

namespace gar {

   class GarG4Particle{

    public:

       GarG4Particle(){}
       GarG4Particle(const simb::MCParticle& mcp, int parentPdg, int progenitorPdg, int progenitorTrackId,
                     const vector<pair<TLorentzVector,TLorentzVector>>& positions,
    		     const vector<pair<TLorentzVector,TLorentzVector>>& momenta, const vector<int>& regions,
                     const vector<size_t>& nptsPerRegion );


       garana::G4Particle GetGaranaObj() { return *fG4P; };

    private:

           garana::G4Particle* fG4P = nullptr;

   };//class
}//namespace gar

#endif
