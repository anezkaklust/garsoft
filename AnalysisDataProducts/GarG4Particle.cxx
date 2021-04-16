//
//  GarG4Particle.cxx
//
//  Written by C. Hilgenberg

#include "garsoft/AnalysisDataProducts/GarG4Particle.h"

using namespace gar;

GarG4Particle::GarG4Particle (const simb::MCParticle& mcp, int parentPdg, int progenitorPdg, int progenitorTrackId,
                              const vector<pair<TLorentzVector,TLorentzVector>>& positions,
    		              const vector<pair<TLorentzVector,TLorentzVector>>& momenta, const vector<int>& regions,
                              const vector<size_t>& nptsPerRegion )
{

    fG4P = new garana::G4Particle(
            mcp.NumberTrajectoryPoints(), mcp.PdgCode(), parentPdg, progenitorPdg, 
            mcp.TrackId(), mcp.Mother(), progenitorTrackId,
            gar::ProcessNameToCode(mcp.Process()),
            gar::ProcessNameToCode(mcp.EndProcess()),
            positions, momenta, regions, nptsPerRegion
          );

}//constructor
