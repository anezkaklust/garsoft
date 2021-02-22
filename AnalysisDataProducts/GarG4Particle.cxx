//
//  GarG4Particle.cxx
//
//  Written by C. Hilgenberg

#include "garsoft/AnalysisDataProducts/GarG4Particle.h"

using namespace gar;

GarG4Particle::GarG4Particle (const simb::MCParticle& mcp, int parentPdg, int progenitorPdg, int progenitorTrackId)
{

    fG4P = new garana::G4Particle(
            mcp.PdgCode(), parentPdg, progenitorPdg, 
            mcp.TrackId(), mcp.Mother(), progenitorTrackId,
            gar::ProcessNameToCode(mcp.Process()),
            gar::ProcessNameToCode(mcp.EndProcess()),
            mcp.Position(), mcp.EndPosition(),
            mcp.Momentum(), mcp.EndMomentum() 
          );

}//constructor
