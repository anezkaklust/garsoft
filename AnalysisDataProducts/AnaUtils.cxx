//
//    AnaUtils.h
//
//    Created by C. Hilgenberg

//#ifndef GAR_ANAUTILS_cxx
//#define GAR_ANAUTILS_cxx

#include "AnalysisDataProducts/AnaUtils.h"
//#include "nusimdata/SimulationBase/MCParticle.h"
//#include "nug4/ParticleNavigation/ParticleList.h"

//#include <unordered_map>
#include <iostream>

//using sim::ParticleList;
//using sim::ParticleHistory;
//using simb::MCParticle;

namespace gar{

    //////////////////////////////////////////////
    const garana::CaloCluster MakeAnaCalCluster(const rec::Cluster& clust) {

       TLorentzVector pos( clust.Position()[0], clust.Position()[1], clust.Position()[2], clust.Time() );
       const float* eigs = clust.EigenVectors(); 
       vector<TVector3> vecs;
       for(int i=0; i<3; i++){
           TVector3 eig( eigs[i], eigs[i+1], eigs[i+2] );
           vecs.push_back(eig);
       }

       garana::CaloCluster outclust(pos, clust.Energy(), clust.EnergyError(), clust.TimeDiffFirstLast(),
                                    clust.Shape(), clust.ITheta(), clust.IPhi(), vecs );
        
       return outclust;
    }//

    ///////////////////////////////////////////////
    const garana::Vee         MakeAnaVee(const rec::Vee& vee) {

       TLorentzVector vtx( vee.Position()[0], vee.Position()[1], vee.Position()[2], vee.Time() );
       vector<TLorentzVector> vecs;
       for(int i=0; i<3; i++) {
           vecs.push_back(vee.FourMomentum(i));
       }

       garana::Vee outvee( vtx, vecs, vee.Chisq(), vee.VertexCov() );
        
       return outvee;
    }//
    ///////////////////////////////////////////////
    const garana::Vertex      MakeAnaVtx(const rec::Vertex& vtx) {

      TLorentzVector pos( vtx.Position()[0], vtx.Position()[0], vtx.Position()[0], vtx.Time() );
      garana::Vertex outvtx(pos, vtx.CovMat());

      return outvtx;
    }//

    ////////////////////////////////////////////////
    std::unordered_map<std::string,int> processMap{
        {"unknown"                   , 0 },
        {"primary"                   , 1 },
        {"compt"                     , 2 },
        {"phot"                      , 3 },
        {"annihil"                   , 4 },
        {"eIoni"                     , 5 },
        {"eBrem"                     , 6 },
        {"conv"                      , 7 },
        {"muIoni"                    , 8 },
        {"muMinusCaptureAtRest"      , 9 },
        {"neutronInelastic"          , 10},
        {"nCapture"                  , 11},
        {"hadElastic"                , 12},
        {"Decay"                     , 13},
        {"CoulombScat"               , 14},
        {"muPairProd"                , 15},
        {"muBrems"                   , 16},
        {"muPairProd"                , 17},
        {"PhotonInelastic"           , 18},
        {"hIoni"                     , 19},
        {"protonInelastic"           , 20},
        {"pi+Inelastic"              , 21},
        {"CHIPSNuclearCaptureAtRest" , 22},
        {"pi-Inelastic"              , 23},
        {"hBertiniCaptureAtRest"     , 24},
        {"photonNuclear"             , 25},
        {"muonNuclear"               , 26},
        {"kaon-Inelastic"            , 27},
        {"kaon+Inelastic"            , 28}, 
        {"Transportation"            , 29},
        {"nKiller"                   , 30},
        {"ionIoni"                   , 31},
        {"msc"                       , 32}
    };


    //convert process name from G4 to numeric value for
    // easy storage in TTrees
    int ProcessNameToCode(std::string const& p)
    {
        if(processMap.find(p)==processMap.end()){
            std::cout << "WARNING(AnaUtils::ProcessNameToCode): "
                      << "unknown process name, '" << p << "'" << std::endl;
            return -9999;  
        }
        
        else 
            return processMap[p];
  
    }

    //=================================================================
    /*void FillParticleList( ParticleList* list, const vector<MCParticle*>& mcps){

        for( auto const& mcp : mcps) {
            list.Add(mcp);
        }

    }

    //==================================================================
    // Given a particle's trackId, find the primary particle (progenitor)
    // that produced it and return its trackID
    int GetProgenitorId(ParticleList* list, int trackId){

        int progenitorId = -9999;

        // Look for the track in the particle list.
        ParticleList::const_iterator search = list->find( trackID );
        
        // While we're still finding particles in the chain...
        while ( search != list->end() ){

            const simb::MCParticle* particle = (*search).second;
            
            // If this is a primary particle, we're done.
            int tmpId = particle->TrackId();
            if ( list->IsPrimary( id ) ) {
                progenitorId = tmpId;
                break;
            }
            
            // Now look for the parent of this particle.
            int parentID = particle->Mother();
            search = list->find( parentID );
            
        } // while we're finding particles in the chain

        return progenitorId;

    }*/

}//namespace

//#endif

