#ifndef GAR_ANAUTILS_H
#define GAR_ANAUTILS_H

#include <string>
#include <unordered_map>

#include "garsoft/ReconstructionDataProducts/Cluster.h"
#include "garana/DataProducts/CaloCluster.h"

#include "garsoft/ReconstructionDataProducts/Vertex.h"
#include "garana/DataProducts/Vertex.h"

#include "garsoft/ReconstructionDataProducts/Vee.h"
#include "garana/DataProducts/Vee.h"

namespace gar{

    //extern std::unordered_map<std::string,int> processMap;
    int ProcessNameToCode(std::string const& p);

    const garana::CaloCluster MakeAnaCalCluster(const rec::Cluster& clust, const std::vector<std::pair<int,float>>& edeps);
    const garana::Vee         MakeAnaVee(const rec::Vee& vee);
    const garana::Vertex      MakeAnaVtx(const rec::Vertex& vtx);
    

}

#endif
