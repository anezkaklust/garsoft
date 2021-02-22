
#include "AnalysisDataProducts/MCDisplayTrack.h"
#include <iostream>

 using namespace gar::adp;

 size_t MCDisplayTrack::NPoints() const{
     return fPositions.size();
 }

 float MCDisplayTrack::X(const size_t i) const{

     if( i >= this->NPoints() || i < 0){
         std::cout << "WARNING MCDisplayTrack::X():"
                   << " passed index out of range...default to last traj point"
                   << std::endl;
         return fPositions.back().X();
     }

     return fPositions[i].X();
 }

 float MCDisplayTrack::Y(const size_t i) const{

     if( i >= this->NPoints() || i < 0){
         std::cout << "WARNING MCDisplayTrack::Y():" 
                   << " passed index out of range...default to last traj point"
                   << std::endl;
         return fPositions.back().Y();
     }

     return fPositions[i].Y();
 }

 float MCDisplayTrack::Z(const size_t i) const {

     if( i >= this->NPoints() || i < 0){
         std::cout << "WARNING MCDisplayTrack::Z():" 
                   << " passed index out of range...default to last traj point"
                   << std::endl;
         return fPositions.back().Z();
     }

     return fPositions[i].Z();
 }

 float MCDisplayTrack::T(const size_t i) const{

     if( i >= this->NPoints() || i < 0){
         std::cout << "WARNING MCDisplayTrack::T():"
                   << " passed index out of range...default to last traj point"
                   << std::endl;
         return fPositions.back().T();
     }

     return fPositions[i].T();
 }

 //we pass by copy to avoid modifying original traj object
 void MCDisplayTrack::Sparsify(simb::MCTrajectory traj) {

     traj.Sparsify(0.1); //configurable margin for points displaced from linear interp

     //loop over trajectory points
     for(size_t i=0; i<traj.size(); i++){
         fPositions.push_back(traj.Position(i));
     }

 }


