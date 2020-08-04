//
//  TPCCluster.cxx
//
//  Created by Brian Rebel on 10/6/16.
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "TMath.h"

namespace gar {
  namespace rec {
   
    //--------------------------------------------------------------------------
    TPCCluster::TPCCluster()
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();
      return;
    }



    bool TPCCluster::operator==(const TPCCluster& rhs) const {
      return (this->fIDnumero == rhs.fIDnumero);
    }

    bool TPCCluster::operator!=(const TPCCluster& rhs) const {
      return (this->fIDnumero != rhs.fIDnumero);
    }

    gar::rec::IDNumber TPCCluster::getIDNumber() const {return fIDnumero;}



    //--------------------------------------------------------------------------
    TPCCluster::TPCCluster( const float        sig,
                            const float       *pos,
                            const float        startT,
                            const float        endT,
                            const float        Time,
                            const float        RMS,
                            const float       *cov)
      : fSignal   (sig   )
      , fTime     (Time  )
      , fStartTime(startT)
      , fEndTime  (endT  )
      , fRMS      (RMS   )
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();
      
      fPosition[0] = pos[0];
      fPosition[1] = pos[1];
      fPosition[2] = pos[2];
      for (size_t i=0; i<6; ++i)
        {
          fCovMat[i] = cov[i];
        }
      
      return;
    }
    
    //--------------------------------------------------------------------------
    // Users should be careful with this method.  It will just add TPCClusters 
    // without checking if the TPCClusters are overlapping in time or close
    // to each other in time
    void TPCCluster::operator+=(gar::rec::TPCCluster const& h)
    {
      
      // add h to this TPCCluster.  Set the new position, etc to be the weighted
      // average of the individual TPCClusters
      float totSig = fSignal + h.Signal();
      
      if(totSig == 0.){
        MF_LOG_WARNING("TPCCluster")
          << "attempting to add two TPCClusters and neithr has any signal, bail.";
        return;
      }
 
      for(size_t i = 0; i < 3; ++i)
        fPosition[i] = (fPosition[i] * fSignal + h.Position()[i] * h.Signal()) / totSig;

      // don't do a weighted average but just make a wider TPCCluster
      
      //fStartTime = (fStartTime * fSignal + h.StartTime() * h.Signal()) / totSig;
      //fEndTime   = (fEndTime   * fSignal + h.EndTime()   * h.Signal()) / totSig;
      
      fStartTime = TMath::Min(fStartTime,h.StartTime());
      fEndTime   = TMath::Max(fEndTime,h.EndTime());

      float avgtime = (fTime * fSignal + h.Time() * h.Signal()) / totSig;

      fRMS = TMath::Sqrt(  (fSignal*(TMath::Sq(fTime-avgtime)+TMath::Sq(fRMS))
                            + h.Signal()*(TMath::Sq(h.Time()-avgtime)+TMath::Sq(h.RMS())))/totSig );
                          
      fTime      = avgtime;

      fSignal += h.Signal();
      
      return;
    }


    //--------------------------------------------------------------------------
    std::ostream& operator<< (std::ostream& o, gar::rec::TPCCluster const& h)
    {
      const float* cov;
      cov = h.CovMatPacked();
      o << "TPCCluster "
        << "\n\tID number = "
        << h.getIDNumber()
        << "\n\tposition = ("
        << h.Position()[0]
        << ", "
        << h.Position()[1]
        << ", "
        << h.Position()[2]
        << ")"
        << "\n\tcovariance matrix = \n\t"
        << cov[0] << "\t" << cov[1] << "\t" << cov[2] << "\n"
        <<         "\t\t" << cov[3] << "\t" << cov[4] << "\n"
        <<                         "\t\t\t" << cov[5]
        << "\n\tsignal = "
        << h.Signal()
        << "\n\tdrift direction rms = "
        << h.RMS()
        << "\n\tstart time: "
        << h.StartTime()
        << "\n\tend time: "
        << h.EndTime();
      
      return o;
    }

  } // rec
} // gar
