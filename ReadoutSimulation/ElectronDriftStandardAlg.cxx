//
//  ElectronDriftStandardAlg.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGauss.h"

#include "ReadoutSimulation/ElectronDriftStandardAlg.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "DetectorInfo/DetectorProperties.h"
#include "SimulationDataProducts/EnergyDeposit.h"

namespace gar {
  namespace rosim {
    
    //--------------------------------------------------------------------------
    ElectronDriftStandardAlg::ElectronDriftStandardAlg(CLHEP::HepRandomEngine      & engine,
                                                       fhicl::ParameterSet    const& pset)
    : ElectronDriftAlg(engine, pset)
    {
      fElectronsPerCluster = pset.get<int>("ElectronsPerCluster");
      fMinClusters         = pset.get<size_t>("MinClusters");

      // get service providers
      fGeo   = gar::providerFrom<geo::GeometryGAr>();
      fTime  = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();
      fClock = fTime->TPCClock();
      
    }
    
    //------------------------------------------------------------------------
    ElectronDriftStandardAlg::~ElectronDriftStandardAlg()
    {
    }

    //--------------------------------------------------------------------------
    void ElectronDriftStandardAlg::DriftElectronsToReadout(gar::sdp::EnergyDeposit       const& dep,
                                                           gar::rosim::ElectronDriftInfo      & driftInfo)
    {
      // Reset the Ionization and Scintillation calculator for this step
      gar::rosim::IonizationAndScintillation::Instance()->Reset(&dep);

      // figure out how many electrons were ionized in the step
      // if none, just return the empty vector
      const int electrons = gar::rosim::IonizationAndScintillation::Instance()->NumberIonizationElectrons();

      if (electrons <= 0) return;

      // We need a random number thrower for the diffusion, etc
      CLHEP::RandGauss GaussRand(fEngine);
      
      double g4time = dep.Time();
      
      // the XYZ position of the midpoint has to be drifted to the readout,
      // use the diffusion coefficients to decide where the electrons end up
      float xyz[3] = {dep.X(),
                      dep.Y(),
                      dep.Z()};

      MF_LOG_DEBUG("ElectronDriftStandardAlg::DriftElectronsToReadout") << "energy deposition: " 
        << dep.X() << " " << dep.Y() << " " << dep.Z() << std::endl;

      // The geometry has x = TPCXCent() at the cathode plane. Compute the drift time for a two-drift-volume
      // TPC.  If we only have one drift volume, need to get that information here.

      //double driftD       = std::max(0.0,26.35-xyz[0]); 
      double driftD       = std::max(0.0,xyz[0]-3.85);  // position of the OROC is at 3.85
      //std::cout<<"DriftD " << driftD <<std::endl;
      if (driftD > 26.35){ // the drift region between the cathode and the OROC
        //std::cout<< "Outside of the drift region, not drifted towards the OROC." <<std::endl;
        return;
      }
      //double driftD       = std::abs(std::abs(xyz[0] - fGeo->TPCXCent()) - fGeo->TPCLength()/2.0);  // assume cathode is at x=TPXCent
      //if (fGeo->TPCNumDriftVols() == 1)
	//{
	  //driftD = std::max(0.0, fGeo->TPCLength()/2.0 - (xyz[0] - fGeo->TPCXCent()));
	//}
      double driftT       = driftD * fInverseVelocity;  // in cm
      double sqrtDriftD   = std::sqrt(driftD);
      double lifetimeCorr = std::exp(driftT / fLifetimeCorrection);
      
      // how many electrons do we expect to make it to the readout?
      float  nElectrons = electrons * lifetimeCorr;
      size_t nClusters  = (size_t)std::ceil(nElectrons / fElectronsPerCluster);
      nClusters = std::min( (size_t) std::ceil(nElectrons),std::max(nClusters,fMinClusters));
      int ourelectronspercluster = nearbyint(nElectrons/nClusters);

      // drift them in sets.  The X(Y)(Z)Diff vectors contain the additional
      // distance to add to the step midpoint to account for diffusion.
      std::vector<double> XDiff(nClusters, 0.);
      std::vector<double> YDiff(nClusters, 0.);
      std::vector<double> ZDiff(nClusters, 0.);
      std::vector<double> TDiff(nClusters, 0.);
      std::vector<int   > nElec(nClusters, ourelectronspercluster);
      
      // the last cluster may have fewer electrons than the configured amount
      nElec.back() = nearbyint(nElectrons - (nClusters - 1) * ourelectronspercluster);
     
      double longDiffSigma = sqrtDriftD * fLongDiffConst;
      double transDiffSigma = sqrtDriftD * fTransDiffConst;

      GaussRand.fireArray(nClusters, &XDiff[0], 0., longDiffSigma);
      GaussRand.fireArray(nClusters, &YDiff[0], 0., transDiffSigma);
      GaussRand.fireArray(nClusters, &ZDiff[0], 0., transDiffSigma);
      
      // set the times and the positions of each cluster
      for(size_t c = 0; c < nClusters; ++c){
        TDiff[c]  = g4time + driftT + XDiff[c] * fInverseVelocity;
        XDiff[c]  = xyz[0];
        YDiff[c] += xyz[1];
        ZDiff[c] += xyz[2];

	  MF_LOG_DEBUG("ElectronDriftStandardAlg::DriftElectronsToReadout")
        << "g4 time: "
        << g4time
        << " drift time "
        << driftT
        << " diffusion  x:"
        << XDiff[c]
        << " y "
        << YDiff[c]
        << " z "
        << ZDiff[c]
        << " t "
        << TDiff[c]
        << " tick period "
        << fTime->TPCClock().TickPeriod();

      }
      
      // reset the ElectronDriftInfo object with this information
      driftInfo.Reset(XDiff, YDiff, ZDiff, TDiff, nElec);
      
      return;
    }
   
    //--------------------------------------------------------------------------

  }
}
