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
      
      // get service providers
      fGeo   = gar::providerFrom<geo::Geometry>();
      fTime  = gar::providerFrom<detinfo::DetectorClocksService>();
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

      if(electrons <= 0) return;

      // We need a random number thrower for the diffusion, etc
      CLHEP::RandGauss GaussRand(fEngine);
      
      double g4time = dep.Time();
      
      // the XYZ position of the midpoint has to be drifted to the readout,
      // use the diffusion coefficients to decide where the electrons end up
      float xyz[3]     = {dep.X(),
                          dep.Y(),
                          dep.Z()};
      
      // The geometry has x = 0 at the readout, so x is conveniently the
      // drift distance for this step and the drift time is x / fDriftVelocity
      float driftT       = xyz[0] * fInverseVelocity;
      float sqrtDriftT   = std::sqrt(driftT);
      float lifetimeCorr = std::exp(driftT / fLifetimeCorrection);
      
      // how many electrons do we expect to make it to the readout?
      float  nElectrons = electrons * lifetimeCorr;
      size_t nClusters  = (size_t)std::ceil(nElectrons / fElectronsPerCluster);

      // drift them in sets.  The X(Y)(Z)Diff vectors contain the additional
      // distance to add to the step midpoint to account for diffusion.
      std::vector<double> XDiff(nClusters);
      std::vector<double> YDiff(nClusters);
      std::vector<double> ZDiff(nClusters);
      std::vector<double> TDiff(nClusters);
      std::vector<int   > nElec(nClusters, fElectronsPerCluster);
      
      // the last cluster may have fewer electrons than the configured amount
      nElec.back() = nElectrons - (nClusters - 1) * fElectronsPerCluster;
     
      float transDiffSigma = sqrtDriftT * fTransDiffConst;
      
      GaussRand.fireArray(nClusters, &XDiff[0], 0., sqrtDriftT * fLongDiffConst);
      GaussRand.fireArray(nClusters, &YDiff[0], 0., transDiffSigma);
      GaussRand.fireArray(nClusters, &ZDiff[0], 0., transDiffSigma);
      
      // set the times and the positions of each cluster
      for(size_t c = 0; c < nClusters; ++c){
        TDiff[c]  = g4time + driftT + (float)XDiff[c] * fInverseVelocity;
        XDiff[c]  = xyz[0];
        YDiff[c] += xyz[1];
        ZDiff[c] += xyz[2];
      }
      
      // reset the ElectronDriftInfo object with this information
      driftInfo.Reset(XDiff, YDiff, ZDiff, TDiff, nElec);
      
      return;
    }
   
    //--------------------------------------------------------------------------

  }
}