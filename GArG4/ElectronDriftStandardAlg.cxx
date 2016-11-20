//
//  ElectronDriftStandardAlg.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "GArG4/ElectronDriftStandardAlg.h"
#include "DetectorInfo/DetectorProperties.h"
#include "SimulationDataProducts/SimChannel.h"

namespace gar {
  namespace garg4 {
    
    //--------------------------------------------------------------------------
    ElectronDriftStandardAlg::ElectronDriftStandardAlg(CLHEP::HepRandomEngine      & engine,
                                                       fhicl::ParameterSet    const& pset)
    : EletronDriftAlg(engine, pset)
    {
      fElectronsPerCluster     = pset.get<int>("ElectronsPerCluster"    );
      fMinimumElectronClusters = pset.get<int>("MinimumElectronClusters");
    }
    
    //------------------------------------------------------------------------
    ElectronDriftStandardAlg::~ElectronDriftStandardAlg()
    {
    }

    //--------------------------------------------------------------------------
    std::vector<sdp::TDCIDE> ElectronDriftStandardAlg::DriftElectronsToReadout(const G4Step* step)
    {
      std::vector<sdp::TDCIDE> tdcIDEs;

      // figure out how many electrons were ionized in the step
      // if none, just return the empty vector
      const int electrons = garg4::IonizationAndScintillation::Instance()->NumberIonizationElectrons();

      if(electrons <= 0) return tdcIDEs;

      // we know we have a step worth recording, so now get the amount of energy deposited
      const double energy = garg4::IonizationAndScintillation::Instance()->EnergyDeposit();

      // We need a random number thrower for the diffusion, etc
      CLHEP::RandGauss GaussRand(fEngine);
      
      G4ThreeVector midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
                                      step->GetPostStepPoint()->GetPosition() ) / CLHEP::cm;
      
      double g4time = step->GetPreStepPoint()->GetGlobalTime();
      
      // the XYZ position of the midpoint has to be drifted to the readout,
      // use the diffusion coefficients to decide where the electrons end up
      double x = midPoint[0];
      double y = midPoint[1];
      double z = midPoint[2];
      
      // The geometry has x = 0 at the readout, so x is conveniently the
      // drift distance for this step and the drift time is x / fDriftVelocity
      double driftT       = x / fDriftVelocity;
      double sqrtDriftT   = std::sqrt(driftT);
      double lifetimeCorr = std::exp(driftT / fLifetimeCorrection);
      
      // how many electrons do we expect to make it to the readout?
      double nElectrons = electrons * lifetimeCorr;
      size_t nClusters  = (size_t)std::ceil(nElectrons / fElectronsPerCluster);

      // drift them in sets
      std::vector<double> XDiff(nClusters);
      std::vector<double> YDiff(nClusters);
      std::vector<double> ZDiff(nClusters);
      std::vector<int   > nElec(nClusters, fElectronsPerCluster);
      
      // the last cluster may have fewer electrons than the configured amount
      nElec.back() = nElectrons - (nClusters - 1) * fElectronsPerCluster;
     
      double transDiffSigma = sqrtDriftT * fTransDiffConst);
      
      GausRand.fireArray(nClusters, &XDiff[0], 0., sqrtDriftT * fLongDiffConst);
      GausRand.fireArray(nClusters, &YDiff[0], y,  transDiffSigma);
      GausRand.fireArray(nClusters, &ZDiff[0], z,  transDiffSigma);
      
      sdp::TDCIDE tdcIDE;
      
      for(size_t c = 0; c < nClusters; ++c){
        
        

        // We aren't getting the track ID from the step because we may
        // have had to add an offset to it if we have multiple neutrinos or
        // cosmic rays going through the detector in the same event

        tdcIDE.fTDC = 0;
        tdcIDE.fIDEs.emplace_back(ParticleListAction::GetCurrentTrackID(),
                                   ,
                                   x,
                                   y,
                                   z,
                                   energy);
      }
      
      return tdcIDEs;
    }
    
  }
}