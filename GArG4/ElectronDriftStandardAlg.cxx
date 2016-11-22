//
//  ElectronDriftStandardAlg.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4Step.hh"
#include "CLHEP/Random/RandGauss.h"

#include "GArG4/ElectronDriftStandardAlg.h"
#include "GArG4/IonizationAndScintillation.h"
#include "GArG4/ParticleListAction.h"
#include "DetectorInfo/DetectorProperties.h"
#include "SimulationDataProducts/SimChannel.h"

namespace gar {
  namespace garg4 {
    
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
    void ElectronDriftStandardAlg::DriftElectronsToReadout(const G4Step* step)
    {
      // Reset the Ionization and Scintillation calculator for this step
      gar::garg4::IonizationAndScintillation::Instance()->Reset(step);

      // figure out how many electrons were ionized in the step
      // if none, just return the empty vector
      const int electrons = garg4::IonizationAndScintillation::Instance()->NumberIonizationElectrons();

      if(electrons <= 0) return;

      // we know we have a step worth recording, so now get the amount of energy deposited
      const double energy = garg4::IonizationAndScintillation::Instance()->EnergyDeposit();

      // We need a random number thrower for the diffusion, etc
      CLHEP::RandGauss GaussRand(fEngine);
      
      G4ThreeVector midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
                                      step->GetPostStepPoint()->GetPosition() ) / CLHEP::cm;
      
      double g4time = step->GetPreStepPoint()->GetGlobalTime();
      
      // the XYZ position of the midpoint has to be drifted to the readout,
      // use the diffusion coefficients to decide where the electrons end up
      float xyz[3]     = {(float)midPoint[0],
                          (float)midPoint[1],
                          (float)midPoint[2]};
      float xyzDiff[3] = {0.};
      
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
      std::vector<int  > nElec(nClusters, fElectronsPerCluster);
      
      // the last cluster may have fewer electrons than the configured amount
      nElec.back() = nElectrons - (nClusters - 1) * fElectronsPerCluster;
     
      float transDiffSigma = sqrtDriftT * fTransDiffConst;
      
      GaussRand.fireArray(nClusters, &XDiff[0], 0., sqrtDriftT * fLongDiffConst);
      GaussRand.fireArray(nClusters, &YDiff[0], 0., transDiffSigma);
      GaussRand.fireArray(nClusters, &ZDiff[0], 0., transDiffSigma);
      
      // Get the track ID for the current particle.
      // We aren't getting the track ID from the step because we may
      // have had to add an offset to it if we have multiple neutrinos or
      // cosmic rays going through the detector in the same event
      auto trackID = ParticleListAction::GetCurrentTrackID();
      
      sdp::TDCIDE tdcIDE;
      auto        channel = std::numeric_limits<unsigned int>::max();

      for(size_t c = 0; c < nClusters; ++c){
        
        xyzDiff[0] = xyz[0];
        xyzDiff[1] = xyz[1] + (float)YDiff[c];
        xyzDiff[2] = xyz[2] + (float)ZDiff[c];
        
        // first get the channel number for this cluster
        try{
          channel = fGeo->NearestChannel(xyzDiff);
          
          // see if we can find an entry in the set for this one
          // already
          
          sdp::IDE ide(trackID,
                       nElec[c],
                       xyz[0],
                       xyz[1],
                       xyz[2],
                       energy * nElec[c] / nElectrons);

          sdp::SimChannel tempChan(channel);
          tempChan.AddIonizationElectrons((unsigned int)fClock.Ticks(fTime->G4ToElecTime(g4time +
                                                                                         driftT +
                                                                                         (float)XDiff[c] * fInverseVelocity)
                                                                     ),
                                          ide
                                          );
          
          auto itr = fChannels.find(tempChan);
          
          // TODO: this process of merging and erasing can probably be done more efficiently
          if(itr != fChannels.end()){
            // found the channel, now add these electrons to it
            // use the temporary one to add to, and then replace
            // the one in the set with the temporary one to get
            // around const-ness issues
            tempChan.MergeSimChannel(*itr, 0);
            fChannels.erase(itr);
          }
          
          // no need for an else statement to the effect
          // that this is a new channel, just put it in the
          // vector (we erased any channel that was already there
          // anyway)
          fChannels.emplace(std::move(tempChan));
          
        }
        catch(cet::exception &e){
          LOG_WARNING("ElectronDriftStandardAlg")
          << "Unable to locate channel relating to position ("
          << xyzDiff[0]
          << ","
          << xyzDiff[1]
          << ","
          << xyzDiff[2]
          << ") corresponding to "
          << nElec[c]
          << " electrons";
        }
      }
      
      return;
    }
    
  }
}