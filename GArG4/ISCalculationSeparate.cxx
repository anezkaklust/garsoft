////////////////////////////////////////////////////////////////////////
/// \file  ISCalcuationSeparate.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4EmSaturation.hh"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "GArG4/ISCalculationSeparate.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

namespace gar {
  namespace garg4{
    
      //----------------------------------------------------------------------------
    ISCalculationSeparate::ISCalculationSeparate(CLHEP::HepRandomEngine&)
    {
    }
    
      //----------------------------------------------------------------------------
    ISCalculationSeparate::~ISCalculationSeparate()
    {
    }
    
    //----------------------------------------------------------------------------
    void ISCalculationSeparate::Initialize()
    {
      //const detinfo::GArProperties*      garp    = gar::providerFrom<detinfo::GArPropertiesService>();
      const detinfo::DetectorProperties* detprop = gar::providerFrom<detinfo::DetectorPropertiesService>();
      
//      double density       = detprop->Density(detprop->Temperature());
      fEfield              = detprop->Efield();
//      fGeVToElectrons      = lgpHandle->GeVToElectrons();
      
      // the recombination coefficient is in g/(MeVcm^2), but
      // we report energy depositions in MeV/cm, need to divide
      // Recombk from the GArG4Parameters service by the density
      // of the argon we got above.
//      fRecombA             = lgpHandle->RecombA();
//      fRecombk             = lgpHandle->Recombk()/density;
//      fModBoxA             = lgpHandle->ModBoxA();
//      fModBoxB             = lgpHandle->ModBoxB()/density;
//      fUseModBoxRecomb     = lgpHandle->UseModBoxRecomb();
      
      // Use Birks Correction in the Scintillation process
      fEMSaturation = G4LossTableManager::Instance()->EmSaturation();
      
      // TODO: make maxsize a configurable parameter
      double maxsize = 1. * CLHEP::cm;
      
      fStepSize = 0.1 * maxsize;
      
      return;
    }
    
    //----------------------------------------------------------------------------
    // fNumIonElectrons returns a value that is not corrected for life time effects
    void ISCalculationSeparate::Reset()
    {
      fEnergyDeposit   = 0.;
      fNumScintPhotons = 0.;
      fNumIonElectrons = 0.;
      
      return;
    }
    
    //----------------------------------------------------------------------------
    // fNumIonElectrons returns a value that is not corrected for life time effects
    void ISCalculationSeparate::CalculateIonizationAndScintillation(const G4Step* step)
    {
      fEnergyDeposit = step->GetTotalEnergyDeposit()/CLHEP::MeV;
      
      // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
      // R = A/(1 + (dE/dx)*k)
      // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
      // from GeV/voxel width
      // A = 0.800 +/- 0.003
      // k = (0.097+/-0.001) g/(MeVcm^2)
      // the dx depends on the trajectory of the step
      // k should be divided by the density as our dE/dx is in MeV/cm,
      // the division is handled in the constructor when we set fRecombk
      // B.Baller: Add Modified Box recombination - ArgoNeuT result submitted to JINST
      
      G4ThreeVector totstep = step->GetPostStepPoint()->GetPosition();
      totstep -= step->GetPreStepPoint()->GetPosition();
      
      double dx     = totstep.mag()/CLHEP::cm;
      double recomb = 0.;
      double dEdx   = fEnergyDeposit/dx;
      
      // Guard against spurious values of dE/dx. Note: assumes density of LAr
      if(dEdx < 1.) dEdx = 1.;
      
      if(fUseModBoxRecomb) {
        if (dx){
          double Xi = fModBoxB * dEdx / fEfield;
          recomb = log(fModBoxA + Xi) / Xi;
        }
        else
          recomb = 0;
      }
      else{
        recomb = fRecombA/(1. + dEdx * fRecombk / fEfield);
      }
      
      
      // 1.e-3 converts fEnergyDeposit to GeV
      fNumIonElectrons = fGeVToElectrons * 1.e-3 * fEnergyDeposit * recomb;
      
      LOG_DEBUG("ISCalculationSeparate")
      << " Electrons produced for " << fEnergyDeposit
      << " MeV deposited with "     << recomb
      << " recombination: "         << fNumIonElectrons;
      
      // Now do the scintillation
      G4MaterialPropertiesTable* mpt = step->GetTrack()->GetMaterial()->GetMaterialPropertiesTable();
      if( !mpt)
        throw cet::exception("ISCalculationSeparate")
        << "Cannot find materials property table"
        << " for this step! "
        << step->GetTrack()->GetMaterial();
      
      // TODO: fix this for gaseous argon
      fNumScintPhotons = std::numeric_limits<float>::min();
      
      return;
    }
    
  }// namespace
} // gar
