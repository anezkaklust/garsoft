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
      const detinfo::GArProperties*      garp    = gar::providerFrom<detinfo::GArPropertiesService>();
      const detinfo::DetectorProperties* detprop = gar::providerFrom<detinfo::DetectorPropertiesService>();
      
      double density       = detprop->Density(detprop->Temperature());
      fEfield              = detprop->Efield();
      fGeVToElectrons      = garp->GeVToElectrons();
      
      // the recombination coefficient is in g/(MeVcm^2), but
      // we report energy depositions in MeV/cm, need to divide
      // Recombk from the GArG4Parameters service by the density
      // of the argon we got above.
      fRecombA = garp->RecombA();
      fRecombk = garp->Recombk()/density;
      
      // Use Birks Correction in the Scintillation process
      fEMSaturation = G4LossTableManager::Instance()->EmSaturation();
      
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
      fEnergyDeposit = step->GetTotalEnergyDeposit() / CLHEP::MeV;
      
      // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
      // R = A/(1 + (dE/dx)*k)
      // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
      // from GeV/voxel width
      // A = 0.800 +/- 0.003
      // k = (0.097+/-0.001) g/(MeVcm^2)
      // the dx depends on the trajectory of the step
      // k should be divided by the density as our dE/dx is in MeV/cm,
      // the division is handled in the constructor when we set fRecombk
      
      G4ThreeVector totstep = (step->GetPostStepPoint()->GetPosition() -
                               step->GetPreStepPoint()->GetPosition()   ) / CLHEP::cm;
      
      double dx     = totstep.mag();
      double recomb = 0.;
      double dEdx   = fEnergyDeposit / dx;
      
      // Guard against spurious values of dE/dx. Note: assumes density of GAr
      if(dEdx < 1.) dEdx = 1.;
      
      recomb = fRecombA / (1. + dEdx * fRecombk / fEfield);
      
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
