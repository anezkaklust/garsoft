////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNEST.h
/// \brief Interface to algorithm class for a specific calculation of 
///        ionization electrons and scintillation photons using NEST
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GARG4ISCALCULATIONNEST_H
#define GARG4ISCALCULATIONNEST_H

#include "GArG4/ISCalculation.h"
#include "GArG4/NestAlg.h"

// forward declaration
namespace CLHEP { class HepRandomEngine; } 

namespace gar {
  namespace garg4 {
    
    class ISCalculationNEST : public ISCalculation {
      
    public:
      
      ISCalculationNEST(CLHEP::HepRandomEngine& engine);
      virtual ~ISCalculationNEST();
      
      void   Initialize();
      void   Reset();
      void   CalculateIonizationAndScintillation(const G4Step* step);
      double EnergyDeposit()              const { return fEnergyDeposit;   }
      int    NumberIonizationElectrons()  const { return fNumIonElectrons; }
      int    NumberScintillationPhotons() const { return fNumScintPhotons; }
      double StepSizeLimit()              const { return fStepSize;        }
      
    private:
      
      NestAlg* fNest;                  ///< the fast optical simulation process
      double   fStepSize;              ///< maximum step to take
      CLHEP::HepRandomEngine& fEngine; ///< random engine
    };
  }
} // gar

#endif // GARG4ISCALCULATIONNEST_H

