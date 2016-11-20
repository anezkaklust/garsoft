////////////////////////////////////////////////////////////////////////
// \file GArProperties.h
//
// \brief pure virtual base interface for LAr properties
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_GArProperties_H
#define DETINFO_GArProperties_H

// C/C++ standard libraries
#include <map>
#include <string>


///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    
    class GArProperties {
    public:
      
      GArProperties(const GArProperties &)              = delete;
      GArProperties(GArProperties &&)                   = delete;
      GArProperties& operator = (const GArProperties &) = delete;
      GArProperties& operator = (GArProperties &&)      = delete;
      virtual ~GArProperties()                          = default;
      
      virtual double RadiationLength()  const = 0;
      virtual double Argon39DecayRate() const = 0;
      
      /// Atomic number of the gas
      virtual double AtomicNumber()     const = 0;
      /// Atomic mass of the gas (g/mol)
      virtual double AtomicMass()       const = 0;
      /// Mean excitation energy of the gas (eV)
      virtual double ExcitationEnergy() const = 0;
      
    protected:
      GArProperties() = default;
      
    }; // class GArProperties
  } //namespace detinfo
} // gar

#endif // DETINFO_IGArProperties_H
