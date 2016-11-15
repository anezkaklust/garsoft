/**
 * @file    LArPropertiesStandard.h
 * @brief   Service provider with utility LAr functions
 * @see     LArPropertiesStandard.cxx LArPropertiesStandardTestHelpers.h
 * 
 * The provider detinfo::LArProperiesStandard supports simple setup for testing
 * environment, by including in your test:
 *     
 *     #include "lardata/DetectorInfo/LArPropertiesStandardTestHelpers.h"
 *     
 */

#ifndef DETECTORINFO_GARPROPERTIESSTANDARD_H
#define DETECTORINFO_GARPROPERTIESSTANDARD_H


// LArSoft libraries
#include "DetectorInfo/GArProperties.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>
#include <set>


namespace gar {
  namespace detinfo {
    /**
     * @brief Properties related to liquid argon environment in the detector
     *
     * This class can access databases via DatabaseUtil service.
     *
     * @note Some of the database connection properties are established before
     * the beginning of the job and if they change this service will not be
     * aware of it. These properties petrain, so far, only the connection mode
     * and not any content of the databases themselves.
     * @note 2: the database connection features for this base class have been removed
     */
    class GArPropertiesStandard : public GArProperties {
    public:
      GArPropertiesStandard();
      explicit GArPropertiesStandard(fhicl::ParameterSet const& pset,
                                     std::set<std::string>      ignore_params = {});
      GArPropertiesStandard(GArPropertiesStandard const&) = delete;
      virtual ~GArPropertiesStandard() = default;
      
      /**
       * @brief Configures the provider
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       *
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      bool   Configure(fhicl::ParameterSet const& pset,
                       std::set<std::string>      ignore_params = {});
      bool   Update(uint64_t ts=0);
      
      virtual double RadiationLength()  	     const override { return fRadiationLength; } ///< g/cm^2
      
      virtual double Argon39DecayRate()              const override { return fArgon39DecayRate; }  // decays per cm^3 per second
      
        /// Ar atomic number
      virtual double AtomicNumber() const override { return fZ; }
        /// Ar atomic mass (g/mol)
      virtual double AtomicMass() const override { return fA; }
        /// Ar mean excitation energy (eV)
      virtual double ExcitationEnergy() const override { return fI; }
      
      
      
      void SetRadiationLength(double rl) { fRadiationLength = rl; }
      void SetArgon39DecayRate(double r) { fArgon39DecayRate = r;}
      void SetAtomicNumber(double z) { fZ = z;}
      void SetAtomicMass(double a) { fA = a;}
      void SetMeanExcitationEnergy(double e) { fI = e;}
      
    private:
    protected:
      
        /// structure with all configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Atom<double> RadiationLength
        { Name("RadiationLength" ), Comment("radiation length [g/cm^2]") };
        fhicl::Atom<double> AtomicNumber
        { Name("AtomicNumber"    ), Comment("atomic number (yes, yes, it's 18...)") };
        fhicl::Atom<double> AtomicMass
        { Name("AtomicMass"      ), Comment("atomic mass [g/mol]") };
        fhicl::Atom<double> MeanExcitationEnergy
        { Name("ExcitationEnergy"), Comment("mean excitation energy [eV]") };
        
        fhicl::Atom<double> Argon39DecayRate
        { Name("Argon39DecayRate"), Comment("decays/(cm^3 s)") };
        
      }; // Configuration_t
      
      
      bool   fIsConfigured;
      
      double fRadiationLength;  ///< g/cm^2
      double fArgon39DecayRate; ///<  decays per cm^3 per second
      
      // Following parameters are for use in Bethe-Bloch formula for dE/dx.
      
      double fZ;                ///< Ar atomic number
      double fA;                ///< Ar atomic mass (g/mol)
      double fI;                ///< Ar mean excitation energy (eV)
      
    }; // class GArPropertiesStandard
  } //namespace detinfo
} // gar
#endif // GARPROPERTIES_H
