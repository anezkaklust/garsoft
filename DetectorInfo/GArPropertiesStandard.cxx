////////////////////////////////////////////////////////////////////////
//
//  GArProperties
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <cmath>
#include <iostream>

// GArSoft includes
#include "DetectorInfo/GArPropertiesStandard.h"
#include "CoreUtils/ProviderUtil.h" // gar::IgnorableProviderConfigKeys()

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

//-----------------------------------------------
gar::detinfo::GArPropertiesStandard::GArPropertiesStandard()
: fIsConfigured(false)
{
}

//-----------------------------------------------
gar::detinfo::GArPropertiesStandard::GArPropertiesStandard(fhicl::ParameterSet   const& pset,
                                                      std::set<std::string>        ignore_params /* = {} */)
: GArPropertiesStandard()
{
  this->Configure(pset, ignore_params);
}

#if 0
//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Configure(fhicl::ParameterSet const& pset)
{  
  this->SetRadiationLength  (pset.get< double >("RadiationLength" ));
  this->SetAtomicNumber     (pset.get< double >("AtomicNumber"));
  this->SetAtomicMass       (pset.get< double >("AtomicMass"));
  this->SetMeanExcitationEnergy (pset.get< double >("ExcitationEnergy"));

  this->SetArgon39DecayRate  (pset.get< double >("Argon39DecayRate"));

  this->SetFastScintEnergies (pset.get< std::vector<double> >("FastScintEnergies"));
  this->SetFastScintSpectrum (pset.get< std::vector<double> >("FastScintSpectrum"));
  this->SetSlowScintEnergies (pset.get< std::vector<double> >("SlowScintEnergies"));
  this->SetSlowScintSpectrum (pset.get< std::vector<double> >("SlowScintSpectrum"));
  this->SetAbsLengthEnergies (pset.get< std::vector<double> >("AbsLengthEnergies"));
  this->SetAbsLengthSpectrum (pset.get< std::vector<double> >("AbsLengthSpectrum"));
  this->SetRIndexEnergies    (pset.get< std::vector<double> >("RIndexEnergies"   ));
  this->SetRIndexSpectrum    (pset.get< std::vector<double> >("RIndexSpectrum"   ));
  this->SetRayleighEnergies  (pset.get< std::vector<double> >("RayleighEnergies" ));
  this->SetRayleighSpectrum  (pset.get< std::vector<double> >("RayleighSpectrum" ));

  this->SetScintResolutionScale (pset.get<double>("ScintResolutionScale"));
  this->SetScintFastTimeConst   (pset.get<double>("ScintFastTimeConst"));
  this->SetScintSlowTimeConst   (pset.get<double>("ScintSlowTimeConst"));
  this->SetScintBirksConstant   (pset.get<double>("ScintBirksConstant"));
  this->SetScintByParticleType  (pset.get<bool>("ScintByParticleType"));
  this->SetScintYield   (pset.get<double>("ScintYield"));
  this->SetScintPreScale(pset.get<double>("ScintPreScale"));
  this->SetScintYieldRatio      (pset.get<double>("ScintYieldRatio"));

  if(ScintByParticleType()){
    this->SetProtonScintYield(pset.get<double>("ProtonScintYield"));
    this->SetProtonScintYieldRatio   (pset.get<double>("ProtonScintYieldRatio"));
    this->SetMuonScintYield  (pset.get<double>("MuonScintYield"));
    this->SetMuonScintYieldRatio     (pset.get<double>("MuonScintYieldRatio"));
    this->SetPionScintYield  (pset.get<double>("PionScintYield"));
    this->SetPionScintYieldRatio     (pset.get<double>("PionScintYieldRatio"));
    this->SetKaonScintYield  (pset.get<double>("KaonScintYield"));
    this->SetKaonScintYieldRatio     (pset.get<double>("KaonScintYieldRatio"));
    this->SetElectronScintYield      (pset.get<double>("ElectronScintYield"));
    this->SetElectronScintYieldRatio (pset.get<double>("ElectronScintYieldRatio"));
    this->SetAlphaScintYield (pset.get<double>("AlphaScintYield"));
    this->SetAlphaScintYieldRatio    (pset.get<double>("AlphaScintYieldRatio"));
  }
  
  this->SetEnableCerenkovLight  (pset.get<bool>("EnableCerenkovLight"));
  
  this->SetReflectiveSurfaceNames(pset.get<std::vector<std::string> >("ReflectiveSurfaceNames"));
  this->SetReflectiveSurfaceEnergies(pset.get<std::vector<double> >("ReflectiveSurfaceEnergies"));
  this->SetReflectiveSurfaceReflectances(pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceReflectances"));
  this->SetReflectiveSurfaceDiffuseFractions(pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceDiffuseFractions"));

  fIsConfigured = true;


  return true;
}
#endif // 0

//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Configure(fhicl::ParameterSet   const& pset,
                                                    std::set<std::string>        ignore_params /* = {} */)
{
  std::set<std::string> ignorable_keys = gar::IgnorableProviderConfigKeys();
  ignorable_keys.insert(ignore_params.begin(), ignore_params.end());
  
  // validation happens here:
  fhicl::Table<Configuration_t> config_table { pset, gar::IgnorableProviderConfigKeys() };
  Configuration_t const& config = config_table();
  
  SetRadiationLength      (config.RadiationLength()     );
  SetAtomicNumber         (config.AtomicNumber()        );
  SetAtomicMass           (config.AtomicMass()          );
  SetMeanExcitationEnergy (config.MeanExcitationEnergy());
  SetArgon39DecayRate     (config.Argon39DecayRate()    );
  SetFanoFactor           (config.FanoFactor()          );
  fIsConfigured = true;

  return true;
}

//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Update(uint64_t ts) 
{
  if (ts == 0) return false;

  return true;
}

