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

