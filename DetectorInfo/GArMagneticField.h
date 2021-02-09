/// \file MPDMagneticField.h
/// \brief Describe the magnetic field structure of a detector
///
/// \version $Id: MPDMagneticField.h,v 1.3 2019-11-06 02:46:38 trj Exp $
/// \author trj@fnal.gov
//////////////////////////////////////////////////////////////////////////
/// \namespace mag
/// A namespace for simulated magnetic fields
//////////////////////////////////////////////////////////////////////////

#ifndef MAG_GARMAGNETICFIELD_H
#define MAG_GARMAGNETICFIELD_H

// nug4 libraries
#include "nug4/MagneticField/MagneticField.h"

//ROOT includes
#include "TH1.h"

namespace fhicl { class ParameterSet; }
namespace mag
{

// Specifies the magnetic field over all space
//
// The default implementation, however, uses a nearly trivial,
// non-physical hack.
class GArMagneticField: public MagneticField {

public:
    
    explicit GArMagneticField(fhicl::ParameterSet const& pset);
    GArMagneticField(GArMagneticField const&) = delete;
    virtual ~GArMagneticField() = default;

    void reconfigure(fhicl::ParameterSet const& pset);

    std::vector<MagneticFieldDescription> const& Fields() const override { return fFieldDescriptions; }
    size_t NumFields() const override { return fFieldDescriptions.size(); }
    MagFieldMode_t const& UseField(size_t f) const override { return fFieldDescriptions[f].fMode; }
    std::string const& MagnetizedVolume(size_t f) const override { return fFieldDescriptions[f].fVolume; }

    // return the field at a particular point
    G4ThreeVector const FieldAtPoint(G4ThreeVector const& p = G4ThreeVector(0)) const override;

    // This method will only return a uniform field based on the input
    // volume name.  If the input volume does not have a uniform field
    // caveat emptor
    G4ThreeVector const UniformFieldInVolume(std::string const& volName) const override;

private:
    void MakeCheckPlots();
    void ReadRZFile(const std::string& filename, RZFieldMap& rzmap);
    void ReadXYZFile(const std::string& filename, XYZFieldMap& xyzmap);

    G4ThreeVector CalcRZField(G4ThreeVector const& p, RZFieldMap const& rzmap) const;
    G4ThreeVector CalcXYZField(G4ThreeVector const& p, XYZFieldMap const& xyzmap) const;

    // bool fEnableField;
    // bool fUseUniformField;

    float fGlobalScaleFactor;

    std::vector<MagneticFieldDescription> fFieldDescriptions; ///< Descriptions of the fields

    // parameters for check plots

    std::vector<int> fNBinsXCheckPlots;
    std::vector<int> fNBinsYCheckPlots;
    std::vector<int> fNBinsZCheckPlots;
    std::vector<float> fXLowCheckPlots;
    std::vector<float> fXHighCheckPlots;
    std::vector<float> fYLowCheckPlots;
    std::vector<float> fYHighCheckPlots;
    std::vector<float> fZLowCheckPlots;
    std::vector<float> fZHighCheckPlots;
    std::vector<float> fXCentCheckPlots;
    std::vector<float> fYCentCheckPlots;
    std::vector<float> fZCentCheckPlots;
    std::vector<TH1*> fCheckPlots;
};

class Interpolator
{
public:
    Interpolator();

    float interpolate(const float* point, const std::vector<std::vector<std::vector<float>>>& g,
                      const float* delta, const float* offset) const;
    float interpolate(float x, float y, float z,
                      const std::vector<std::vector<std::vector<float>>>& g, float hx, float hy,
                      float hz, float xo, float yo, float zo) const;

private:
    float conv_kernel(float s) const;
};

} // namespace mag

#endif // MAG_GARMAGNETICFIELD_H
