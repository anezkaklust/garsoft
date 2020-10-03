#ifndef UTILS_H
#define UTILS_H 1

#include "TVector3.h"
#include "TRandom3.h"

class Utils
{
public:

    /* Default constructor */
    Utils();

    /* Copy constructor */
    Utils(const Utils &) = default;

    /* Destructor */
    ~Utils();

    /* Set the seed */
    void SetSeed(int seed);

    /* Set the origin */
    void SetOrigin(double *origin);

    /* Check if MCP point is in fiducial */
    bool PointInFiducial(const TVector3 &point);

    /* Check if MCP point is in TPC Volume (can be outside Fid) */
    bool PointInTPC(const TVector3 &point);

    /* Check if MCP point is in the ECAL */
    bool PointInCalo(const TVector3 &point);

    /* Check if MCP point is in between the TPC and the ECAL */
    bool PointStopBetween(const TVector3 &point);

    /* Check if MCP point is not in the fiducial and not in the ECAL */
    bool isThroughCalo(const TVector3 &point);

    /* Check if MCP decayed in calo */
    bool hasDecayedInCalo(const TVector3 &point);

    /* Check if MCP is in the Barrel region */
    bool isBarrel(const TVector3 &point);

    /* Check if MCP is in the Endcap region */
    bool isEndcap(const TVector3 &point);

    /* Check if the mcp is a backscatter */
    bool isBackscatter(const TVector3 &spoint, const TVector3 &epoint);

    /* Check if it is a Bremsstrahlung photon */
    bool isBremsstrahlung(const TVector3 &spoint, const int& pdg, const int& motherpdg);

    float GetRamdomNumber() { return _rando->Rndm(); }

    float GaussianSmearing(const float& mean, const float& sigma) { return _rando->Gaus(mean, sigma); }

    double* GetOrigin() { return &_origin[0]; }

private:
    double _origin[3];                   ///< coordinates of the origin
    unsigned long int _seed;         ///< seed
    TRandom3 *_rando;                ///< random generator

    /* HARDCODED */
    double _TPCFidRadius = 222.5;
    double _TPCFidLength = 215.;
    double _TPCRadius = 273.;
    double _TPCLength = 259.;
    double _ECALInnerRadius = 278.;
    double _ECALOuterRadius = 321.;
    double _ECALStartX = 364.;
    double _ECALEndX = 406.;
};

#endif
