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
    void SetOrigin(float *origin);

    /* Set the origin with hardcoded values for the CDR production geometry */
    void SetTPCOrigin();

    /* Check if MCP started in tracker */
    bool hasOriginInTracker(TVector3 spoint);

    /* Check if MCP decayed in calo */
    bool hasDecayedInCalo(TVector3 epoint);

    /* Check if the mcp is a backscatter */
    bool isBackscatter(TVector3 spoint, TVector3 epoint);

    /* Check if it is a Bremsstrahlung photon */
    bool isBremsstrahlung(TVector3 spoint, int pdg, int motherpdg);

    float GetRamdomNumber() { return _rando->Rndm(); }

    float GaussianSmearing(float mean, float sigma) { return _rando->Gaus(mean, sigma); }

    float* GetOrigin() { return &_origin[0]; }

    float* GetOriginTPC() { return &_originTPC[0]; }

private:
    float _origin[3];                   ///< coordinates of the origin
    float _originTPC[3];                ///< coordinates of the origin of TPC used for the CDR production
    unsigned long int _seed;         ///< seed
    TRandom3 *_rando;                ///< random generator
};

#endif
