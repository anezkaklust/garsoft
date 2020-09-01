# ParamSim

[![Build Status](https://travis-ci.org/ebrianne/ParamSim.svg?branch=master)](https://travis-ci.org/ebrianne/ParamSim)

The parametrized simulation, currently a module in GArSoft v02_04_00 and GArSoft development branch is complementary to GArSoft's full reconstruction approach. GArSoft is Multipurpose Detector's (MPD) official simulation, reconstruction, and analysis software framework (https://cdcvs.fnal.gov/redmine/projects/garsoft/wiki), and v02_04_00 of GArSoft is the version tagged for DUNE near detector Conceptual Design Report. In GArSoft development branch, you may find this module under the Ana/ParamSim directory.

Following is the reconstruction approach taken to parametrize the detector response:

1) Gluckstern relation (https://www.sciencedirect.com/science/article/pii/S0168900209009474) for reconstructing long tracks curling up in magnetic field  
2) Particle range for reconstructing short tracks that stop in the high pressure gas argon TPC (HPgTPC) component of the MPD and do not curl up in magnetic field

In addition to reconstructing tracks, a dE/dx-based PID is implemented in the module. This is based on Tom Junk's parametrization of PEP-4's dE/dx curve: https://home.fnal.gov/~trj/mpd/dedx_sep2019/ (PID matrices are in pid.root file)

*Units are in cm, GeV, ns.*
*For vectors that are not filled, a default value of -1 is set in order to keep the vector length the same (nFSP)*

The module is designed to take GArSoft's analysis tree, anatree as input and produce a so called "cafanatree" ntuple as output. A description of cafanatree tree variables are as the following:

* Event: an art event (not neutrino-interaction event rather the spill number)

* Run: run #

* SubRun: subrun #

- Generator-level information (more information about generator-level variables can be found at: https://internal.dunescience.org/doxygen/classsimb_1_1MCNeutrino.html)

  * nFSP: number of monte-carlo particle in the event, not only primaries but all particles (includes shower particles)

  * mode: interaction mode

  * q2: momentum transfer

  * w: hadronic invariant mass

  * y: inelasticity

  * x: Bjorken-x variable

  * theta: angle between incoming and outgoing lepton

  * t: a variable used in coherent pion analysis

  * ntype: neutrino PDG

  * ccnc: charged current (ccnc == 0) or neutral current (ccnc == 1)

  * gint: a variable used in coherent pion analysis

  * tgtpdg: target PDG (nucleus-level variable)

  * weight: number of protons-on-target

  * gt_t: a variable used in coherent pion analysis

  * intert: interaction type

  * mcnupx: x-component momentum of neutrino in GeV/c

  * mcnupy: y-component momentum of neutrino in GeV/c

  * mcnupz: z-component momentum of neutrino in GeV/c

  * vertx: vertex of primary mc particle in x in cm

  * verty: vertex of primary mc particle in y in cm

  * vertz: vertex of primary mc particle in z in cm

- GEANT4-level information:

  * pdgmother: pdg code of the particle that created the particle under consideration

  * mctrkid: number created in G4 to track the particles (relations mother <-> daughters). The original neutrino has a track id of -1 and then it is incremented by G4.

  * motherid: id number associated with the mother mc particle -> returns the index in the MCP vector of the mother of this particle

  * mctime: detector response time/time information of the particle with respect to neutrino interaction time (which is 0 nano seconds)

  * MCPStartX: starting position of the particle along the x-direction (starting position of the track) in cm

  * MCPStartY: starting position of the particle along the y-direction (starting position of the track) in cm

  * MCPStartZ: starting position of the particle along the z-direction (starting position of the track) in cm

  * truepx: truth-level x-component momentum of the particle in GeV/c

  * truepy: truth-level y-component momentum of the particle in GeV/c

  * truepz: truth-level z-component momentum of the particle in GeV/c

  * MCPEndX: end position of the particle along the x-direction (end position of the track) in cm

  * MCPEndY: end position of the particle along the y-direction (end position of the track) in cm

  * MCPEndZ: end position of the particle along the z-direction (end position of the track) in cm

  * MCProc: a vector containing a string of the GEANT4 process that created a particle (e.g. for particles emerging from a neutrino interaction, the GEANT4 process is "primary")

  * MCEndProc: a vector containing a string of the GEANT4 end process of a particle (e.g. for e+e- pairs emerging from a photon conversion, the end process is "conv")

  * angle: truth-level angle with respect to the beam direction (along the z)

  * truep: truth-level momentum of the particle in GeV/c

  * truepdg: truth-level PDG code of the particle

- Reconstruction-level information (note that this is different from GArSoft reconstruction -- see above for more information about the parametrized reconstruction in the ParamSim module):

  * recopid: reconstructed PDG code of the particle using the parametrized dE/dx

  * trkLen: length of the track

  * trkLenPerp: length of the track in the direction transverse to the beam

  * preco: reconstructed particle momentum

  * anglereco: reconstructed angle of the particle with respect to the beam direction (wrt the z-direction)

  * erecon: reconstructed energy of neutral and charged particles that reach the ECAL (when track length != 0 and endpoint is in the calo). Particles that are not reconstructed in the ECAL have erecon of 0. Calculated from ECAL parametrisation using E = sqrt( p*p + m*m ) - m (kinetic energy)

  * prob_arr: array of PID scores/probabilities (using log-likelihood) for reconstructed particles

  * recopidecal: reconstructed PDG code of the neutral and charged particles with the ECAL. A value of 0 is considered if the particle does not reach ECAL.

  * detected: whether the particle is seen by the ECAL

  * etime: smeared mc time by 1 ns, reco by the ECAL (if reaches or seen in ECAL only)

- Geometry information (simple flags 0 or 1)

  * isFidStart/End: Check if the particle start/end point is in the fiducial volume of the TPC

  * isTPCStart/End: Check if the particle start/end point is in the volume of the TPC (includes fiducial and out of fiducial)

  * isCaloStart/End: Check if the particle start/end point is in the ECAL

  * isThroughCaloStart/End: Check if the particle start/end point is not in the TPC and the ECAL (assumes went through the ECAL)

  * isInBetweenStart/End: Check if the particle start/end point is in between the TPC and the ECAL

  * isBarrelStart/End: Flag to tell if the particle start/end point is in the *Barrel region*

  * isEndcapStart/End: Flag to tell if the particle start/end point is in the *Endcap region*

Check out the test directory for an example macro on how to read the cafanatree analysis ntuples that are produced as outputs of running the ParamSim module.   

**Hardcoded values!!**

```C++
  double _TPCFidRadius = 222.5;
  double _TPCFidLength = 215.;
  double _TPCRadius = 273.;
  double _TPCLength = 259.;
  double _ECALInnerRadius = 278.;
  double _ECALOuterRadius = 321.;
  double _ECALStartX = 364.;
  double _ECALEndX = 406.;
```

The volumes are defined as the following:

```C++
bool Utils::PointInFiducial(TVector3 point)
{
    bool isInFiducial = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > _TPCFidRadius ) isInFiducial = false;
    if( r_point < _TPCFidRadius && std::abs(point.X()) > _TPCFidLength ) isInFiducial = false;

    return isInFiducial;
}
```
```C++
bool Utils::PointInTPC(TVector3 point)
{
    if(PointInFiducial(point)) return true;
    bool isInTPC = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > _TPCRadius ) isInTPC = false;
    if( r_point < _TPCRadius && std::abs(point.X()) > _TPCLength ) isInTPC = false;

    return isInTPC;
}
```
```C++
bool Utils::PointInCalo(TVector3 point)
{
    bool isInCalo = false;
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    //in the Barrel
    if( r_point > _ECALInnerRadius && r_point < _ECALOuterRadius && std::abs(point.X()) < _ECALStartX ) isInCalo = true;
    //in the Endcap
    if( r_point < _ECALInnerRadius && std::abs(point.X()) > _ECALStartX && std::abs(point.X()) < _ECALEndX ) isInCalo = true;

    return isInCalo;
}
```
```C++
bool Utils::PointStopBetween(TVector3 point)
{
    //Barrel Radius 278 cm
    //Endcap starts at 364 cm
    bool isStopBetween = false;
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    //in the Barrel
    if( r_point < _ECALInnerRadius && r_point > _TPCRadius && std::abs(point.X()) < _TPCLength ) isStopBetween = true;
    //in the Endcap
    if( r_point < _ECALInnerRadius && std::abs(point.X()) > _TPCLength && std::abs(point.X()) < _ECALStartX ) isStopBetween = true;

    return isStopBetween;
}
```
```C++
bool Utils::isThroughCalo(TVector3 point)
{
    return !PointInTPC(point) && !PointStopBetween(point) && !PointInCalo(point);
}
```
```C++
bool Utils::isBarrel(TVector3 point)
{
  bool isBarrel = false;
  float theta = std::atan(_ECALInnerRadius / _ECALStartX ); //angle for barrel/endcap transition
  float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
  float theta_point = std::atan(r_point / point.X() ); //angle for barrel/endcap transition for the point

  if( theta_point > theta ) isBarrel = true;
  return isBarrel;
}
```
```C++
bool Utils::isEndcap(TVector3 point)
{
    bool isEndcap = false;
    if( !isBarrel(point) ) isEndcap = true;
    return isEndcap;
}
```
