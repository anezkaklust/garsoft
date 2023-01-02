
Here are the particle lists for the particles arriving at the LArIAT location with 4 different magnet configurations (+100 A, - 100 A, +60 A, - 60A).

The text file has the following variables in this order:
EventID, SpillID, TrackID, PDGid, Px, Py, Pz, x, y, z

I need
StatusCode PDGid FirstMotherEntry SecondMotherEntry FirstDaughterEntry SecondDaughterEntry Px Py Pz Energy Mass x y z ProductionTime

NOTE! 
- StatusCode PDGid FirstMotherEntry SecondMotherEntry FirstDaughterEntry SecondDaughterEntry must be integers!

- mass in GeV
- energy in GeV/c
- momentum in GeV

- files from Supraja are in MeV
E^2  = m^2 + p^2
p^2 = px^2 + py^2 + pz^2

// For example, if you want to simulate a single muon with a 5 GeV energy 
// moving only in the z direction, the entry would be
//
// 0 1
// 1 13 0 0 0 0 0. 0. 5.0 5.0011 0.105 1.0 1.0 1.0 0.0
