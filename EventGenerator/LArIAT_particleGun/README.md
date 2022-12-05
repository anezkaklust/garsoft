
Here are the particle lists for the particles arriving at the LArIAT location with 4 different magnet configurations (+100 A, - 100 A, +60 A, - 60A).

The text file has the following variables in this order:
EventID, SpillID, TrackID, PDGid, Px, Py, Pz, x, y, z

I need
StatusCode PDGid FirstMotherEntry SecondMotherEntry FirstDaughterEntry SecondDaughterEntry Px Py Pz Energy Mass x y z ProductionTime

// For example, if you want to simulate a single muon with a 5 GeV energy 
// moving only in the z direction, the entry would be
//
// 0 1
// 1 13 0 0 0 0 0. 0. 5.0 5.0011 0.105 1.0 1.0 1.0 0.0
