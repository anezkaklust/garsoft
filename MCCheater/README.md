The BackTracker service must be, as of the ides of April 2020, instantiated with a const\_cast:

`	`cheat::BackTrackerCore const\* const\_bt = gar::providerFrom<cheat::BackTracker>();

`	`cheat::BackTrackerCore\*             bt = const\_cast<cheat::BackTrackerCore\*>(const\_bt);

It also must be instantiated in the analyze method of the *art* plugin derived from EDAnalyzer, rather than the beginJob method.


The base class, BackTrackerCore, has a set of stdlib structures that are filled upon art’s sPreProcessEvent.

Then an example use is:

`	`art::Handle< std::vector<rec::Cluster> > RecoClusterHandle;

`	`if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {

`		`throw cet::exception("anatest") << " No rec::Cluster branch."

`		`<< " Line " << \_\_LINE\_\_ << " in file " << \_\_FILE\_\_ << std::endl;

`	`}

`	`for ( rec::Cluster cluster : \*RecoClusterHandle ) {

`		`std::vector<std::pair<simb::MCParticle\*,float>> whatMatches;

`		`whatMatches = bt->ClusterToMCParticles(&cluster);

`		`std::cout << "\nCluster No. " << cluster.getIDNumber() << " is made of MCParticles: " << std::endl;

`		`for (auto itr = whatMatches.begin(); itr!=whatMatches.end(); ++itr) {

`			`std::cout << "G4 track number " << itr->first->TrackId() << "\thas PDG code " <<

`				`itr->first->PdgCode() << "\tand its mother is G4 track " << itr->first->Mother()

`				`<< "\tand energy fraction " << 100\*(itr->second) << "%\n";

`		`}

`	`}





*The following should always work* ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

explicit BackTrackerCore(fhicl::ParameterSet const& pset);			But of course.  Default destructor.

~BackTrackerCore();

void AdoptEveIdCalculator(sim::EveIdCalculator\* ec);				One of the functions of the ParticleList class is to find Eve, the mother of all

`														`mother particles.  The algorithm for determining Eve is declared with this

`														`method; the only algorithm that exists now is the nutools EmEveIdCalculator

`														`which moves up the tree until it finds a particle not produced by one of the

`														`GEANT processes {conv, LowEnConversion, Pair, compt, Compt, Brem,

`														`phot, Photo, Ion, Annihil}.  Other functions of ParticleList are less useful

`														`for us.





*The following should work as long as there are* MCParticle *and* MCTruth *data products in the event, with an* art::Assns *between them* ---------------------------------------------

sim::ParticleList\* GetParticleList() const;						The backtracker uses the nutools ParticleList class (more recently, in nug4)

`														`as an augmentation of std::map<int, simbMCParticle\*>, where the int

`														`is the GEANT trackID.  This method returns a pointer to the instance of the ParticleList 

`														`class which is in BackTrackerCore and is filled by the RebuildNoSC method.

simb::MCParticle\* const

`	`TrackIDToParticle(int const& id) const;						From a GEANT TrackID to an MCParticle.  The reverse process is just 

`														`MCParticle::TrackID().  Negative GEANT TrackIDs are believed to be

`														`EM shower particles.  Will return nullptr on error condition.

simb::MCParticle\* const

`	`FindMother(simb::MCParticle\* const p) const;					From an MCParticle to its mother particle.  The MCParticle class has a 

`														`Mother() method, but it returns a track ID not an MCParticle.

simb::MCParticle\* const

`	`FindEve(simb::MCParticle\* const p) const;					From an MCParticle to its Eve MCParticle, using the adopted EveIdCalculator.

`														`Will return nullptr on error condition in TrackIDToParticle.

simb::MCParticle\* const										From an MCParticle to the progenitor MCParticle which exited the TPC or is a primary.

`	`FindTPCEve(simb::MCParticle\* const p);						Uses the GeometryCore::PointInMPD to find if one end of the particle is in the TPC and 

`														`the other is not.  Will return nullptr on error condition in TrackIDToParticle.

simb::MCParticle\* const										Signature using trackID mostly for call in CellIDToCalIDEs but may as well make it

`	`FindTPCEve(int const trackID);							public.  This is the method that prevents BackTrackerCore from being const.

bool IsForebearOf(simb::MCParticle\* const p, simb::MCParticle\* const c) const;

`														`Is MCParticle c in the descent tree of MCParticle p? Does NOT use the

`														`adopted EveIdCalculator; just walks up the tree.  MCParticle equality is 

`														`tested by comparison of the MCParticle::TrackID() fields.  Returns true

`														`if p==c.

art::Ptr<simb::MCTruth> const									From an MCParticle to the MCTruth that it originates from.

`	`ParticleToMCTruth(simb::MCParticle\* const p) const;





*The following should work as long as there are also* RawDigit *and*  EnergyDeposit *data products in the event, with an* art::Assns *between them* -------------------------------

std::vector<HitIDE>

`	`HitToHitIDEs(art::Ptr<rec::Hit> const& hit) const;			This is the HitToTrackID method used in tpcpatreccheat\_module.cc,

std::vector<HitIDE>											but we change the name for clarity.  The returned cheat::HitIDE is a

`	`HitToHitIDEs(         rec::Hit  const& hit) const;			struct defined in BackTrackerCore.h and it includes the GEANT track ID,

`														`a fraction of the RawDigit energy from the particle with that track ID, and an

`														`“ionization” which seems to be charge in that channel but outside a certain time

`														`window.

std::vector<art::Ptr<rec::Hit>> const							Returns all the TPC hits that have some energy from a particle.  If you want to find the

`	`ParticleToHits(simb::MCParticle\* const p,					fraction of the energy in a hit from that particle, you need a HitIDE, found with

`				  `std::vector<art::Ptr<rec::Hit>> const& allhits,	HitToHitIDEs, above and this method only returns hits with that fraction over

`				  `checkNeutrals=false) const;					the value given by fMinHitEnergyFraction, which is specified in BackTracker.fcl;

`														`the default value of 0.1 is probably fine.  I guess.  If checkNeutrals==false, will

`														`return immediately for input particle is a photon, neutron, pi0 or either a mu or e

`														`neutrino thereby returning an empty vector promptly.

std::vector<HitIDE>

`	`ChannelToHitIDEs(raw::Channel\_t const& channel,				The above 3 public methods are wrappers to this private one.  It works by looping

`				  `double const start, double const stop) const);	through all the EnergyDeposits for the channel as given by fChannelToEDepCol,

`														`which in turn is filled by RebuildNoSC.

std::pair<double,double>										Input is a single MCParticle and a collection of hits in the TPC; the return value is

`	`HitPurity(simb::MCParticle\* const p,						the fraction, and the binomial error on the fraction, of the hits that are from the MCParticle.

`              `std::vector<art::Ptr<rec::Hit>> const& hits,			Warning: the standard binomial error formula is dubious when used with weights.

`              `bool weightByCharge=false) const;

std::pair<double,double>										Input is a single MCParticle and two collections of hits in the TPC.  The allhits

`	`HitEfficiency(simb::MCParticle\* const p,					collection should be all the hits in the event, or maybe all the plausible hits in the event.

`                  `std::vector<art::Ptr<rec::Hit>> const& hits,		The intent is that allhits is all the hits that MCParticle created.  The hits collection

`				  `std::vector<art::Ptr<rec::Hit>> const& allhits,	might be, say, all the hits in a reconstructed track from that MCParticle.  The returned

`                  `bool weightByCharge=false) const;				pair of doubles is the efficiency of track reconstruction and the binomial uncertainty.  Again,

`														`the standard binomial error formula is dubious when used with weights.





*The following should work as long as there are also* CaloRawDigit *and*  CaloDeposit *data products in the event, with an* art::Assns *between them* ---------------------------



*First we’ll need to define an analog to* HitsIDE:

struct CalIDE {

`      `int   trackID;      ///< Geant4 supplied trackID

`      `float energyFrac;   ///< fraction of CaloRawDigit energy from the particle with this trackID

`      `CalIDE()

`      `{}

`      `CalIDE(int id, float ef)

`      `: trackID(id), energyFrac(ef)

`      `{}

};



*Second, create a* typedef raw::CellID *for* long long int*.  It’s in* CaloRawDigit.h*.  Including this is a small but ugly addition to our dependencies.  We should collect*

*all these typedefs in one place.  I’ve put that on the task list.*



std::vector<CalIDE>											Similar to TPCHitToHitIDEs except that the CaloDeposit::TrackID should be the

`	`CaloHitToCalIDEs(art::Ptr<rec::CaloHit> const& hit); 			track ID of the particle coming into the calorimeter.  That chase up the tree has to be done in

std::vector<CalIDE>											CellIDToCalIDEs.  There is reportedly a similar functionality in edepsim.

`	`CaloHitToCalIDEs(         rec::CaloHit  const& hit);



std::vector<art::Ptr<rec::CaloHit>>							Returns all the calorimeter hits that have some energy from TPC-leaving particle (see

`	`ParticleToCaloHits(simb::MCParticle\* const p					FindTPCEve).  Input allhits is all plausible hits in the ECAL calorimeter that

`				    `std::vector<art::Ptr<rec::CaloHit>> const& allhits) const;

`														`this method will search through.

std::vector<CalIDE>&										The above 3 methods are wrappers to this private one.  Works by looping through all the 

`	`CellIDToCalIDEs(raw::CellID\_t const& cellID, float const time);	CaloDeposits for the CellID as given by fCellIDToCaloDepCol.  Uses FindTPCEve.

std::pair<double,double>										Corresponding to HitCollectionPurity, except on the calorimeter side.

`	`CaloHitPurity(simb::MCParticle\* const p,

`                   `vector<art::Ptr<rec::CaloHit>> const& hits,

`                   `bool weightByCharge=false) const;

std::pair<double,double>										Corresponding to HitCollectionEfficiency, except on the calorimeter side.

`	`CaloHitEfficiency(simb::MCParticle\* const p,					

`				   `vector<art::Ptr<rec::CaloHit>> const& hits,

`				      `std::vector<art::Ptr<rec::CaloHit>> const& allhits,

`				   `bool weightByCharge=false) const;				                      





*The following should work as long as there are also* Track *data products in the event, with (presumably chained)* art::Assns *back to the* Hits* --------------------------------------

std::vector<art::Ptr<rec::Hit>> const							Uses the chained associations to get a collection of Hits for this Track.

`	`TrackToHits(rec::Track\* const t);

std::vector<std::pair<simb::MCParticle\*,float>>					Uses TrackToHits and then HitToHitIDEs to get energy fractions for each 

`	`TrackToMCParticles(rec::Track\* const t); 					MCParticle that might match the track; and returns that list of candidate MCParticles

`														`sorted by their contribution.  Electrons that have parents that are one of the other 

`														`candidate MCParticles get their ionization included in with that parent.  The float

`														`in the returned vector is the fraction of the reconstructed track’s energy attributed to

`														`each MCParticle.

vector<art::Ptr<rec::Track>>									The inverse operation to TrackToMCParticles.  Just loops over the input array of

`	`MCParticleToTracks(simb::MCParticle\* const p,				reconstructed tracks and see which tracks have more than fTrackFracMCP of their

`				    `vector<art::Ptr<rec::Track>> const& tracks);	ionization energy from the input MCParticle.  However, because of the merging

`														`of electrons in TrackToMCParticles, will return a zero-length vector if the argument 

`														`is one of those electrons.




*The following should work as long as there are also* Cluster *data products in the event, with (presumably chained)* art::Assns *back to the* CaloHits* ---------------------------

std::vector<art::Ptr<rec::CaloHit>> const						Uses the chained associations to get a collection of CaloHits for this Cluster.	ClusterToCaloHits(rec::Cluster\* const c);

std::vector<std::pair<simb::MCParticle\*,float>>					Uses ClusterToHits and then CaloHitToCalIDEs to get energy fractions for each  

`	`ClusterToMCParticles(rec::Cluster\* const c);					MCParticle that might match the cluster; and returns that list of candidate MCParticles

`														`sorted by their contribution.    Gets the MCParticle entering the calorimeter using

`														`FoundTPCEve.

vector<art::Ptr<rec::Cluster>>

`	`MCParticleToClusters(simb::MCParticle\* const p,				The inverse operation to ClusterToMCParticles.

`					  `vector<art::Ptr<rec::Cluster>> const& clusters);

bool ClusterCreatedMCParticle(simb::MCParticle\* const p,			Uses MCPartsInCluster to get a collection of MCParticles in the cluster.  Then uses 

`                              `rec::Cluster\* const c);				IsDescendedFrom to see if MCParticle p is descended from any of those MCParticles.

`														`If p’s ionization created c, then the answer is yes!



bool MCParticleCreatedCluster(simb::MCParticle\* const p,			Uses MCPartsInCluster to get a collection of MCParticles in the cluster.  Then uses 

`                              `rec::Cluster\* const c);				IsDescendedFrom to see if any of those MCParticles are descended from MCParticle p.

`														`If p’s ionization created c, then the answer is yes!



std::vector<simb::MCParticle\*> MCPartsInCluster(rec::Cluster\* const c);

`														`A private method used by ClusterCreatedMCParticle and MCParticleCreatedCluster

`														`to use ClusterToCaloHits, CaloHitToCalIDEs and TrackIDToParticle to get a collection 

`                    										`of MCParticles in the cluster.
