This area contains a framework to rapidly analyze GArSoft anatrees; the example 
given here will read a large number of anatree files and produce a smaller file,
a DST, containing events where there was a charged current interaction of an
electron neutrino in the TPC of the MPD.

The general scheme is:
	1) Write something like makeDST.cxx, which contains 3 "callback" methods and
		an enum.
	2) Modify the Makefile so as to be able to build it with the AnaTreeFrame.
	3) Run & look at your output.

Let's go through this file by file, using the example of wanting to create some 
analysis code of your own.  Call it meSoSmart.  Start by copying this whole area 
to some directory in your local work area.  Then delete CMakeLists.txt and
Makefile_please_save.  You won't need them.



AnaTreeFrame.cxx	This is the main executable.  It will call the methods in 
	makeDST.cxx, using inputs and creating output files as specified in the
	command line.  Command line arguments are described with a -h flag in your
	own code, after you have written and built it.  The framework allows for 2
	sets of input files, labeled "signal" and "background".  You are not 
	required to use both - the makeDST example does not - and the 
	signal/background terminolosgy is arbitrary.  After the files are found and
	constructed into a TChain, AnaTreeFrame creates a PlotMap object (see below)
	and an output anatree and then processes the input TChain.


Makefile			This Makefile is configured to make debug, optimized, and 
	profiling versions of the example makeDST here.  To add your meSoSmart 
	analysis code, define a new variable CODE_2 right after the definition of
	CODE_1.  Then add a line to the all: target that says (make meSoSmart_OPT) -
	be aware that the first symbol here has to be a tab, rather than any other
	kind of whitespace.  Then copy the 3 rules for how to make the 3 versions of
	makeDST into the file down a little lower and change CODE_1 to CODE_2 in the
	copied lines.  Change also the names of the targets to something like
	meSoSmart_DEB, meSoSmart_OPT and meSoSmart_PRO.  Now you will be able to 
	build, assuming you know about the make command.  Most linux systems use
	GNUmake, actually.

tools				A directory of things you want to know about before you get
	started writing analysis code.
	
	filehamna.c		To tell if a file exists or not.
	globbing.h		To find all the files that match a string (with wildcards)
						and return the names of those files as a std::vector of
						std::string.
	concat.h		To take 2 input std::vector<std::string>, say a list of 
						signal file names an-------------------------------d a list of background file names,
						and merge them into a single list of filenames, such as
						you might want to create a TChain.  For example.
	PlotMap.h		This class is a map<string,TH1*> where the string is the name
						of a TH1D, TH2D or TProfile (these classes are all
						derived from TH1).  You want to make a more complicated
						plot?  I question your presentation skills.  The []
						operator lets you get directly to the histogram.  There
						are methods to add histograms, output one to a specific
						file or to output all of them in the TDir which is
						already open.
	TreeReader.h	This set of templates simplify reading a Tbranch from a 
						TTree.  In particular, many of the TBranch in 
						GArAnaTree are std::vector, and TreeReader.h provides
						handy accessors such as operator() and also makes sure
						that you don't read the same TBranch in repeatedly if
						you don't need to.  Its use will become apparent as we
						proceed.
	old_examp.cxx	This is an example from an earlier generation of this 
						analysis framework.  It probably won't compile any more
						but is an example of how one can copy a std::vector from
						TreeReader.h into a vector in your code which will be
						in your output tree by virtue of a TTree::Branch call.
						Then you can do all sorts of manipulations on that
						vector before writing it out with a TTree::Fill().
	PDGutils.h		These two files are codebits for common analysis tasks; one
	analysisFuncs.h		is related to MC particle identities and the other is...
						not.
	PlotRatio1D.C	Macros for plotting data/MC ratios.  These are intended for
	PlotRatio2D.C		interactive use.  Both binomial and sum-in-quadrature
						uncertainty propagation is allowed for.  By the way I am
						*pretty* sure the uncertainties calculated by ROOT v6
						in the case of binomial uncertainties on weighted events
						is wrong.  Not wildly wrong, though.
	Plot2DProf.C	Plots a 2D histogram and overlays it with a profile.  Not
						always a good idea, on account of overflow & underflow
						in the histogram.
	includeROOT.h	This is a set of includes and codebits to simplify the
						writing of your callback functions.
						
makeDST.cxx			This a the "callback" code of routines and enums which you
	will modify to make your own code (after you copy it into meSoSmart.cxx).
	We shall walk through it step by step below the line.



------------------------------- This is the line -------------------------------
The first section of your meSoSmart file is 5 includes of files described above.
You might want to read includeROOT.h at least before plunging in.  For one
thing, it will tell you to stop typing std:: all over the place.

The second section is the enum CounterTag.  It's common to want to have counters
in an analysis, to see how often execution reaches a certain line in the code.
There's a vector called counters, and this enum lets you define the counters.
So, the first element of the enum might well be nEventsRead which will index you
into counters at some location which be an int to store how many events have
been read.  It is necessary that the last element in the enum is nCounterTag.  
The analysis framework has a consistency check that needs nCounterTag.

The third section is the routine printCounters(std::vector<int>,std::ofstream*).
You don't need to figure out the arguments in these routines.  That's the 
framework's job.  Instead, you write the print the counters with whatever 
messages you want.  Did you read includeROOT.h, like I told you to?  Good.  
Then you are aware of the binerr, stringEffErr and stringEffErrMicro methods, 
which compute and printout binomial uncertainties in a spiff format, as long 
as your terminal has UTF-8 encoding.  Not only does a text file get written 
(with a name specified by the -o option) but the contents of that file are 
also sent to your terminal.

The fourth section is definePlots(PlotMap&).  Again, you don't need to figure
out the arguments.  Instead, use the methods of the PlotMap.h header to add 
TH1D, TH2D or TProfile histograms.  These will appear in the anaplots TDir of
your output file.

The fifth section is 
		processTChain(TChain* inChain, PlotMap& outPlots, TTree* outTree, 
					  vector<int>& counters, int nSigFiles);
					  
The first few lines are an internal consistency check in the framework; the next
are the computation of 3 variables needed in most analyses: nReadThisFile,
nFilesRead, fileBeingRead and sigFile.  As it happens, sigFile is not needed by
the makeDST example, as it does not make use of the "background" files and 
therefor never needs to check if it is looking at an event in the "signal" 
file sample.

The TreeReader.h templates are used to access certain branches in the input
GArAnaTree.  The TPC_X, (Y,Z) scalars in GArAnaTree are the location of the
center of the MPD's TPC in whatever coordinate system the gdml file gave at
the GENIE, GEANT, readoutsim and reco stages.  The PDG branch is a vector of 
the PDG particle ID codes of all the MCParticles in the event, and the
PDGMother branch is the same for the mother particle to each particle.  If the 
MCParticle has no mother, i.e. is a primary particle from the neutrino 
interaction, PDGMother will be 0.

This analysis will read vents from many files, all generated as run number 1.
To tag a specific event for further study, I'd like to replace the run number
with a number taken from a file.  This requires doing a SetBranch of the "Run"
TBranch in the output file.  At the start of the loop over the entries of the 
input TChain, fileBeingRead is used to get an int which is the file number and
this is placed in fFileNumber so that it will be there at the TTree output.

The next few lines start a loop over the entries of the input TChain and the
calculation of nReadThisFile, nFilesRead, fileBeingRead and sigFile.  The 1st
thing is the increment of a counter for number of events read; then fFileNumber
is set, as discussed above.  On the first entry, the center of the TPC is saved.
Then we enter a loop over all the MCParticles in the event.  If the MCParticle
has PDG particle code not 0, i.e. the MCParticle has a mother particle, then
we leave this loop, with hasNUEinTPC == false.  Otherwise, we find the location
of the starting point of the MCParticle relative to the center of the TPC
and fill some TH2D histograms from PlotMap with that information.  We increment
a counter for each such primary particle, and for each such primary particle
which is a muon.  Then, if it is not a primary electron, we continue going 
through the loop of MCParticles.  If if is a primary electron, again increment
a counter and fill some histograms, and then use insideTPC(double,double,double)
to see if this is a primary electron from an interaction in the TPC gas.  When
working in the art framework, as in reco, there is a large geometry package
available for this kind of thing, but it is not possible to access that 
information from here.  Finally, if hasNUEinTPC was set, output this entry to
the output GArAnaTree.







