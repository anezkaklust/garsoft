#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
using std::cout;		using std::endl;
using std::string;		using std::vector;

#include "tools/filehamna.h"
#include "tools/globbing.h"
#include "tools/concat.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "tools/PlotMap.h"





//==============================================================================
//==============================================================================
// For each analysis, create a ---.cxx file with these methods
void definePlots(PlotMap& outPlots);

void processTChain(TChain* inChain, PlotMap& outPlots, TTree* outTree, 
					  vector<int>& counters, int nSigFiles);

void printCounters(vector<int> counters, std::ofstream* textFileOut);





//========================================================================S======
//==============================================================================
int main(int argc , const char* argv[]){



	// Parse the input arguments
	bool haveInFile = false;
	string globListFile;
	string outName = "";
    bool haveOutFile = false;

	int	nSigFiles, nBkgFiles, nFiles;
	std::vector<string> sigList;	std::vector<string> bkgList;   std::vector<string> fileList;

	int iArg = 1;
	while ((iArg < argc) && (argv[iArg][0] == '-')){
		switch (argv[iArg][1]) {
		case 'i':
			haveInFile = true;
			++iArg;
			globListFile = argv[iArg];
			break;
		case 'I':
			haveInFile = true;
			++iArg;
			fileList.push_back(argv[iArg]);
			break;
		case 'o':
			// Cleverly, argv comes to us with environmental variables expanded
			// and we don't need wordexp.h.
			outName = argv[iArg+1];
			haveOutFile = true;			
			++iArg;
			break;

		case 'h':
		default:
			cout << "-I filename Process just the file with this name" << endl;
			cout << "-i filename File listing input globbing patterns.  1st line in" << endl;
			cout << "               file is globbing pattern for \"signal\" files; 2nd" << endl;
			cout << "               line is for \"data\" files; if either line is e.g." << endl;
			cout << "               \"dCache */OverlaidDST_*.root\" then will use xrootd" << endl;
			cout << "               to get all the files matching */OverlaidDST_*.root" << endl;
			cout << "               from /pnfs/dune/persistent/users/bellanto/.  Edit" << endl;
			cout << "               the value of dCacheDIR near line 1290 of" << endl;
			cout << "               AnaTreeFrame.cxx to pull from a different directory." << endl;
			cout << "               If the 2nd line is just \"end\" then only \"signal\"" << endl;
			cout << "               files will be processed." << endl;
            cout << "[-o name]   Name for output files.  Must be specified even if there" << endl;
			cout << "                will be nothing in it.  The files created will be" << endl;
			cout << "               <name>.root and <name>.txt.  The former will have a" << endl;
			cout << "               TDir named \"anatree\" containing the produced" << endl;
			cout << "               GArAnaTree and another TDir named \"anaplots\" " << endl;
			cout << "               containing histograms defined in the callback method" << endl;
			cout << "               definePlots(PlotMap&).  Plots and GArAnaTree are" << endl;
			cout << "               filled in the callback routine" << endl;
			cout << "               processTChain(TChain*,PlotMap&,TTree,std::vector<int>,int)" << endl;
			cout << "               The <name>.txt file contains the values in the counter" << endl;
			cout << "               array, as defined with the CounterTag enum definition," << endl;
			cout << "               filled in the callback processTChain and printed in the" << endl;
			cout << "               callback method printCounters(vector<int>, std::ofstream*)" << endl;
			cout << "               Environmental variables in <name> will be expanded." << endl;
			cout << "[-h]        This message" << endl;
			exit(0);
		}
		++iArg;
    }
	if (!haveInFile) {
		cout << "You need an -i or -I argument in the command line" << endl;
		cout << "Use -h for quick help" << endl;
		exit(1);
	}
	if (!haveOutFile) {
		cout << "You need an -o argument in the command line" << endl;
		cout << "Use -h for quick help" << endl;
		exit(1);
	}





	// Get the list of input files
	string dCacheDir  = "/pnfs/dune/persistent/users/bellanto/";
	if (fileList.size()==0) {
		string lookFor;
		std::ifstream globListInStream(globListFile);

		globListInStream >> lookFor;
		if (lookFor == "dCache") {
			globListInStream >> lookFor;
			lookFor = dCacheDir + lookFor;
			sigList = globbing(lookFor);
			string newAccess = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
			std::vector<string>::iterator it = sigList.begin();
			for (; it<sigList.end(); ++it) {
				it->erase(0,5);
				*it = newAccess + *it;
			}
		} else {
			sigList = globbing(lookFor);
		}
		cout << ( nSigFiles = (int)sigList.size() ) << " \"signal\" file(s) found.\n";

		globListInStream >> lookFor;
		if (lookFor == "dCache") {
			globListInStream >> lookFor;
			lookFor = dCacheDir + lookFor;
			bkgList = globbing(lookFor);
			string newAccess = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
			std::vector<string>::iterator it = bkgList.begin();
			for (; it<bkgList.end(); ++it) {
				it->erase(0,5);
				*it = newAccess + *it;
			}
		} else {
			string endLower = lookFor;
			for (size_t i=0; i<endLower.length(); ++i) endLower[i]=tolower(endLower[i]);
			if (endLower != "end") {
				bkgList = globbing(lookFor);
			}
		}
		globListInStream.close();
		cout << ( nBkgFiles = (int)bkgList.size() ) << " \"background\" file(s) found.\n";
		fileList = concat<string>(sigList,bkgList);
	} else {
		// Processing only one file ( -I option ); call it signal.
		nSigFiles = (int)fileList.size();
	}

	cout << ( nFiles = (int)fileList.size() ) << " total files found.\n";
	if (nFiles==0) exit(2);
	cout << endl;





	// Output counters defined, output plots definition called back
	vector<int> counters;
	counters.resize(32);
	PlotMap outPlots;
	definePlots(outPlots);





	// Build a TChain for input.  Because the TTrees are copied using CloneTree(0) 
	// and outTree->Fill() it is necessary to create the input file and TTree
	// (well TChain) before creating the output file and TTree.  Because ROOT is
	// so great.
	bool hosed = false;
	cout << "Globbing the list of input files." << endl;
	TChain* inChain = new TChain("anatree/GArAnaTree");
	for (int iFile=0; iFile<nFiles; ++iFile) {
		if (filehamna(fileList[iFile].c_str())) {
			if (!hosed) {
				cout << fileList[iFile] << " does not exist!  Skipping it." << endl;
			} else {
				// Only allow 2 errors in getting the input file list
				cout << fileList[iFile] << " also does not exist!  Crashing, burning." << endl;
				exit(3);
			}
		}
        cout << "Adding file : " << fileList[iFile].c_str() << endl;
		inChain->Add(fileList[iFile].c_str());
     }





    // Output file made and opened here
	string RootOutName  = outName + ".root";
    if (!filehamna(RootOutName)) {
        cout << "Already existing " << RootOutName << endl;
        exit(4);
    }
	TFile* outFile = TFile::Open(RootOutName.c_str(),"NEW");
    if (outFile->IsZombie()) {cout << "Eat yo brains!" << endl; return 1;};
    outFile->mkdir("anatree");
    outFile->cd("anatree");
    TTree* outTree = new TTree("GArAnaTree","GArAnaTree");
 	outTree = inChain->CloneTree(0);		// 0 argument prevents copy of entries.





    // Go figger!
	time_t unixT;	struct tm* localT;
	time(&unixT);	localT = localtime(&unixT);
	cout << "\nStarting at " << asctime(localT) << endl;

	processTChain(inChain, outPlots, outTree, counters, nSigFiles);

	localT = localtime(&unixT);
	cout << "Ending at " << asctime(localT) << endl;





	// Ship those puppies out the door & split
	outFile->Write();

	string TextOutName = outName;
	TextOutName += ".txt";
	std::ofstream textFileOut(TextOutName.c_str());

	printCounters(counters, &textFileOut);

    textFileOut.close();
	string stdOut = "cat " + TextOutName;
	system( stdOut.c_str() );

	outFile->cd("/");
	outFile->mkdir("anaplots");
	outFile->cd("anaplots");
    outPlots.writePlots();
	outFile->Close();

    exit(0);
}
