#include "tools/includeROOT.h"
#include "tools/TreeReader.h"
#include "tools/PlotMap.h"
#include "tools/analysisFuncs.h"
#include "tools/PDGutils.h"





//==============================================================================
//==============================================================================
enum CounterTag   {nEventsRead,	   nPrimaries, 
				   nPrimaryMuons,  nPrimaryElectrons, 
				   nEventsWritten,

nCounterTag}; // Must be last enum





//==============================================================================
//==============================================================================
void printCounters(vector<int> counters, std::ofstream* textFileOut) {



    *textFileOut << counters[nEventsRead]       << "\tevents read" << endl;
    *textFileOut << counters[nPrimaries]        << "\tprimary particles found" << endl;
    *textFileOut << counters[nPrimaryMuons]     << "\tprimary muons found" << endl;
    *textFileOut << counters[nPrimaryElectrons] << "\tprimary electrons found" << endl;
    *textFileOut << counters[nEventsWritten]	<< "\tevents written" << endl;



	return;
}




//==============================================================================
//==============================================================================
void definePlots(PlotMap& outPlots) {



	Int_t nBinsX,nBinsY;		Double_t loX,loY, hiX, hiY;
	nBinsX = nBinsY = 150;		loX = loY = -500.0;		hiX = hiY = +500.0;
	
	outPlots.addTH2D("MCPriAllXY", "XY, All MC primary vertices",	 nBinsX,loX,hiX,nBinsY,loY,hiY);
	outPlots.addTH2D("MCPriAllXZ", "XZ, All MC primary vertices",	 nBinsX,loX,hiX,nBinsY,loY,hiY);
	outPlots.addTH2D("MCPriAllZY", "ZY, All MC primary vertices",	 nBinsX,loX,hiX,nBinsY,loY,hiY);
	outPlots.addTH2D("MCPriCCeXY", "XY, CC nue MC primary vertices", nBinsX,loX,hiX,nBinsY,loY,hiY);
	outPlots.addTH2D("MCPriCCeXZ", "XZ, CC nue MC primary vertices", nBinsX,loX,hiX,nBinsY,loY,hiY);
	outPlots.addTH2D("MCPriCCeZY", "ZY, CC nue MC primary vertices", nBinsX,loX,hiX,nBinsY,loY,hiY);



	return;
}






//==============================================================================
//==============================================================================
void processTChain(TChain* inChain, PlotMap& outPlots, TTree* outTree, 
					  vector<int>& counters, int nSigFiles) {

	// Gotta check size of counter array before going any further
	if (nCounterTag>32) {
		cout << "Too many tags for counters; max is 32" << endl;
		exit(5);
	}

	// These variables needed for almost any analysis.
	int nReadThisFile    =  0;		// means, entries read so far in this TFile
	int nFilesRead       = -1;		// means, number of TFiles completely read
	bool                 sigFile;	// True for a "signal" file.
	string fileBeingRead = "";
    Long64_t             iEntry;





	// Create handy tools to get into TBranches of the input TChain
	scalarFromTree<Float_t>  TPC_X    (inChain,"TPC_X",     &iEntry);
	scalarFromTree<Float_t>  TPC_Y    (inChain,"TPC_Y",     &iEntry);
	scalarFromTree<Float_t>  TPC_Z    (inChain,"TPC_Z",     &iEntry);
	vectorFromTree<Int_t>	 PDG      (inChain,"PDG",       &iEntry);
	vectorFromTree<Int_t>	 PDGMother(inChain,"PDGMother", &iEntry);
	vectorFromTree<Float_t>  MCPStartX(inChain,"MCPStartX", &iEntry);
	vectorFromTree<Float_t>  MCPStartY(inChain,"MCPStartY", &iEntry);
	vectorFromTree<Float_t>  MCPStartZ(inChain,"MCPStartZ", &iEntry);





	// In this job, I wish to reset the run number to the number of
	// the input file number.
	Float_t Tulsa[3] = {0,0,0};			
    int fFileNumber;
	outTree->SetBranchAddress("Run",&fFileNumber);





	// Loop over the entries of the input TChain, setting the variables needed for
	// almost any analysis: nReadThisFile, , nFilesRead, sigFile
    Long64_t nEntry = inChain->GetEntries();
    for (iEntry=0; iEntry<nEntry; ++iEntry) {
		inChain->GetEntry(iEntry);
		bool sigFile = (inChain->GetTreeNumber() < nSigFiles);
		string nameTFile = inChain->GetFile()->GetName();
		if (nameTFile != fileBeingRead) {
			if (iEntry!=0) {
				cout << " found " << nReadThisFile << " entries." << endl;
			}
			cout << "Reading " << nameTFile << "...";
			fileBeingRead = nameTFile;
			nReadThisFile = 0;
			++nFilesRead;
		}
		++nReadThisFile;





		++counters[nEventsRead];

		if (nReadThisFile == 1) {
			// Replace the Run number with the file number.  Assumes file number is
			// proceeded by the only underscore in the full file name, including
			// directory.  Get the input file number into inFileNumber, move it to
			// the memory location set for the TBranch "Run".
			size_t strtLoc = fileBeingRead.find("_") +1;
			size_t lenLoc  = fileBeingRead.find(".root") -strtLoc;
			string numberAsString = fileBeingRead.substr(strtLoc,lenLoc);
			int inFileNumber;
			sscanf(numberAsString.c_str(), "%d", &inFileNumber);
			fFileNumber = inFileNumber;
		}

		if (iEntry==0) {
			Tulsa[0] = TPC_X();	Tulsa[1] = TPC_Y();	Tulsa[2] =  TPC_Z();
		}

		int nMCParts = PDG.size();
		bool hasNUEinTPC = false;
		for (int iMCPart=0; iMCPart<nMCParts; ++iMCPart) {
			Int_t MC_Mother     = PDGMother(iMCPart);
			if (abs(MC_Mother) != 0)           break;

			Float_t XvertMC = MCPStartX(iMCPart) -Tulsa[0];
			Float_t YvertMC = MCPStartY(iMCPart) -Tulsa[1];
			Float_t ZvertMC = MCPStartZ(iMCPart) -Tulsa[2];

			++counters[nPrimaries];
			outPlots["MCPriAllXY"]->Fill(XvertMC,YvertMC);
			outPlots["MCPriAllXZ"]->Fill(XvertMC,ZvertMC);
			outPlots["MCPriAllZY"]->Fill(ZvertMC,YvertMC);

			Int_t MC_Me         = PDG(iMCPart);
			if (abs(MC_Me)     == muonPDG)   ++counters[nPrimaryMuons];
			if (abs(MC_Me)     != electronPDG) continue;

			++counters[nPrimaryElectrons];
			outPlots["MCPriCCeXY"]->Fill(XvertMC,YvertMC);
			outPlots["MCPriCCeXZ"]->Fill(XvertMC,ZvertMC);
			outPlots["MCPriCCeZY"]->Fill(ZvertMC,YvertMC);


			if (!insideTPC(XvertMC,YvertMC,ZvertMC)) continue;
			// Now you have a primary electron in the TPC.  Stop looking through
			// the MCParticle list
			hasNUEinTPC = true;
			break;
		}

		// If you found a primary electron in the MCParticle list, write the entry
		if (hasNUEinTPC) {
			outTree->Fill();
			++counters[nEventsWritten];
		}
    }	// End of loop over iEntry





	cout << " found " << nReadThisFile << " entries." << endl << endl;
    return;
}
