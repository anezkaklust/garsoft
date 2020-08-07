void CheckSizeBranches()
{
    // CaliceStyle();

    //Read sizes of the caf tree vectors and check that the sizes are consistent 

    //if you are reading only one file, for debugging  purposes: 1) comment out lines 8 through 39 
    //and then again lines 145 and 146, 2) un-comment lines 41 and 42.   
    const char *dirname="/pnfs/dune/persistent/users/ND_GAr/2020_06_21/neutrino/";   
    
    const char *ext=".root";
    std::vector<const char*> lista;
    lista.clear();

    TChain *chain = new TChain("caf");

    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    

    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile*)next())) {
    fname = file->GetName();
    if (file->IsDirectory() && fname.EndsWith(ext)) continue;

    char result[255];   // array to hold the result.
    strcpy(result,dirname); // copy string one into the result.
    strcat(result,fname.Data()); // append string two to the result.
    
    lista.push_back(result);

    for(std::vector<const char*>::iterator it = lista.begin(); it != lista.end(); it++)
    {

    std::cout << *it << std::endl;
    string a(*it);
    TFile * tf = new TFile ( Form("%s",a.c_str() ) );
    chain->Add( Form("%s/caf", a.c_str()) );
    
    //TChain *chain = new TChain("caf");
    //chain->Add("/pnfs/dune/persistent/users/ND_GAr/2020_06_21/neutrino/neutrino.nd_hall_mpd_only.volTPCGas.Ev336000.Ev336999.11141.caf.root"); 

    //std::cout << *it << std::endl;

    int evt;
    std::vector<double>* mctime = 0;
    std::vector<int> *_nFSP = 0;
    std::vector<double> *theta = 0;
    std::vector<int> *detected = 0;
    std::vector<int> *pdgmother = 0;
    std::vector<int> *truepdg = 0;
    std::vector<int> * _MCPStartX = 0;
    std::vector<int> *_MCPStartY = 0;
    std::vector<int> *_MCPStartZ = 0;
    std::vector<int> *_MCPEndX = 0;
    std::vector<int> *_MCPEndY = 0;
    std::vector<int> *_MCPEndZ = 0;
    std::vector<std::string> *_MCProc = 0;
    std::vector<std::string> *_MCEndProc = 0;
    std::vector<float> *trkLen = 0;
    std::vector<float> *trkLenPerp = 0;
    std::vector<float> *truep = 0;
    std::vector<float> *truepx = 0;
    std::vector<float> *truepy = 0;
    std::vector<float> *truepz = 0;
    std::vector<float> *_angle = 0;
    //Reco values
    std::vector<int> *recopid = 0;
    std::vector<float> *prob_arr = 0;
    std::vector<float> *partereco = 0;
    std::vector<float> *anglereco = 0;
    std::vector<float> *_preco = 0;
    std::vector<float> *erecon = 0;

    chain->SetBranchAddress("Event", &evt);
    chain->SetBranchAddress("theta", &theta);
    chain->SetBranchAddress("nFSP", &_nFSP);
    chain->SetBranchAddress("detected", &detected);
    chain->SetBranchAddress("pdgmother", &pdgmother);
    chain->SetBranchAddress("MCPTime", &mctime);
    chain->SetBranchAddress("MCPStartX", &_MCPStartX);
    chain->SetBranchAddress("MCPStartY", &_MCPStartY);
    chain->SetBranchAddress("MCPStartZ", &_MCPStartZ);
    chain->SetBranchAddress("truepx", &truepx);
    chain->SetBranchAddress("truepy", &truepy);
    chain->SetBranchAddress("truepz", &truepz);
    chain->SetBranchAddress("MCPEndX", &_MCPEndX);
    chain->SetBranchAddress("MCPEndY", &_MCPEndY);
    chain->SetBranchAddress("MCPEndZ", &_MCPEndZ);
    chain->SetBranchAddress("MCProc", &_MCProc);
    chain->SetBranchAddress("MCEndProc", &_MCEndProc);
    chain->SetBranchAddress("angle", &_angle);
    chain->SetBranchAddress("truep", &truep);
    chain->SetBranchAddress("truepdg", &truepdg);
    //Reco info
    chain->SetBranchAddress("recopid", &recopid);
    chain->SetBranchAddress("trkLen", &trkLen);
    chain->SetBranchAddress("trkLenPerp", &trkLenPerp);
    chain->SetBranchAddress("preco", &_preco);
    chain->SetBranchAddress("anglereco", &anglereco);
    chain->SetBranchAddress("erecon", &erecon);
    chain->SetBranchAddress("prob_arr", &prob_arr);

    for(int itree = 0; itree < chain->GetEntries(); itree++)
    {
        chain->GetEntry(itree);

        std::cout << "Event " << evt << std::endl;
        std::cout << "theta " << theta->size() << std::endl;
        std::cout << "Size pdgmother: " << pdgmother->size() << std::endl;
        std::cout << "Size MCPTime: " << mctime->size() << std::endl;
        std::cout << "Size MCPStartX: " << _MCPStartX->size() << std::endl;
        std::cout << "Size MCPStartY: " << _MCPStartY->size() << std::endl;
        std::cout << "Size MCPStartZ: " << _MCPStartZ->size() << std::endl;
        std::cout << "Size truepx: " << truepx->size() << std::endl;
        std::cout << "Size truepy: " << truepy->size() << std::endl;
        std::cout << "Size truepz: " << truepz->size() << std::endl;
        std::cout << "Size MCPEndX: " << _MCPEndX->size() << std::endl;
        std::cout << "Size MCPEndY: " << _MCPEndY->size() << std::endl;
        std::cout << "Size MCPEndZ: " << _MCPEndZ->size() << std::endl;
        std::cout << "Size _MCProc: " << _MCProc->size() << std::endl;
        std::cout << "Size _MCEndProc: " << _MCEndProc->size() << std::endl;
        std::cout << "Size _angle: " << _angle->size() << std::endl;
        std::cout << "Size truep: " << truep->size() << std::endl;
        std::cout << "Size truepdg: " << truepdg->size() << std::endl;
        std::cout << "Size recopid: " << recopid->size() << std::endl;
        std::cout << "Size trkLen: " << trkLen->size() << std::endl;
        std::cout << "Size trkLenPerp: " << trkLenPerp->size() << std::endl;
        std::cout << "Size _preco: " << _preco->size() << std::endl;
        std::cout << "Size anglereco: " << anglereco->size() << std::endl;
        std::cout << "Size erecon: " << erecon->size() << std::endl;
        std::cout << "Size prob_arr: " << prob_arr->size() << std::endl;
        std::cout << "Size _nFSP: " << _nFSP->size() << " value " << _nFSP->at(0) << std::endl;
        std::cout << "Size detected: " << detected->size() << std::endl;
    	
	for(int i = 0; i < theta->size(); i++)
        {
   	
		std::cout << "value of theta is: " << theta->at(i) << std::endl;	
	
	} 
   }
}
}
}
