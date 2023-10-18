CXX       = g++
HEADFLAGS = -I$(MRB_SOURCE)/garsoft -I$(NUSIMDATA_DIR)/source/ 
ROOTFLAGS = `root-config --cflags --glibs`
DEBFLAGS = -g -Wall -fPIC
OPTFLAGS = -O3
PROFLAGS = -pg -O3



FRAMEWORK = AnaTreeFrame.cxx

CODE_1 = makeDST



all:
	(make makeDST_OPT)



makeDST_DEB: $(CODE_1).cxx $(FRAMEWORK)
	$(CXX) -c $(CODE_1).cxx $(ROOTFLAGS) $(DEBFLAGS) $(HEADFLAGS)
	$(CXX) -o $(CODE_1).exe $(FRAMEWORK) $(CODE_1).o $(ROOTFLAGS) $(DEBFLAGS) $(HEADFLAGS)
	rm $(CODE_1).o
makeDST_OPT: $(CODE_1).cxx $(FRAMEWORK)
	$(CXX) -c $(CODE_1).cxx $(ROOTFLAGS) $(OPTFLAGS) $(HEADFLAGS)
	$(CXX) -o $(CODE_1).exe $(FRAMEWORK) $(CODE_1).o $(ROOTFLAGS) $(OPTFLAGS) $(HEADFLAGS)
	rm $(CODE_1).o
makeDST_PRO: $(CODE_1).cxx $(FRAMEWORK)
	$(CXX) -c $(CODE_1).cxx                          $(ROOTFLAGS) $(PROFLAGS) $(HEADFLAGS)
	$(CXX) -o $(CODE_1).exe $(FRAMEWORK) $(CODE_1).o $(ROOTFLAGS) $(PROFLAGS) $(HEADFLAGS)
	rm $(CODE_1).o
