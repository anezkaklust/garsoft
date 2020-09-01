#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"

#include "Utils.h"

class CAF {

    typedef std::pair<float, std::string> P;

public:

    /* Default Constructor */
    CAF();

    /* Constructor */
    CAF( std::string infile, std::string filename, int correct4origin, double *originTPC);

    /* Default Copy Constructor */
    CAF(const CAF &) = default;

    /* Destructor */
    ~CAF();

    /* Method to book the output TFile */
    bool BookTFile();

    /* Method to fill the output TTree */
    void FillTTree();

    /* Method to write the output TTree */
    void WriteTTree();

    /* Method to close the output TTree */
    void CloseTFile();

    /* Method to clear the event variables */
    void ClearVectors();

    /* Method to loop */
    void loop();

    /* Method to check vector size */
    bool CheckVectorSize();

private:

    TFile * cafFile; ///< The output TFile pointer */
    TTree * cafMVA; ///< The output TTree pointer */

    TFile* _intfile;
    TTree* _inttree;

    Utils *_util;

    std::string _inputfile;
    std::string _outputFile; ///< The output TFile name */
    unsigned int _correct4origin;     ///< sets the string for the coordinates origins (World or TPC)

    //Event-wise values
    unsigned int _Run, _Event, _SubRun;
    //Generator values
    std::vector<int> mode, ccnc, ntype, gint, weight, tgtpdg, gt_t, intert, detected;
    std::vector<double> q2, w, y, x, theta, t, mctime, mcnupx, mcnupy, mcnupz, vertx, verty, vertz;
    //MC Particle Values, with motherid added
    std::vector<unsigned int> _nFSP;
    std::vector<int> mctrkid, motherid, pdgmother, truepdg, _MCPStartX, _MCPStartY, _MCPStartZ, _MCPEndX, _MCPEndY, _MCPEndZ;
    std::vector<std::string> _MCProc, _MCEndProc;
    std::vector<double> trkLen, trkLenPerp, truep, truepx, truepy, truepz, _angle;
    //Reco values
    std::vector<int> recopid, recopidecal;
    std::vector<double> prob_arr, partereco, anglereco, _preco, erecon, etime;
    //Geometry
    std::vector<unsigned int> isFidStart, isTPCStart, isCaloStart, isInBetweenStart, isThroughCaloStart;
    std::vector<unsigned int> isFidEnd, isTPCEnd, isCaloEnd, isInBetweenEnd, isThroughCaloEnd;
    std::vector<unsigned int> isBarrelStart, isEndcapStart, isBarrelEnd, isEndcapEnd;
};

#endif
