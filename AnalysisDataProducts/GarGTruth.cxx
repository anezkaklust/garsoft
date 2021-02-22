
#include "GarGTruth.h"

using namespace gar::adp;

GarGTruth::GarGTruth() {}

GarGTruth::GarGTruth(const simb::GTruth& gt) {

    FillGT(gt);

}

garana::GTruth GarGTruth::GetGaranaObj(){
    return fGT;
}

void GarGTruth::FillGT(const simb::GTruth& gt) {

    fGT.fVertex = gt.fVertex;
    fGT.fweight = gt.fweight;
    fGT.fprobability = gt.fprobability;
    fGT.fXsec = gt.fXsec;
    fGT.fDiffXsec = gt.fDiffXsec;
    fGT.fGPhaseSpace = gt.fGPhaseSpace;
 
    fGT.fProbePDG = gt.fProbePDG;
    fGT.fProbeP4 = gt.fProbeP4;
    fGT.fTgtP4 = gt.fTgtP4;

    fGT.ftgtZ = gt.ftgtZ;
    fGT.ftgtA = gt.ftgtA;
    fGT.ftgtPDG = gt.ftgtPDG;   
    fGT.fHitNucPDG = gt.fHitNucPDG;
    fGT.fHitQrkPDG = gt.fHitQrkPDG;
    fGT.fIsSeaQuark = gt.fIsSeaQuark;
    fGT.fHitNucP4 = gt.fHitNucP4;
    fGT.fHitNucPos = gt.fHitNucPos;

    fGT.fGscatter = gt.fGscatter; 
    fGT.fGint = gt.fGint;

    fGT.fgQ2 = gt.fgQ2;
    fGT.fgq2 = gt.fgq2;
    fGT.fgW = gt.fgW ;
    fGT.fgT = gt.fgT;
    fGT.fgX = gt.fgX;
    fGT.fgY = gt.fgY;
    fGT.fFSleptonP4 = gt.fFSleptonP4; 
    fGT.fFShadSystP4 = gt.fFShadSystP4; 

    fGT.fIsCharm = gt.fIsCharm;
    fGT.fCharmHadronPdg = gt.fCharmHadronPdg; 
    fGT.fIsStrange = gt.fIsStrange;
    fGT.fStrangeHadronPdg = gt.fStrangeHadronPdg; 
    fGT.fNumProton = gt.fNumProton;
    fGT.fNumNeutron = gt.fNumNeutron;
    fGT.fNumPi0 = gt.fNumPi0;
    fGT.fNumPiPlus = gt.fNumPiPlus;  
    fGT.fNumPiMinus = gt.fNumPiMinus;
    fGT.fResNum = gt.fResNum;
    fGT.fDecayMode = gt.fDecayMode;
}
