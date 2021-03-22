
#include "AnalysisDataProducts/AnaTrack.h"
#include <iostream>

using namespace gar::adp;

garana::Track AnaTrack::GaranaTrack() {

    const float* tmp = fTrack.Vertex();
    TLorentzVector vtx(tmp[0], tmp[1], tmp[2], fTrack.Time());

    tmp = fTrack.End();
    TLorentzVector end(tmp[0], tmp[1], tmp[2], fTrack.Time()); // eventually track end times should be different

    tmp = fTrack.VtxDir();      
    TVector3 vtxDir(tmp[0], tmp[1], tmp[2]);

    tmp = fTrack.EndDir();
    TVector3 endDir(tmp[0], tmp[1], tmp[2]);

    //const int chgbeg = fTrack.ChargeBeg();
    //const int chgend = fTrack.ChargeEnd();

    garana::Track trk(fTrack.LengthForward(), fTrack.LengthBackward(), fTrack.Momentum_beg(), fTrack.Momentum_end(),
			vtx, end, vtxDir, endDir, fTrack.ChisqForward(), fTrack.ChisqBackward(), fTrack.NHits(),
			fTrack.TrackParBeg(), fTrack.TrackParEnd(), fTrack.CovMatBegPacked(), fTrack.CovMatEndPacked(),
                        /*fTrack.ChargeBeg(), fTrack.ChargeEnd(),*/ fChgF, fChgB, fPidF, fPidB, fIonF, fIonB );


    return trk;

}
