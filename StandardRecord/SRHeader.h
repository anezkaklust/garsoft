#ifndef CAFSRHEADER_H
#define CAFSRHEADER_H

namespace caf
{
  class SRHeader
  {
  public:
    SRHeader();

    // global event info
    int   event;      ///< number of the event being processed
    int   run;        ///< number of the run being processed
    int   subrun;     ///< number of the sub-run being processed

    float TPC_X;      ///< center of TPC stored as per-event & compressed by root
    float TPC_Y;
    float TPC_Z;
  };
}

#endif
