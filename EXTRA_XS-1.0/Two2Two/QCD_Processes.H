#ifndef QCD_Process_H
#define QCD_Process_H

#include "XS_Group.H"

namespace EXTRAXS {

  class QCD_Processes: public XS_Group {
    double aS;
    void Fill4qmodes();
    void Fill2q2gmodes();
    void Fill4gmodes();
    void CreateFSRChannels();
  public:
    QCD_Processes(PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		  ATOOLS::Flavour *,ATOOLS::Selector_Data *,int,int,double);
  };
}

#endif
