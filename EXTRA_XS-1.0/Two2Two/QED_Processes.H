#ifndef QED_Process_H
#define QED_Process_H

#include "Multi_Channel.H"
#include "XS_Group.H"

namespace EXTRAXS {

  class QED_Processes: public XS_Group {
    void CreateFSRChannels();
  public:
    QED_Processes(PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		  ATOOLS::Flavour *,ATOOLS::Selector_Data *,int,int,double);
  };
}

#endif
