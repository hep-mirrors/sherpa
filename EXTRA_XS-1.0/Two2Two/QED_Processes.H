#ifndef QED_Process_H
#define QED_Process_H

#include "XS_Group.H"

namespace AMEGIC { class  Process_Group; }

namespace EXTRAXS {

  class QED_Processes: public XS_Group {
  public:
    QED_Processes();
    QED_Processes(int);
    void Initialize(ISR::ISR_Handler *, BEAM::Beam_Handler *,
		    APHYTOOLS::Selector_Data *, AMEGIC::Process_Group *);
  };
}

#endif
