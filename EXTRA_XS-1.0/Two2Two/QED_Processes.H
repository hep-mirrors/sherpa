#ifndef QED_Process_H
#define QED_Process_H

#include "XS_Group.H"

namespace EXTRAXS {

  class QED_Processes: public XS_Group {
  private:

    void CreateFSRChannels();

  public:

    QED_Processes(PDF::ISR_Handler *const isr,BEAM::Beam_Spectra_Handler *const beam,
		  const ATOOLS::Flavour *fl,ATOOLS::Selector_Data *const seldata,
		  const PHASIC::scl::scheme scalescheme,const int kfactorscheme); 

  };

}

#endif
