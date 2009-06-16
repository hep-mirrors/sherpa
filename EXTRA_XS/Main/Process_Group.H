#ifndef EXTRA_XS_Main_Process_Group_H
#define EXTRA_XS_Main_Process_Group_H

#include "PHASIC++/Process/Process_Group.H"

namespace EXTRAXS {

  class Process_Group: public PHASIC::Process_Group {
  public :

    // constructors
    Process_Group();
    ~Process_Group();

    // member functions
    PHASIC::Process_Base *GetProcess(const PHASIC::Process_Info &pi) const;
    bool Initialize(PHASIC::Process_Base *const proc);
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
  };

}// end of namespace EXTRAXS

#endif


