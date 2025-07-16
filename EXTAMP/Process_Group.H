#ifndef EXTAMP_Main_Process_Group_H
#define EXTAMP_Main_Process_Group_H

#include "PHASIC++/Process/Process_Group.H"

namespace EXTAMP {

  class Process_Group: public PHASIC::Process_Group {


    /* Pure virtuals from PHASIC::Process_Group */

    /* Instantiate a PARTONIC process (i.e. pi specifies a
       'Single_Process'). This being a method of a 'Process_Group'
       does not make any sense, but need comply with legacy structure.
       Just call the appropriate method of the external ME
       interface. */
    
    PHASIC::Process_Base *GetProcess(const PHASIC::Process_Info &pi) const {
      /* This returns NULL if the partonic channel does not exist, as
	 required by PHASIC::Process_Group. */

      PHASIC::Process_Info cpi(pi);
      if(!External_ME_Interface::PartonicProcessExists(cpi)) return NULL;
      return External_ME_Interface::InstantiatePartonicProcess(cpi);
    }

    bool Initialize(Process_Base *const proc)
    {
      proc->SetParent((PHASIC::Process_Base*)this);
      return true;
    };

  };

  
}
#endif
