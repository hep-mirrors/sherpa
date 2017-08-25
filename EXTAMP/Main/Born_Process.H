#ifndef EXTAMP_Born_Process_H
#define EXTAMP_Born_Process_H

#include "EXTAMP/Main/Process.H"

namespace PHASIC {
  class Tree_ME2_Base;
  class Process_Info;
}

namespace EXTAMP {

  class Born_Process : public Process {

  public :

    Born_Process(const PHASIC::Process_Info& pi);
    
    ~Born_Process() {};

    double Partonic(const ATOOLS::Vec4D_Vector &p,
		    const int mode);

  private:

    PHASIC::Tree_ME2_Base* p_born_me;

  };

}

#endif
