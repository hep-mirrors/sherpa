#ifndef COMIX_Main_Process_Group_H
#define COMIX_Main_Process_Group_H

#include "COMIX/Main/Process_Base.H"
#include "PHASIC++/Process/Process_Group.H"

namespace COMIX {

  class Single_Process;

  class Process_Group: 
    public PHASIC::Process_Group, public Process_Base {
  protected:

    void InitPSGenerator(const size_t &ismode);

  public :

    // constructors
    Process_Group(); 
    Process_Group(Model *const model);

    // member functions
    bool Initialize(std::map<std::string,std::string> *const pmap,
		    std::vector<Single_Process*> *const procs);
    bool MapProcess();

    PHASIC::Process_Base *GetProcess(const PHASIC::Process_Info &pi) const;
    bool Initialize(PHASIC::Process_Base *const proc);

    void ConstructPSVertices(PS_Generator *ps);
    
    bool GeneratePoint();

    bool Tests();

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

  };// end of class Process_Group

}// end of namespace COMIX

#endif


