#include "PHASIC++/EWSudakov/Sudakov.H"
#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Sudakov::Sudakov(Process_Base *const proc)
{
  p_proc=proc;
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetNIn(p_proc->NIn());
  p_ampl->SetProc(p_proc);
  for(int i(0);i<p_proc->NIn()+p_proc->NOut();++i)
    if (i<p_proc->NIn()) p_ampl->CreateLeg(Vec4D(),p_proc->Flavours()[i].Bar());
    else p_ampl->CreateLeg(Vec4D(),p_proc->Flavours()[i]);
  p_ci  = new Comix_Interface(p_proc,p_ampl);
}

Sudakov::~Sudakov()
{
  if(p_ci) delete p_ci;
}

double Sudakov::EWSudakov(ATOOLS::Vec4D_Vector &mom)
{
  // do something with p_ci
  return 0.;
}
