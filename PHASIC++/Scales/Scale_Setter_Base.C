#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Scale_Setter_Base
#define PARAMETER_TYPE PHASIC::Scale_Setter_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

#include "ATOOLS/Org/Smart_Pointer.C"

namespace ATOOLS { template class SP(Scale_Setter_Base); }

std::ostream &operator<<(std::ostream &ostr,const stp::id &scl)
{
  switch (scl) {
  case stp::ren: return ostr<<"ren";
  case stp::fac: return ostr<<"fac";
  case stp::size: return ostr<<"<error>";
  }
  return ostr<<"<unknown>";
}

Scale_Setter_Base::~Scale_Setter_Base()
{
}

double Scale_Setter_Base::CalculateScale2(const std::vector<ATOOLS::Vec4D> &p)
{
  return m_scale[stp::fac];
}

void Scale_Setter_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available scale choices\n\n";
  Scale_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

Vec4D Scale_Setter_Base::Momentum(const size_t &i) const
{
  if (i>p_proc->NIn()+p_proc->NOut()) THROW(fatal_error,"Momentum index too large");
  return p_proc->Integrator()->PSHandler()->LabPoint()[i];
}

double Scale_Setter_Base::HT() const
{
  SP(Phase_Space_Handler) psh(p_proc->Integrator()->PSHandler());
  double ht(0.0);
  for (size_t i(0);i<p_proc->NOut();++i) 
    ht+=psh->LabPoint()[i+p_proc->NIn()].PPerp();
  return sqrt(ht);
}
