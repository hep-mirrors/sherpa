#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::KFactor_Setter_Base
#define PARAMETER_TYPE PHASIC::KFactor_Setter_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

#include "ATOOLS/Org/Smart_Pointer.C"

namespace ATOOLS { template class SP(KFactor_Setter_Base); }

KFactor_Setter_Base::KFactor_Setter_Base(Process_Base *const proc): 
  p_proc(proc), m_on(true) {}

KFactor_Setter_Base::~KFactor_Setter_Base()
{
}

double KFactor_Setter_Base::KFactor2(const bool diff)
{
  if (diff) {
    m_kfkey<<ATOOLS::UNDEFINED_WEIGHT;
    return KFactor();
  }
  return m_kfkey.Weight();
}

void KFactor_Setter_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available kfactor choices\n\n";
  KFactor_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}
