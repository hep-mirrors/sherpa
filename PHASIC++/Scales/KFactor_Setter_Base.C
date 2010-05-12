#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::KFactor_Setter_Base
#define PARAMETER_TYPE PHASIC::KFactor_Setter_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

#include "ATOOLS/Org/Smart_Pointer.C"

namespace ATOOLS { template class SP(KFactor_Setter_Base); }

KFactor_Setter_Base::KFactor_Setter_Base
(const KFactor_Setter_Arguments &args): 
  p_proc(args.p_proc), m_weight(0.0), m_on(true) {}

KFactor_Setter_Base::~KFactor_Setter_Base()
{
}

void KFactor_Setter_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available kfactor choices\n\n";
  KFactor_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}
