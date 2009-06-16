#include "PHASIC++/Process/ME_Generator_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::ME_Generator_Base
#define PARAMETER_TYPE PHASIC::ME_Generator_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;

ME_Generator_Base::~ME_Generator_Base()
{
}

double ME_Generator_Base::Mass(const ATOOLS::Flavour &fl) const
{
  return fl.Mass();
}

void ME_Generator_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  ME_Generator_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}
