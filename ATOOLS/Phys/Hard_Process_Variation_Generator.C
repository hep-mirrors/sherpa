#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "PHASIC++/Process/Process_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Hard_Process_Variation_Generator_Base
#define PARAMETER_TYPE ATOOLS::Hard_Process_Variation_Generator_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace ATOOLS;

void Hard_Process_Variation_Generator_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available hard-process variation generators\n\n";
  Hard_Process_Variation_Generator_Getter_Function::
    PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}
