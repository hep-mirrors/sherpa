#include "PDF/Main/POWHEG_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::POWHEG_Base
#define PARAMETER_TYPE PDF::POWHEG_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PDF;
using namespace ATOOLS;

POWHEG_Base::POWHEG_Base(const std::string &name):
  m_name(name), p_shower(NULL), m_kt2min(-1.0) {}

POWHEG_Base::~POWHEG_Base() 
{
}

void POWHEG_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  POWHEG_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}

