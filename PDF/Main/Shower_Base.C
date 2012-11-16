#include "PDF/Main/Shower_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::Shower_Base
#define PARAMETER_TYPE PDF::Shower_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

#include "PDF/Main/Cluster_Definitions_Base.H"

using namespace PDF;
using namespace ATOOLS;

Shower_Base::Shower_Base(const std::string &name):
  p_cluster(NULL), m_name(name), m_weight(1.0), m_on(1) {}

Shower_Base::~Shower_Base() 
{
  if (p_cluster) delete p_cluster;
}

void Shower_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  Shower_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}

