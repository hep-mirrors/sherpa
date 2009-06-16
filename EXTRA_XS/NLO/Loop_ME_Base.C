#include "EXTRA_XS/NLO/Loop_ME_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE EXTRAXS::Loop_ME_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HELICITIES;

Loop_ME_Base::Loop_ME_Base(const Process_Info& pi,
                           const Flavour_Vector& flavs) :
  m_pinfo(pi),
  m_flavs(flavs),
  m_res(flavs, DivArrC(vector<Complex>(6, Complex(0.0,0.0))))
{
}

Loop_ME_Base::~Loop_ME_Base()
{
}

typedef Getter_Function<Loop_ME_Base, Process_Info> Loop_ME_Getter;

Loop_ME_Base* Loop_ME_Base::GetME(const PHASIC::Process_Info& pi)
{
  Loop_ME_Getter::Getter_List glist(Loop_ME_Getter::GetGetters());
  for (Loop_ME_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Loop_ME_Base* me=(*git)->GetObject(pi);
    if (me) return me;
  }
  return NULL;
}

Loop_ME_Base* Loop_ME_Base::GetME(const std::string& tag,const Process_Info& pi)
{
  Loop_ME_Base* me=Loop_ME_Getter::GetObject(tag, pi);
  if (me==NULL) {
    THROW(fatal_error, "Did not find ME "+tag);
  }
  else return me;
}
