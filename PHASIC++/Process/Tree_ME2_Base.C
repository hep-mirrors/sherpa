#include "PHASIC++/Process/Tree_ME2_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Tree_ME2_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

Tree_ME2_Base::Tree_ME2_Base(const Process_Info &pi,
                             const Flavour_Vector &flavs):
  m_pinfo(pi), m_flavs(flavs)
{
}

Tree_ME2_Base::~Tree_ME2_Base() {}

typedef ATOOLS::Getter_Function
<Tree_ME2_Base,PHASIC::Process_Info> Tree_ME2_Getter;

Tree_ME2_Base* Tree_ME2_Base::GetME2(const PHASIC::Process_Info& pi)
{
  Tree_ME2_Getter::Getter_List glist(Tree_ME2_Getter::GetGetters());
  for (Tree_ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Tree_ME2_Base *me2=(*git)->GetObject(pi);
    if (me2) return me2;
  }
  return NULL;
}

Tree_ME2_Base *Tree_ME2_Base::GetME2(const std::string& tag,
				     const Process_Info& pi)
{
  Tree_ME2_Base* me2=Tree_ME2_Getter::GetObject(tag, pi);
  if (me2==NULL) {
    THROW(fatal_error, "Did not find ME^2 "+tag);
  }
  else return me2;
}

double Trivial_Tree::Calc(const ATOOLS::Vec4D_Vector &p)
{
  return 0.0;
}
