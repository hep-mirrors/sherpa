#include "PHASIC++/Process/Tree_ME2_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Tree_ME2_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

Tree_ME2_Base::Tree_ME2_Base(const Process_Info &pi,
                             const Flavour_Vector &flavs):
  m_pinfo(pi), m_flavs(flavs), p_aqcd(NULL), p_aqed(NULL)
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

void Tree_ME2_Base::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  p_aqcd=cpls->Get("Alpha_QCD");
  p_aqed=cpls->Get("Alpha_QED");
}

double Tree_ME2_Base::AlphaQCD() const
{
  return p_aqcd ? p_aqcd->Default()*p_aqcd->Factor() : s_model->ScalarFunction("alpha_S");
}

double Tree_ME2_Base::AlphaQED() const
{
  return p_aqed ? p_aqed->Default()*p_aqed->Factor() : s_model->ScalarFunction("alpha_QED");
}




double Trivial_Tree::Calc(const ATOOLS::Vec4D_Vector &p)
{
  return 0.0;
}
