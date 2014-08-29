#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Virtual_ME2_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

Virtual_ME2_Base::Virtual_ME2_Base(const Process_Info& pi,
                             const Flavour_Vector& flavs) :
  m_pinfo(pi), m_flavs(flavs),
  m_res(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  m_mur2(1.0), m_mode(0), m_drmode(0), m_colmode(0),
  m_born(0.0), m_norm(1.0)
{
}

Virtual_ME2_Base::~Virtual_ME2_Base()
{
}

double Virtual_ME2_Base::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  //MSbar
  return 4.0*M_PI;
}

double Virtual_ME2_Base::ScaleDependenceCoefficient(const int i)
{
  THROW(not_implemented,"Invalid call");
  return 0.0;
}

void Virtual_ME2_Base::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  p_aqcd=p_aqed=NULL;
  if (cpls.find("Alpha_QCD")!=cpls.end()) p_aqcd=cpls.Get("Alpha_QCD");
  if (cpls.find("Alpha_QED")!=cpls.end()) p_aqed=cpls.Get("Alpha_QED");
}

double Virtual_ME2_Base::AlphaQCD() const
{
  return p_aqcd ? p_aqcd->Default()*p_aqcd->Factor() : s_model->ScalarFunction("alpha_S");
}

double Virtual_ME2_Base::AlphaQED() const
{
  return p_aqed ? p_aqed->Default()*p_aqed->Factor() : s_model->ScalarFunction("alpha_QED");
}

typedef ATOOLS::Getter_Function<Virtual_ME2_Base, PHASIC::Process_Info>
Virtual_ME2_Getter;

Virtual_ME2_Base* Virtual_ME2_Base::GetME2(const PHASIC::Process_Info& pi)
{
  Virtual_ME2_Getter::Getter_List glist(Virtual_ME2_Getter::GetGetters());
  for (Virtual_ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Virtual_ME2_Base* me2=(*git)->GetObject(pi);
    if (me2) return me2;
  }
  return NULL;
}

Virtual_ME2_Base* Virtual_ME2_Base::GetME2(const std::string& tag,
                                           const Process_Info& pi)
{
  Virtual_ME2_Base* me2=Virtual_ME2_Getter::GetObject(tag, pi);
  if (me2==NULL) {
    THROW(fatal_error, "Did not find ME^2 "+tag);
  }
  else return me2;
}
