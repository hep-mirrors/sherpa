#include "PHASIC++/Process/Color_Correlated_ME2.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Color_Correlated_ME2
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"


using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;


Color_Correlated_ME2::Color_Correlated_ME2(const Process_Info &pi,
                             const Flavour_Vector &flavs):
  m_pinfo(pi), m_flavs(flavs), p_aqcd(NULL), p_aqed(NULL) { }


void Color_Correlated_ME2::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  p_aqcd=cpls.Get("Alpha_QCD");
  p_aqed=cpls.Get("Alpha_QED");
}


double Color_Correlated_ME2::AlphaQCD() const
{
  return p_aqcd ? p_aqcd->Default()*p_aqcd->Factor() : s_model->ScalarConstant("alpha_S");
}


double Color_Correlated_ME2::AlphaQED() const
{
  return p_aqed ? p_aqed->Default()*p_aqed->Factor() : s_model->ScalarConstant("alpha_QED");
}


typedef ATOOLS::Getter_Function
<Color_Correlated_ME2,PHASIC::Process_Info> Color_Correlated_ME2_Getter;

Color_Correlated_ME2* Color_Correlated_ME2::GetME2(const PHASIC::Process_Info& pi)
{
  Color_Correlated_ME2_Getter::Getter_List glist(Color_Correlated_ME2_Getter::GetGetters());
  for (Color_Correlated_ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Color_Correlated_ME2 *me2=(*git)->GetObject(pi);
    if (me2) return me2;
  }
  return NULL;
}


Color_Correlated_ME2 *Color_Correlated_ME2::GetME2(const std::string& tag,
						   const Process_Info& pi)
{
  Color_Correlated_ME2* me2=Color_Correlated_ME2_Getter::GetObject(tag, pi);
  if (me2==NULL) THROW(fatal_error, "Did not find correlated ME "+tag);
  return me2;
}
