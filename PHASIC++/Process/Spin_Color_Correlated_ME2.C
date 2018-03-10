#include "PHASIC++/Process/Spin_Color_Correlated_ME2.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Spin_Color_Correlated_ME2
#define PARAMETER_TYPE PHASIC::Correlator_Args
#include "ATOOLS/Org/Getter_Function.C"


using namespace PHASIC;
using namespace ATOOLS;


Spin_Color_Correlated_ME2::Spin_Color_Correlated_ME2(const Correlator_Args& args):
  m_flavs(args.m_flavs), p_aqcd(NULL), p_aqed(NULL) { }


void Spin_Color_Correlated_ME2::SetCouplings(MODEL::Coupling_Data* p_rqcd,
					     MODEL::Coupling_Data* p_rqed)
{
  p_aqcd=p_rqcd;
  p_aqed=p_rqed;
}


double Spin_Color_Correlated_ME2::AlphaQCD() const
{
  return p_aqcd ? p_aqcd->Default()*p_aqcd->Factor() : MODEL::s_model->ScalarConstant("alpha_S");
}


double Spin_Color_Correlated_ME2::AlphaQED() const
{
  return p_aqed ? p_aqed->Default()*p_aqed->Factor() : MODEL::s_model->ScalarConstant("alpha_QED");
}


typedef Getter_Function
<Spin_Color_Correlated_ME2,Correlator_Args> Spin_Color_Correlated_ME2_Getter;

Spin_Color_Correlated_ME2* Spin_Color_Correlated_ME2::GetME2(const Correlator_Args& args)
{
  Spin_Color_Correlated_ME2_Getter::Getter_List glist(Spin_Color_Correlated_ME2_Getter::GetGetters());
  for (Spin_Color_Correlated_ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Spin_Color_Correlated_ME2 *me2=(*git)->GetObject(args);
    if (me2) return me2;
  }
  return NULL;
}
