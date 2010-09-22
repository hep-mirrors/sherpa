#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Virtual_ME2_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

Virtual_ME2_Base::Virtual_ME2_Base(const Process_Info& pi,
                             const Flavour_Vector& flavs) :
  m_pinfo(pi), m_flavs(flavs),
  m_res(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  m_mur2(1.0), m_mode(0), m_drmode(0)
{
}

Virtual_ME2_Base::~Virtual_ME2_Base()
{
}

double Virtual_ME2_Base::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  //MSbar
   return 2.*M_PI*m_mur2/(mom[0]*mom[1]);
}

bool Virtual_ME2_Base::SetColours(const Vec4D_Vector& mom)
{
  THROW(fatal_error, "Virtual function called.");
  return false;
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
