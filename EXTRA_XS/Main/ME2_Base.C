#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE EXTRAXS::ME2_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

ME2_Base::ME2_Base(const Process_Info& pi, const Flavour_Vector& flavs) : 
  m_pinfo(pi), m_flavs(flavs), m_oew(99), m_oqcd(99), m_sintt(7)
{
  p_colours = new int*[m_flavs.size()];
  for (size_t i(0);i<m_flavs.size();++i) {
    p_colours[i] = new int[2];
    p_colours[i][0]=p_colours[i][1]=0;
  }
}

ME2_Base::~ME2_Base()
{
  for (size_t i(0);i<4;++i) delete [] p_colours[i];
  delete [] p_colours;
}

bool ME2_Base::SetColours(const Vec4D_Vector& mom)
{
//   THROW(fatal_error, "Virtual function called.");
  return false;
}



typedef ATOOLS::Getter_Function<ME2_Base, PHASIC::Process_Info> ME2_Getter;

ME2_Base* ME2_Base::GetME2(const PHASIC::Process_Info& pi)
{
  ME2_Getter::Getter_List glist(ME2_Getter::GetGetters());
  for (ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    ME2_Base* me2=(*git)->GetObject(pi);
    if (me2) return me2;
  }
  return NULL;
}

ME2_Base* ME2_Base::GetME2(const std::string& tag, const Process_Info& pi)
{
  ME2_Base* me2=ME2_Getter::GetObject(tag, pi);
  if (me2==NULL) {
    THROW(fatal_error, "Did not find ME^2 "+tag);
  }
  else return me2;
}
