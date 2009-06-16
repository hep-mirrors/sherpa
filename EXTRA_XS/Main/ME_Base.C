#include "EXTRA_XS/Main/ME_Base.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE EXTRAXS::ME_Base
#define PARAMETER_TYPE PHASIC::Process_Info
#include "ATOOLS/Org/Getter_Function.C"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

ME_Base::ME_Base(const Process_Info& pi, const Flavour_Vector& flavs, const Particle_Vector& parts) : 
  m_pinfo(pi), m_flavs(flavs), m_res(flavs, Complex(0.0,0.0)), m_colres(parts) 
{

}

ME_Base::~ME_Base()
{
}


typedef ATOOLS::Getter_Function<ME_Base, PHASIC::Process_Info> ME_Getter;

ME_Base* ME_Base::GetME(const PHASIC::Process_Info& pi)
{
  ME_Getter::Getter_List glist(ME_Getter::GetGetters());
  for (ME_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    ME_Base* me=(*git)->GetObject(pi);
    if (me) return me;
  }
  return NULL;
}

ME_Base* ME_Base::GetME(const std::string& tag, const Process_Info& pi)
{
  ME_Base* me=ME_Getter::GetObject(tag, pi);
  if (me==NULL) {
    THROW(fatal_error, "Did not find ME "+tag);
  }
  else return me;
}
