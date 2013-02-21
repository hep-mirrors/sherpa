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
  Tree_ME2_Base(pi,flavs), m_oew(99), m_oqcd(99), m_sintt(7)
{
  m_symfac=pi.m_fi.FSSymmetryFactor();
  m_symfac*=pi.m_ii.ISSymmetryFactor();
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

double ME2_Base::Calc(const ATOOLS::Vec4D_Vector &p)
{
  return (*this)(p)*m_symfac;
}

bool ME2_Base::SetColours(const Vec4D_Vector& mom)
{
//   THROW(fatal_error, "Virtual function called.");
  return false;
}

double ME2_Base::CouplingFactor(const int oqcd,const int oew) const
{
  double fac(1.0);
  if (p_aqcd && oqcd) fac*=pow(p_aqcd->Factor(),oqcd);
  if (p_aqed && oew) fac*=pow(p_aqed->Factor(),oew);
  return fac;
}
