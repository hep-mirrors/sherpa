#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Process/External_ME_Args.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE EXTRAXS::ME2_Base
#define PARAMETER_TYPE PHASIC::External_ME_Args
#include "ATOOLS/Org/Getter_Function.C"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

ME2_Base::ME2_Base(const External_ME_Args& args) : 
  Tree_ME2_Base(args), m_oew(99), m_oqcd(99), m_sintt(7),
  m_sprimemin(-1.), m_sprimemax(-1.), 
  m_hasinternalscale(false), m_internalscale(sqr(rpa->gen.Ecms()))
{
  m_symfac= Flavour::FSSymmetryFactor(args.m_outflavs);
  m_symfac*=Flavour::ISSymmetryFactor(args.m_inflavs);
  m_colours.resize(m_flavs.size());
  for (size_t i(0);i<m_flavs.size();++i) {
    m_colours[i].resize(2);
    m_colours[i][0]=m_colours[i][1]=0;
  }
}

ME2_Base::~ME2_Base() { }

double ME2_Base::Calc(const ATOOLS::Vec4D_Vector &p)
{
  // the symfac is multiplied here to cancel the symfac in the ME2's
  // since the Tree_ME2_Base::Calc function is supposed to return
  // the pure ME2, without sym fac
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

int ME2_Base::OrderQCD(const int &id) const
{
  return m_oqcd;
}

int ME2_Base::OrderEW(const int &id) const
{
  return m_oew;
}

