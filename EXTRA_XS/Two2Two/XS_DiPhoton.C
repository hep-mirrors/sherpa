#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Model.H"

#include "EXTRA_XS/Main/ME2_Base.H"

#define PropID(i,j) ((1<<i)|(1<<j))

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;

namespace EXTRAXS {
  class XS_PP_ffbar : public ME2_Base {
  private:
    int     m_r, m_qcd;
    double  m_cpl, m_m2;
  public:
    XS_PP_ffbar(const External_ME_Args& args);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_PP_SSbar : public ME2_Base {
  private:
    int     m_r, m_qcd;
    double  m_cpl, m_m2;
  public:
    XS_PP_SSbar(const External_ME_Args& args);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_PP_S : public ME2_Base {
  private:
    double  m_mass, m_BR, m_ME2;
  public:
    XS_PP_S(const External_ME_Args& args);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };
}

XS_PP_ffbar::XS_PP_ffbar(const External_ME_Args& args) :
  ME2_Base(args)
{
  m_sintt=2|4;
  m_r=m_flavs[2].IsAnti();
  m_qcd=m_flavs[2].StrongCharge();
  m_cpl=sqr(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_QED"))
        *sqr(sqr(m_flavs[2].Charge()))
        *(m_flavs[2].Strong()?3.0:1.0);
  m_m2=sqr(m_flavs[2].Mass());
  for (short int i=0;i<4;i++) m_colours[i][0] = m_colours[i][1] = 0;
  m_oew=2; m_oqcd=0;
  m_cfls[PropID(0,2)] = Flavour_Vector{};
  m_cfls[PropID(1,2)] = Flavour_Vector{};
  m_cfls[PropID(0,3)] = Flavour_Vector{};
  m_cfls[PropID(1,3)] = Flavour_Vector{};
  m_cfls[PropID(0,2)].push_back(m_flavs[2]);
  m_cfls[PropID(1,2)].push_back(m_flavs[2]);
  m_cfls[PropID(0,3)].push_back(m_flavs[3]);
  m_cfls[PropID(1,3)].push_back(m_flavs[3]);
}

double XS_PP_ffbar::operator()(const Vec4D_Vector& mom)
{
  // checked with Budnev, Ginzburg, Meledin, Serbo, (E.5)
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2+m_r]).Abs2();
  double u=(mom[0]-mom[3-m_r]).Abs2();
  //if (s<m_threshold) return 0.;
  double tp(t-m_m2), up(u-m_m2);
  double mtt(2.0*(tp*(up-2.0*m_m2)-4.0*m_m2*m_m2)/(tp*tp));
  double muu(2.0*(up*(tp-2.0*m_m2)-4.0*m_m2*m_m2)/(up*up));
  double mtu(2.0*m_m2*(s-4.0*m_m2)/(tp*up));
  return m_cpl*CouplingFactor(0,2)*(mtt+muu+2.0*mtu);
}


bool XS_PP_ffbar::SetColours(const Vec4D_Vector& mom)
{
  size_t nc(m_qcd?Flow::Counter():0);
  m_colours[2+m_r][0]=m_colours[3-m_r][1]=nc;
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_PP_ffbar,"XS_PP_ffbar")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_PP_ffbar>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsPhoton() && fl[1].IsPhoton() &&
      fl[2].Charge() &&	fl[3]==fl[2].Bar() &&
      args.m_orders[0]==0 && args.m_orders[1]==2) {
    if (fl[2].IsFermion()) return new XS_PP_ffbar(args);
  }
  return NULL;
}

XS_PP_SSbar::XS_PP_SSbar(const External_ME_Args& args) :
  ME2_Base(args)
{
  m_sintt     = 2|4;
  m_r         = m_flavs[2].IsAnti();
  m_qcd       = m_flavs[2].StrongCharge();
  m_cpl       = ( sqr(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_QED")) *
		  sqr(sqr(m_flavs[2].Charge())) *
		  (m_flavs[2].Strong()?3.0:1.0) );
  m_m2        = sqr(m_flavs[2].Mass());
  for (short int i=0;i<4;i++) m_colours[i][0] = m_colours[i][1] = 0;
  m_oew=2; m_oqcd=0;
  m_cfls[PropID(0,2)] = Flavour_Vector{};
  m_cfls[PropID(1,2)] = Flavour_Vector{};
  m_cfls[PropID(0,3)] = Flavour_Vector{};
  m_cfls[PropID(1,3)] = Flavour_Vector{};
  m_cfls[PropID(0,2)].push_back(m_flavs[2]);
  m_cfls[PropID(1,2)].push_back(m_flavs[2]);
  m_cfls[PropID(0,3)].push_back(m_flavs[3]);
  m_cfls[PropID(1,3)].push_back(m_flavs[3]);
}

double XS_PP_SSbar::operator()(const Vec4D_Vector& mom)
{
  // from Budnev, Ginzburg, Meledin, Serbo, (E.6)
  double s   = (mom[0]+mom[1]).Abs2();
  if (s<4.*m_m2) return 0.;
  double tp  = (mom[0]-mom[2+m_r]).Abs2()-m_m2;
  double up  = (mom[0]-mom[3-m_r]).Abs2()-m_m2;
  double r2  = tp*up/s-m_m2;
  double ME2 = 2.*(sqr(m_m2)+sqr(r2))/(tp*up);
  return m_cpl*CouplingFactor(0,2)*ME2;
}

bool XS_PP_SSbar::SetColours(const Vec4D_Vector& mom)
{
  size_t nc(m_qcd?Flow::Counter():0);
  m_colours[2+m_r][0]=m_colours[3-m_r][1]=nc;
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_PP_SSbar,"XS_PP_SSbar")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_PP_SSbar>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsPhoton() && fl[1].IsPhoton() &&
      fl[2].Charge() &&	fl[3]==fl[2].Bar() &&
      args.m_orders[0]==0 && args.m_orders[1]==2) {
    if (fl[2].IsScalar())  return new XS_PP_SSbar(args);
  }
  return NULL;
}




XS_PP_S::XS_PP_S(const External_ME_Args& args) : ME2_Base(args)
{
  m_mass  = m_flavs[2].Mass();
  m_BR    = 0.;
  switch (m_flavs[2].Kfcode()) {
  case kf_h0:            m_BR = 0.0025;   break;
  case kf_pi:            m_BR = 0.98823;  break;
  case kf_eta:           m_BR = 0.3936;   break;
  case kf_eta_prime_958: m_BR = 0.02307;  break; 
  case kf_eta_c_1S:      m_BR = 0.000166; break; 
  case kf_eta_b:         m_BR = 0.000;    break; 
  default:
    THROW(fatal_error,"Scalar particle not found in list of photonic BRs: "+m_flavs[2].IDName());
  }
  m_ME2   = 8.*M_PI/m_mass*m_BR; 
  for (short int i=0;i<3;i++) m_colours[i][0] = m_colours[i][1] = 0;
  m_oew=2; m_oqcd=0;
}

double XS_PP_S::operator()(const Vec4D_Vector& mom) { return m_ME2; }

bool XS_PP_S::SetColours(const Vec4D_Vector& mom) { return true; }

DECLARE_TREEME2_GETTER(EXTRAXS::XS_PP_S,"XS_PP_S")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_PP_S>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=3) return NULL;
  if (fl[0].IsPhoton() && fl[1].IsPhoton() && fl[2].Charge()==0 &&
      fl[2].IsScalar() &&
      args.m_orders[0]==0 && args.m_orders[1]==2) {
    if (fl[2]==Flavour(kf_h0) ||
	fl[2]==Flavour(kf_pi) || fl[2]==Flavour(kf_eta) || fl[2]==Flavour(kf_eta_prime_958) ||
	fl[2]==Flavour(kf_eta_c_1S) || fl[2]==Flavour(kf_eta_b) )
      return new XS_PP_S(args);
  }
  return NULL;
}

