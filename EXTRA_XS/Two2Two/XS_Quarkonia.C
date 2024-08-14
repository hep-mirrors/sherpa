#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/UFO/UFO_Model.H"

#include "EXTRA_XS/Main/ME2_Base.H"

#define PropID(i,j) ((1<<i)|(1<<j))

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;


namespace EXTRAXS {
  class XS_qg_q1S0 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g1S0 : public ME2_Base {
  private:
    size_t m_S, m_a;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g1S0 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g3S1 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g3S1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}


//////////////////////////////////////////////////////////////////////////
//
// Singlet production of eta states (1S0 = 441, 551)
//
//////////////////////////////////////////////////////////////////////////
DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q1S0,"1XS_qg_q1S0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q1S0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_qg_q1S0(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_qg_q1S0(args);
    }
  }
  return NULL;
}

XS_qg_q1S0::XS_qg_q1S0(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R02 = 0.;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.49;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 4.54;
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 2./9.*sqr(4.*M_PI)*m_R02/m_mass;
}

double XS_qg_q1S0::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(t-m_mass2)-2.*s*u)/(-t*sqr(t-m_mass2));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qg_q1S0::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[m_g][m_a]   = m_colours[5-m_S][m_a] = Flow::Counter();
  m_colours[m_g][1-m_a] = m_colours[m_q][m_a]   = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g1S0,"1XS_qqbar_g1S0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g1S0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_qqbar_g1S0(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_qqbar_g1S0(args);
    }
  }
  return NULL;
}

XS_qqbar_g1S0::XS_qqbar_g1S0(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R02 = 0.;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.49;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 4.54;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 16./27.*sqr(4.*M_PI)*m_R02/m_mass;
}

double XS_qqbar_g1S0::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(s-m_mass2)-2.*t*u)/(s*sqr(s-m_mass2));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qqbar_g1S0::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[5-m_S][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[5-m_S][1] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g1S0,"1XS_gg_g1S0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g1S0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_gg_g1S0(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_eta_c_1S || kfc==kf_eta_b) return new XS_gg_g1S0(args);
    }
  }
  return NULL;
}

XS_gg_g1S0::XS_gg_g1S0(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R02 = 0.;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.49;
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 4.54;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 1./2.*sqr(4.*M_PI)*m_R02/m_mass;
}

double XS_gg_g1S0::operator()(const Vec4D_Vector& mom) 
{
  double s    = (mom[0]+mom[1]).Abs2(), sM = s-m_mass2, sM2 = sqr(s-m_mass2);
  double t    = (mom[0]-mom[2]).Abs2(), tM = t-m_mass2, tM2 = sqr(t-m_mass2);
  double u    = (mom[0]-mom[3]).Abs2(), uM = u-m_mass2, uM2 = sqr(u-m_mass2);
  double stu  = s*t*u, sum2 = sqr(s)+sqr(t)+sqr(u), mass4 = sqr(m_mass2);
  double mix2 = sqr(s*t)+sqr(s*u)+sqr(t*u);           // dim = GeV^8
  double fac  = 1./(stu*sM2*tM2*uM2);                 // dim = 1/GeV^18
  double sbit = ( pow(s,4)*sM2 * (sM2+2.*mass4) -     // dim = GeV^16 
		  stu * (4./3.  * sum2 * sM*tM*uM +   // dim = GeV^6 * GeV^10
			 16./3. * m_mass2 * mix2 +    // dim = GeV^6 * GeV^10
			 28./3. * mass4 * stu) );     // dim = GeV^6 * GeV^10
  double tbit = ( pow(t,4)*sM2 * (tM2+2.*mass4) -
		  stu * (4./3.  * sum2 * sM*tM*uM +
			 16./3. * m_mass2 * mix2 +
			 28./3. * mass4 * stu) );
  double ubit = ( pow(u,4)*uM2 * (uM2+2.*mass4) -
		  stu * (4./3.  * sum2 * sM*tM*uM +
			 16./3. * m_mass2 * mix2 +
			 28./3. * mass4 * stu) );
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*fac*(sbit+tbit+ubit);
}

bool XS_gg_g1S0::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}


//////////////////////////////////////////////////////////////////////////
//
// Singlet production of 3S1 states (3S1 = 443, 553)
//
//////////////////////////////////////////////////////////////////////////
DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3S1,"1XS_gg_g3S1")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3S1>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_J_psi_1S || kfc==kf_Upsilon_1S) return new XS_gg_g3S1(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_J_psi_1S || kfc==kf_Upsilon_1S) return new XS_gg_g3S1(args);
    }
  }
  return NULL;
}

XS_gg_g3S1::XS_gg_g3S1(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R02 = 0.;
  if (fl[m_S].Kfcode()==kf_J_psi_1S)   m_R02 = 1.11;
  if (fl[m_S].Kfcode()==kf_Upsilon_1S) m_R02 = 9.28;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 5./9.*sqr(4.*M_PI)*m_R02/m_mass;
}

double XS_gg_g3S1::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);
  double M2 = m_mass2/(sM2*tM2*uM2) * (sqr(s)*sM2+sqr(t)*tM2+sqr(u)*uM2);
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_gg_g3S1::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}

