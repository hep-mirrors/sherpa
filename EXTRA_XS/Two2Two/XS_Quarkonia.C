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
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qg_q1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3S1 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;
  public:
    
    XS_qg_q3S1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P0 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qg_q3P0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P1 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qg_q3P1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P2 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qg_q3P2(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g1S0 : public ME2_Base {
  private:
    size_t m_S, m_a;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qqbar_g1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3S1 : public ME2_Base {
  private:
    size_t m_S, m_a;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qqbar_g3S1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

   class XS_qqbar_g3P0 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qqbar_g3P0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P1 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qqbar_g3P1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P2 : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_qqbar_g3P2(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g1S0 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_gg_g1S0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g3S1 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref, ldme;

    public:
    XS_gg_g3S1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P0 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_gg_g3P0(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P1 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_gg_g3P1(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P2 : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R12, m_mass, m_mass2, m_pref, ldme;
  public:
    XS_gg_g3P2(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  // OCTETS

  class XS_qg_q1S0_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q1S0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3S1_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q3S1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P0_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q3P0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P1_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q3P1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qg_q3P2_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qg_q3P2_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g1S0_oct : public ME2_Base {
  private:
    size_t m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g1S0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3S1_oct : public ME2_Base {
  private:
    size_t m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g3S1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

   class XS_qqbar_g3P0_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g3P0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P1_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g3P1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P2_oct : public ME2_Base {
  private:
    size_t m_g, m_q, m_S, m_a;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_qqbar_g3P2_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g1S0_oct : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g1S0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };

  class XS_gg_g3S1_oct : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g3S1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P0_oct : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g3P0_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P1_oct : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g3P1_oct(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };

  class XS_gg_g3P2_oct : public ME2_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
  public:
    XS_gg_g3P2_oct(const External_ME_Args& args);
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
      if ((kfc==kf_eta_c_1S && fl[2].Kfcode() == kf_c) || (kfc==kf_eta_b && fl[2].Kfcode() == kf_b)) return new XS_qg_q1S0(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_eta_c_1S && fl[3].Kfcode() == kf_c) || (kfc==kf_eta_b && fl[3].Kfcode() == kf_b)) return new XS_qg_q1S0(args);
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
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.810;
  if (fl[m_S].Kfcode()==kf_eta_b) m_R02 = 6.477;
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
      if ((kfc==kf_eta_c_1S && abs(fl[1].Kfcode()) == kf_c) || (kfc==kf_eta_b && abs(fl[1].Kfcode()) == kf_b)) return new XS_qqbar_g1S0(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_eta_c_1S && abs(fl[1].Kfcode()) == kf_c) || (kfc==kf_eta_b && abs(fl[1].Kfcode()) == kf_b)) return new XS_qqbar_g1S0(args);
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
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.810;
  if (fl[m_S].Kfcode()==kf_eta_b) m_R02 = 6.477;
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
  if (fl[m_S].Kfcode()==kf_eta_c_1S) m_R02 = 0.810;
  if (fl[m_S].Kfcode()==kf_eta_b) m_R02 = 6.477;
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
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
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
      if ((kfc==kf_J_psi_1S || kfc == kf_psi_2S || kfc==kf_Upsilon_1S || kfc == kf_Upsilon_2S)) return new XS_gg_g3S1(args); 
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_J_psi_1S || kfc == kf_psi_2S || kfc==kf_Upsilon_1S || kfc == kf_Upsilon_2S)) return new XS_gg_g3S1(args);
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
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["R02_443"].SetDefault(0.810);
  pss["R02_100443"].SetDefault(0.529);
  pss["R02_553"].SetDefault(6.477);
  pss["R02_100553"].SetDefault(3.234);
  if (fl[m_S].Kfcode()==kf_J_psi_1S)   m_R02 = pss["R02_443"].Get<double>();
  if (fl[m_S].Kfcode()==kf_psi_2S)     m_R02 = pss["R02_100443"].Get<double>();
  if (fl[m_S].Kfcode()==kf_Upsilon_1S) m_R02 = pss["R02_553"].Get<double>();
  if (fl[m_S].Kfcode()==kf_Upsilon_2S) m_R02 = pss["R02_100553"].Get<double>();
  // std::cout << "ldme = " << ldme << std::endl;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");

}

double XS_gg_g3S1::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = s+t+u, sM2 = sqr(s-M2), tM2 = sqr(t-M2),  uM2 = sqr(u-M2);
  m_pref   = (5./9.)*sqr(4.*M_PI)*sqrt(M2)*m_R02;
  double all = sqr(s)/(uM2*tM2)+sqr(t)/(uM2*sM2)+sqr(u)/(sM2*tM2);
  // double prefix = 10./(9.*64.);
  // double x = M2/s;
  // double first = pow(m_alphaS,3)*CouplingFactor(3,0)*sqr(M_PI)/pow((sqrt(M2)),5);
  // double fgg = prefix*256*sqr(x)/(3*sqr(-1+x)*pow((1+x),3))*(2+x+2*sqr(x)-4*pow(x,4)-pow(x,5)+2*sqr(x)*(5+2*x+sqr(x))*log(x));
  // double all2 = first*fgg*ldme;
  // std::cout << "all = " << all << "fgg = "<< fgg << "first = "<< first<< "x = "<< x <<std::endl;
  // return all;
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*all;
}

bool XS_gg_g3S1::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Singlet production of 3P0 states (3P0 = 10441, 10551)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P0,"1XS_qg_q3P0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c0_1P && fl[2].Kfcode() == kf_c) || (kfc==kf_chi_b0_1P && fl[2].Kfcode() == kf_b) ||
      (kfc==kf_chi_b0_2P && fl[2].Kfcode() == kf_b)) return new XS_qg_q3P0(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c0_1P && fl[3].Kfcode() == kf_c) || (kfc==kf_chi_b0_1P && fl[3].Kfcode() == kf_b) || 
      (kfc==kf_chi_b0_2P && fl[3].Kfcode() == kf_b)) return new XS_qg_q3P0(args);
    }
  }
  return NULL;
}

XS_qg_q3P0::XS_qg_q3P0(const External_ME_Args& args):
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
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["R12_10441"].SetDefault(0.075);
  pss["R12_10551"].SetDefault(1.417);
  pss["R12_110551"].SetDefault(1.653);
  if (fl[m_S].Kfcode()==kf_chi_c0_1P) m_R12 = pss["R12_10441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_chi_b0_1P) m_R12 = pss["R12_10551"].Get<double>();
  if (fl[m_S].Kfcode()==kf_chi_b0_2P) m_R12 = pss["R12_110551"].Get<double>();
  //std::cout << "ldme = " << ldme << std::endl;
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 8./9.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qg_q3P0::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(t-3*m_mass2)*(sqr(s)+sqr(u)))/(-t*pow(t-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qg_q3P0::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[m_g][m_a]   = m_colours[5-m_S][m_a] = Flow::Counter();
  m_colours[m_g][1-m_a] = m_colours[m_q][m_a]   = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P0,"1XS_qqbar_g3P0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c0_1P == kf_c && abs(fl[1].Kfcode())== kf_c ) || (kfc==kf_chi_b0_1P && abs(fl[1].Kfcode()) == kf_b)||
      (kfc==kf_chi_b0_2P&& abs(l[1].Kfcode()) == kf_b)) return new XS_qqbar_g3P0(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c0_1P == kf_c && abs(fl[1].Kfcode())== kf_c ) || (kfc==kf_chi_b0_1P && abs(fl[1].Kfcode()) == kf_b)||
      (kfc==kf_chi_b0_2P&& abs(l[1].Kfcode()) == kf_b)) return new XS_qqbar_g3P0(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P0::XS_qqbar_g3P0(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["R12_10441"].SetDefault(0.075);
  pss["R12_10551"].SetDefault(1.417);
  pss["R12_110551"].SetDefault(1.653);
  if (fl[m_S].Kfcode()==kf_chi_c0_1P) m_R12 = pss["R12_10441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_chi_b0_1P) m_R12 = pss["R12_10551"].Get<double>();
  if (fl[m_S].Kfcode()==kf_chi_b0_2P) m_R12 = pss["R12_110551"].Get<double>();
  //std::cout << "ldme = " << ldme << std::endl;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 64./27.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qqbar_g3P0::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = -(sqr(s-3*m_mass2)*(sqr(t)+sqr(u)))/(-s*pow(s-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qqbar_g3P0::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[5-m_S][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[5-m_S][1] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P0,"1XS_gg_g3P0")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P0>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_chi_c0_1P || kfc==kf_chi_b0_1P || kfc==kf_chi_b0_2P) return new XS_gg_g3P0(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_chi_c0_1P || kfc==kf_chi_b0_1P || kfc==kf_chi_b0_2P) return new XS_gg_g3P0(args);
    }
  }
  return NULL;
}

XS_gg_g3P0::XS_gg_g3P0(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["R12_10441"].SetDefault(0.075);
  pss["R12_10551"].SetDefault(1.417);
  pss["R12_110551"].SetDefault(1.653);
  if (fl[m_S].Kfcode()==kf_chi_c0_1P) m_R12 = pss["R12_10441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_chi_b0_1P) m_R12 = pss["R12_10551"].Get<double>();
  if (fl[m_S].Kfcode()==kf_chi_b0_2P) m_R12 = pss["R12_110551"].Get<double>();
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 4.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_gg_g3P0::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double Q = s*t*u;
  double P = s*t+t*u+u*s;
  double M2 = s+t+u; double M4 = sqr(M2); double M8 = sqr(M4);
  double mul = 1/(Q*pow(Q-M2*P,4));
  double all1 = 9.*M4*pow(P,4)*(M8-2.*M4*P+sqr(P))-6.*M2*pow(P,3)*Q*(2.*M8-5.*M4*P+sqr(P));
  double all2 = -sqr(P)*sqr(Q)*(M8+2.*M4*P-sqr(P))+2.*M2*P*pow(Q,3)*(M4-P)+6.*M4*pow(Q,4);
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*mul*(all1+all2);
}

bool XS_gg_g3P0::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter(); //bug fixed. Original: m_colours[0][bit]=...
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Singlet production of 3P1 states (3P1 = 20443, 20553)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P1,"1XS_qg_q3P1")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P1>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c1_1P && fl[2].Kfcode()==kf_c) || (kfc==kf_chi_b1_1P && fl[2].Kfcode()==kf_b) ||
       (kfc==kf_chi_b1_2P && fl[2].Kfcode()==kf_b)) return new XS_qg_q3P1(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c1_1P && fl[3].Kfcode()==kf_c) || (kfc==kf_chi_b1_1P && fl[3].Kfcode()==kf_b) ||
       (kfc==kf_chi_b1_2P && fl[3].Kfcode()==kf_b)) return new XS_qg_q3P1(args);
    }
  }
  return NULL;
}

XS_qg_q3P1::XS_qg_q3P1(const External_ME_Args& args):
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
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;

  if (fl[m_S].Kfcode()==kf_chi_c1_1P)   m_R12 = 3.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b1_1P)   m_R12 = 3.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b1_2P)   m_R12 = 3.*1.653;
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 16./3.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qg_q3P1::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = ((-t)*(sqr(s)+sqr(u))-(4*m_mass2*s*u))/(pow(t-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qg_q3P1::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[m_g][m_a]   = m_colours[5-m_S][m_a] = Flow::Counter();
  m_colours[m_g][1-m_a] = m_colours[m_q][m_a]   = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P1,"1XS_qqbar_g3P1")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
            PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P1>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
  fl[2].IsGluon() && fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c1_1P && abs(fl[1].Kfcode())==kf_c)|| (kfc==kf_chi_b1_1P && abs(fl[1].Kfcode())==kf_b) || 
      (kfc==kf_chi_b1_2P && abs(fl[1].Kfcode())==kf_b)) return new XS_qqbar_g3P1(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
  fl[3].IsGluon() && fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c1_1P && abs(fl[1].Kfcode())==kf_c)|| (kfc==kf_chi_b1_1P && abs(fl[1].Kfcode())==kf_b) || 
      (kfc==kf_chi_b1_2P && abs(fl[1].Kfcode())==kf_b)) return new XS_qqbar_g3P1(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P1::XS_qqbar_g3P1(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  if (fl[m_S].Kfcode()==kf_chi_c1_1P)   m_R12 = 3.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b1_1P)   m_R12 = 3.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b1_2P)   m_R12 = 3.*1.653;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = -8./3.*16./3.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qqbar_g3P1::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = ((-s)*(sqr(t)+sqr(u))-(4*m_mass2*t*u))/(pow(s-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qqbar_g3P1::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[5-m_S][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[5-m_S][1] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P1,"1XS_gg_g3P1")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P1>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_chi_c1_1P || kfc==kf_chi_b1_1P || kfc==kf_chi_b1_2P) return new XS_gg_g3P1(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_chi_c1_1P || kfc==kf_chi_b1_1P || kfc==kf_chi_b1_2P) return new XS_gg_g3P1(args);
    }
  }
  return NULL;
}

XS_gg_g3P1::XS_gg_g3P1(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  if (fl[m_S].Kfcode()==kf_chi_c1_1P)   m_R12 = 3.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b1_1P)   m_R12 = 3.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b1_2P)   m_R12 = 3.*1.653;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 12.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_gg_g3P1::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double Q = s*t*u;
  double P = s*t+t*u+u*s; double P2 = sqr(P);
  double M2 = s+t+u; double M4 = sqr(M2);
  double nom = P2*(M2*P2*(M4-4.*P)+2.*Q*(-sqr(M4)+5.*M4*P+P2)-15.*M2*sqr(Q));
  double dnom = pow(Q-M2*P,4);
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*nom/dnom;
}

bool XS_gg_g3P1::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Singlet production of 3P2 states (3P2 = 445, 555)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P2,"1XS_qg_q3P2")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P2>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c2_1P && fl[2].Kfcode()==kf_c) || (kfc==kf_chi_b2_1P && fl[2].Kfcode())==kf_b || 
      (kfc==kf_chi_b2_2P && fl[2].Kfcode())==kf_b) return new XS_qg_q3P2(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c2_1P && fl[3].Kfcode()==kf_c) || (kfc==kf_chi_b2_1P && fl[3].Kfcode())==kf_b || 
      (kfc==kf_chi_b2_2P && fl[3].Kfcode())==kf_b)return new XS_qg_q3P2(args);
    }
  }
  return NULL;
}

XS_qg_q3P2::XS_qg_q3P2(const External_ME_Args& args):
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
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  if (fl[m_S].Kfcode()==kf_chi_c2_1P)   m_R12 = 5.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b2_1P)   m_R12 = 5.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b2_2P)   m_R12 = 5.*1.653;
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 16./9.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qg_q3P2::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(t-m_mass2)*(sqr(t)+6*sqr(m_mass2))-2*s*u*(sqr(t)-6*m_mass2*(t-m_mass2)))/(-t*pow(t-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qg_q3P2::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[m_g][m_a]   = m_colours[5-m_S][m_a] = Flow::Counter();
  m_colours[m_g][1-m_a] = m_colours[m_q][m_a]   = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P2,"1XS_qqbar_g3P2")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P2>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if ((kfc==kf_chi_c2_1P && abs(fl[1].Kfcode())==kf_c) || (kfc==kf_chi_b2_1P && abs(fl[2].Kfcode())==kf_b) || 
      (kfc==kf_chi_b2_2P && abs(fl[1].Kfcode())==kf_b)) return new XS_qqbar_g3P2(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if ((kfc==kf_chi_c2_1P && abs(fl[1].Kfcode())==kf_c) || (kfc==kf_chi_b2_1P && abs(fl[2].Kfcode())==kf_b) || 
      (kfc==kf_chi_b2_2P && abs(fl[1].Kfcode())==kf_b)) return new XS_qqbar_g3P2(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P2::XS_qqbar_g3P2(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  if (fl[m_S].Kfcode()==kf_chi_c2_1P)   m_R12 = 5.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b2_1P)   m_R12 = 5.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b2_2P)   m_R12 = 5.*1.653;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = -8./3.*16./9.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_qqbar_g3P2::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(s-m_mass2)*(sqr(s)+6*sqr(m_mass2))-2*t*u*(sqr(s)-6*m_mass2*(s-m_mass2)))/(-s*pow(s-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2;
}

bool XS_qqbar_g3P2::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[5-m_S][0] = Flow::Counter(); 
  m_colours[m_a][1]   = m_colours[5-m_S][1] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P2,"1XS_gg_g3P2")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P2>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_chi_c2_1P || kfc==kf_chi_b2_1P || kfc==kf_chi_b2_2P) return new XS_gg_g3P2(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_chi_c2_1P || kfc==kf_chi_b2_1P || kfc==kf_chi_b2_2P) return new XS_gg_g3P2(args);
    }
  }
  return NULL;
}

XS_gg_g3P2::XS_gg_g3P2(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); m_R12 = 0.;
  if (fl[m_S].Kfcode()==kf_chi_c2_1P)   m_R12 = 5.*0.075;
  if (fl[m_S].Kfcode()==kf_chi_b2_1P)   m_R12 = 5.*1.417;
  if (fl[m_S].Kfcode()==kf_chi_b2_2P)   m_R12 = 5.*1.653;
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 4.*sqr(4.*M_PI)*m_R12/pow(m_mass,3);
}

double XS_gg_g3P2::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double Q = s*t*u;
  double P = s*t+t*u+u*s; double P2 = sqr(P);
  double M2 = s+t+u; double M4 = sqr(M2); double M8 = sqr(M4);
  double nom = 12.*M4*sqr(P2)*(M8-2*M4*P+P2)-3.*M2*pow(P,3)*Q*(8.*M8-M4*P+4.*P2)+2.*P2*sqr(Q)*(-7.*M8+43.*M4*P+P2);
  double nom2 = M2*P*pow(Q,3)*(16.*M4-61*P)+12.*M4*pow(Q,4);
  double dnom = Q*pow(Q-M2*P,4);
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*(nom+nom2)/dnom;
}

bool XS_gg_g3P2::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  m_colours[0][bit] = m_colours[1][1-bit]     = Flow::Counter();
  m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
  m_colours[1][bit] = m_colours[5-m_S][bit]   = Flow::Counter();
  return true;
}

//For S
//LDME = 9/(2*M_PI)*M_R02

//For P
//LDME = 3/(4*M_PI)*(2J+1)*m_R12 [1503.08439]

//values from [9503356]


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
//                                OCTETS
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//LDMEs from: https://arxiv.org/pdf/1906.10049, https://arxiv.org/html/2501.15575v1
//LDME rules: https://arxiv.org/pdf/1411.5287


//////////////////////////////////////////////////////////////////////////
//
// Octet production of eta states (1S0 = 9900441, 9900551)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q1S0_oct,"1XS_qg_q1S0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q1S0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S) return new XS_qg_q1S0_oct(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S) return new XS_qg_q1S0_oct(args);
    }
  }
  return NULL;
}

XS_qg_q1S0_oct::XS_qg_q1S0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_990000441"].SetDefault(0.018/3.*1000);
  pss["LDME_990000443"].SetDefault(0.012/3.*1000);
  pss["LDME_990100443"].SetDefault(0.005/3.*1000);
  pss["LDME_990000551"].SetDefault(0.0231/3.*1000);
  pss["LDME_990000553"].SetDefault(0.0151/3.*1000);
  pss["LDME_990100553"].SetDefault(0.0124/3.*1000);
  pss["LDME_990200553"].SetDefault(0.0101/3.*1000);
  if (fl[m_S].Kfcode()==kf_1S0_c_8_eta_c)      LDME = 0.001*pss["LDME_990000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_990000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_psi_2S)     LDME = 0.001*pss["LDME_990100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_eta_b)      LDME = 0.001*pss["LDME_990000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_990000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_990100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_990200553"].Get<double>(); 
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qg_q1S0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2= (sqr(s)+sqr(u))/(t*sqr(t-m_mass2));
  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  return -5./(72.*m_mass)*m_pref*M2*LDME;
}

bool XS_qg_q1S0_oct::SetColours(const Vec4D_Vector& mom) 
{
  if ((m_g+2) == m_S){
    size_t cross = ran->Get()<0.5 ? 0 : 1;
    if (cross==0){
      m_colours[m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_g][1-m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_S][1-m_a] = Flow::Counter();
  
    }
    if (cross==1){
      m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_g][m_a] = Flow::Counter();
    }
  }
  if ((m_g+2) != m_S){
    m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
    m_colours[m_g][m_a] = m_colours[5-m_S][m_a] = Flow::Counter();
    m_colours[1-m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
  }
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g1S0_oct,"1XS_qqbar_g1S0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g1S0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsOctetMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S) return new XS_qqbar_g1S0_oct(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S) return new XS_qqbar_g1S0_oct(args);
    }
  }
  return NULL;
}

XS_qqbar_g1S0_oct::XS_qqbar_g1S0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_990000441"].SetDefault(0.018/3.*1000);
  pss["LDME_990000443"].SetDefault(0.012/3.*1000);
  pss["LDME_990100443"].SetDefault(0.005/3.*1000);
  pss["LDME_990000551"].SetDefault(0.0231/3.*1000);
  pss["LDME_990000553"].SetDefault(0.0151/3.*1000);
  pss["LDME_990100553"].SetDefault(0.0124/3.*1000);
  pss["LDME_990200553"].SetDefault(0.0101/3.*1000);
  if (fl[m_S].Kfcode()==kf_1S0_c_8_eta_c)      LDME = 0.001*pss["LDME_990000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_990000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_psi_2S)     LDME = 0.001*pss["LDME_990100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_eta_b)      LDME = 0.001*pss["LDME_990000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_990000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_990100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_990200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qqbar_g1S0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(t)+sqr(u))/(s*sqr(s-m_mass2));
  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  return 5./(27.*m_mass)*m_pref*M2*LDME;
}

bool XS_qqbar_g1S0_oct::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[3-m_a][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[2+m_a][1] = Flow::Counter();
  m_colours[2][m_a]   = m_colours[3][1-m_a] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g1S0_oct,"1XS_gg_g1S0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g1S0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&

	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S) return new XS_gg_g1S0_oct(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_1S0_c_8_eta_c || kfc==kf_1S0_c_8_J_psi_1S   || kfc==kf_1S0_c_8_psi_2S     ||
          kfc==kf_1S0_b_8_eta_b || kfc==kf_1S0_b_8_Upsilon_1S || kfc==kf_1S0_b_8_Upsilon_2S || kfc==kf_1S0_b_8_Upsilon_3S)return new XS_gg_g1S0_oct(args);
    }
  }
  return NULL;
}

XS_gg_g1S0_oct::XS_gg_g1S0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_990000441"].SetDefault(0.018/3.*1000);
  pss["LDME_990000443"].SetDefault(0.012/3.*1000);
  pss["LDME_990100443"].SetDefault(0.005/3.*1000);
  pss["LDME_990000551"].SetDefault(0.0231/3.*1000);
  pss["LDME_990000553"].SetDefault(0.0151/3.*1000);
  pss["LDME_990100553"].SetDefault(0.0124/3.*1000);
  pss["LDME_990200553"].SetDefault(0.0101/3.*1000);
  if (fl[m_S].Kfcode()==kf_1S0_c_8_eta_c)      LDME = 0.001*pss["LDME_990000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_990000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_c_8_psi_2S)     LDME = 0.001*pss["LDME_990100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_eta_b)      LDME = 0.001*pss["LDME_990000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_990000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_990100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_1S0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_990200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_gg_g1S0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s    = (mom[0]+mom[1]).Abs2(), sM = s-m_mass2, sM2 = sqr(s-m_mass2);
  double t    = (mom[0]-mom[2]).Abs2(), tM = t-m_mass2, tM2 = sqr(t-m_mass2);
  double u    = (mom[0]-mom[3]).Abs2(), uM = u-m_mass2, uM2 = sqr(u-m_mass2);
  double M2   = (sqr(s)*sM2+s*t*u*(m_mass2-2.*s)+sqr(t*u))*(sqr(sqr(s)-m_mass2*s+sqr(m_mass2))-t*u*(2*sqr(t)+3*t*u+2*sqr(u)))/(s*t*u*sM2*tM2*uM2);
  m_pref      = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  return 5./(16.*m_mass)*m_pref*M2*LDME;
}

bool XS_gg_g1S0_oct::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Octet production of 3S1 states (3S1 = 99000443, 99000553, ...)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3S1_oct,"1XS_qg_q3S1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3S1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P ) return new XS_qg_q3S1_oct(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P) return new XS_qg_q3S1_oct(args);
    }
  }
  return NULL;
}

XS_qg_q3S1_oct::XS_qg_q3S1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_991000441"].SetDefault(0.0180*1000);
  pss["LDME_991000443"].SetDefault(0.0013*1000);
  pss["LDME_991100443"].SetDefault(0.0033*1000);
  pss["LDME_991010441"].SetDefault(0.00187*1000);
  // pss["LDME_991020443"].SetDefault(3*0.00187);
  // pss["LDME_991000445"].SetDefault(5*0.00187);

  pss["LDME_991000551"].SetDefault(0.0121*1000);
  pss["LDME_991000553"].SetDefault(0.0477*1000);
  pss["LDME_991100553"].SetDefault(0.121*1000);
  pss["LDME_991200553"].SetDefault(0.307*1000);
  pss["LDME_991010551"].SetDefault(0.1008*1000);
  // pss["LDME_991020553"].SetDefault(3*0.1008);
  // pss["LDME_991000555"].SetDefault(5*0.1008);

  if (fl[m_S].Kfcode()==kf_3S1_c_8_eta_c)      LDME = 0.001*pss["LDME_991000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_991000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_psi_2S)     LDME = 0.001*pss["LDME_991100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c0_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c1_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c2_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 

  if (fl[m_S].Kfcode()==kf_3S1_b_8_eta_b)      LDME = 0.001*pss["LDME_991000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_991000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_991100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_991200553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b0_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b1_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b2_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 

  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qg_q3S1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), s2 = sqr(s);
  double t  = (mom[0]-mom[2]).Abs2(); 
  double u  = (mom[0]-mom[3]).Abs2(), u2 = sqr(u);
  
  double M2 = s+t+u, sM = sqr(s-M2), tM = sqr(t-M2), uM = sqr(u-M2);
  double heq0 = 2.*M2*t;
  double heq1 = ((s2+u2+2.*M2*t)*sM-2.*M2*s*t*u)/(s*u); 
  double dnom = sqr(sM*tM);
  double nom = (4.*(s2+u2)-s*u);
  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  return -1./(108.*pow(m_mass,3))*m_pref*(heq0+heq1)/dnom*nom*LDME;
}

bool XS_qg_q3S1_oct::SetColours(const Vec4D_Vector& mom) 
{
  if ((m_g+2) == m_S){
    size_t cross = ran->Get()<0.5 ? 0 : 1;
    if (cross==0){
      m_colours[m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_g][1-m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_S][1-m_a] = Flow::Counter();
  
    }
    if (cross==1){
      m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_g][m_a] = Flow::Counter();
    }
  }
  if ((m_g+2) != m_S){
    m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
    m_colours[m_g][m_a] = m_colours[5-m_S][m_a] = Flow::Counter();
    m_colours[1-m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
  }
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3S1_oct,"1XS_qqbar_g3S1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3S1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsOctetMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P) return new XS_qqbar_g3S1_oct(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P) return new XS_qqbar_g3S1_oct(args);
    }
  }
  return NULL;
}

XS_qqbar_g3S1_oct::XS_qqbar_g3S1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_991000441"].SetDefault(0.0180*1000);
  pss["LDME_991000443"].SetDefault(0.0013*1000);
  pss["LDME_991100443"].SetDefault(0.0033*1000);
  pss["LDME_991010441"].SetDefault(0.00187*1000);
  // pss["LDME_991020443"].SetDefault(3*0.00187);
  // pss["LDME_991000445"].SetDefault(5*0.00187);

  pss["LDME_991000551"].SetDefault(0.0121*1000);
  pss["LDME_991000553"].SetDefault(0.0477*1000);
  pss["LDME_991100553"].SetDefault(0.121*1000);
  pss["LDME_991200553"].SetDefault(0.307*1000);
  pss["LDME_991010551"].SetDefault(0.1008*1000);
  // pss["LDME_991020553"].SetDefault(3*0.1008);
  // pss["LDME_991000555"].SetDefault(5*0.1008);

  if (fl[m_S].Kfcode()==kf_3S1_c_8_eta_c)      LDME = 0.001*pss["LDME_991000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_991000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_psi_2S)     LDME = 0.001*pss["LDME_991100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c0_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c1_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c2_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 

  if (fl[m_S].Kfcode()==kf_3S1_b_8_eta_b)      LDME = 0.001*pss["LDME_991000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_991000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_991100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_991200553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b0_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b1_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b2_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qqbar_g3S1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(); 
  double t  = (mom[0]-mom[2]).Abs2(); double t2 = sqr(t);
  double u  = (mom[0]-mom[3]).Abs2(); double u2 = sqr(u);

  double alphaS = 0.41;
  double M2 = s+t+u;
  double heq0 = 4.*M2*s;
  double heq1 = (sqr(s)+sqr(M2))*(t2+u2)/(t*u);
  double nom = 4.*(t2+u2)-t*u;
  double dnom = pow(s-M2,4);

  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  return 2./(81.*pow(m_mass,3))*m_pref*(heq0+heq1)/dnom*nom*LDME;
}

bool XS_qqbar_g3S1_oct::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[3-m_a][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[2+m_a][1] = Flow::Counter();
  m_colours[2][m_a]   = m_colours[3][1-m_a] = Flow::Counter();
  return true;
}


DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3S1_oct,"1XS_gg_g3S1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3S1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P) return new XS_gg_g3S1_oct(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3S1_c_8_J_psi_1S  || kfc==kf_3S1_c_8_psi_2S     || kfc==kf_3S1_c_8_eta_c      || 
          kfc==kf_3S1_c_8_chi_c0_1P || kfc==kf_3S1_c_8_chi_c1_1P  || kfc==kf_3S1_c_8_chi_c2_1P  ||
          kfc==kf_3S1_b_8_Upsilon_1S|| kfc==kf_3S1_b_8_Upsilon_2S || kfc==kf_3S1_b_8_Upsilon_3S || kfc==kf_3S1_b_8_eta_b      ||
          kfc==kf_3S1_b_8_chi_b0_1P || kfc==kf_3S1_b_8_chi_b1_1P  || kfc==kf_3S1_b_8_chi_b2_1P) return new XS_gg_g3S1_oct(args);
    }
  }
  return NULL;
}

XS_gg_g3S1_oct::XS_gg_g3S1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_991000441"].SetDefault(0.0180*1000);
  pss["LDME_991000443"].SetDefault(0.0013*1000);
  pss["LDME_991100443"].SetDefault(0.0033*1000);
  pss["LDME_991010441"].SetDefault(0.00187*1000);
  // pss["LDME_991020443"].SetDefault(3*0.00187);
  // pss["LDME_991000445"].SetDefault(5*0.00187);

  pss["LDME_991000551"].SetDefault(0.0121*1000);
  pss["LDME_991000553"].SetDefault(0.0477*1000);
  pss["LDME_991100553"].SetDefault(0.121*1000);
  pss["LDME_991200553"].SetDefault(0.307*1000);
  pss["LDME_991010551"].SetDefault(0.1008*1000);
  // pss["LDME_991020553"].SetDefault(3*0.1008);
  // pss["LDME_991000555"].SetDefault(5*0.1008);

  if (fl[m_S].Kfcode()==kf_3S1_c_8_eta_c)      LDME = 0.001*pss["LDME_991000441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_991000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_psi_2S)     LDME = 0.001*pss["LDME_991100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c0_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c1_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_c_8_chi_c2_1P)  LDME = 0.001*pss["LDME_991010441"].Get<double>(); 

  if (fl[m_S].Kfcode()==kf_3S1_b_8_eta_b)      LDME = 0.001*pss["LDME_991000551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_991000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_991100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_991200553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b0_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b1_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3S1_b_8_chi_b2_1P)  LDME = 0.001*pss["LDME_991010551"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_gg_g3S1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double M2 = s+t+u;
  double heq0 = 2.*s*M2*(sqr(t)+sqr(u))*t*u;
  double heq1 = sqr(s)*(sqr(sM2)+pow(t,4)+pow(u,4)+2.*sqr(M2)*sqr(t*u/s));
  double nom = 27.*(s*t+t*u+u*s)-19.*sqr(M2);
  double dnom = sM2*sM2*tM2*uM2;
  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);
  // std::cout<<"LDME = "<<LDME<< " all = " << -1./(144.*pow(m_mass,3))*m_pref*(heq0+heq1)*nom/dnom*LDME<< std::endl;
  return -1./(144.*pow(m_mass,3))*m_pref*(heq0+heq1)*nom/dnom*LDME;
}

bool XS_gg_g3S1_oct::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Octet production of 3P0 states (3P0 = 9910441, 9910551)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P0_oct,"1XS_qg_q3P0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S|| kfc==kf_3P0_b_8_Upsilon_3S) return new XS_qg_q3P0_oct(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S|| kfc==kf_3P0_b_8_Upsilon_3S)  return new XS_qg_q3P0_oct(args);
    }
  }
  return NULL;
}

XS_qg_q3P0_oct::XS_qg_q3P0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(1 * 3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  
  if (fl[m_S].Kfcode()==kf_3P0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_c_8_psi_2S)     LDME = 0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_992200553"].Get<double>(); 
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = -5./(54.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qg_q3P0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double M2 = (sqr(t-3.*m_mass2)*(sqr(s)+sqr(u)))/(t*pow(t-m_mass2,4));
  return pow(m_alphaS,3)*CouplingFactor(3,0)*m_pref*M2*LDME;
}

bool XS_qg_q3P0_oct::SetColours(const Vec4D_Vector& mom) 
{
  if ((m_g+2) == m_S){
    size_t cross = ran->Get()<0.5 ? 0 : 1;
    if (cross==0){
      m_colours[m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_g][1-m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_S][1-m_a] = Flow::Counter();
  
    }
    if (cross==1){
      m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_g][m_a] = Flow::Counter();
    }
  }
  if ((m_g+2) != m_S){
    m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
    m_colours[m_g][m_a] = m_colours[5-m_S][m_a] = Flow::Counter();
    m_colours[1-m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
  }
  return true;
}


DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P0_oct,"1XS_qqbar_g3P0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsOctetMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S|| kfc==kf_3P0_b_8_Upsilon_3S)  return new XS_qqbar_g3P0_oct(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if  (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S || kfc==kf_3P0_b_8_Upsilon_3S)  return new XS_qqbar_g3P0_oct(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P0_oct::XS_qqbar_g3P0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(1 * 3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  
  if (fl[m_S].Kfcode()==kf_3P0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_c_8_psi_2S)     LDME = 0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_992100553"].Get<double>();
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_992200553"].Get<double>();  
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 20./(81.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qqbar_g3P0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();

  double M2 = sqr(s-3.*m_mass2)*(sqr(t)+sqr(u))/(s*pow(s-m_mass2,4));
  return CouplingFactor(3,0)*m_pref*M2*LDME;
}

bool XS_qqbar_g3P0_oct::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[3-m_a][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[2+m_a][1] = Flow::Counter();
  m_colours[2][m_a]   = m_colours[3][1-m_a] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P0_oct,"1XS_gg_g3P0_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P0_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S|| kfc==kf_3P0_b_8_Upsilon_3S)  return new XS_gg_g3P0_oct(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P0_c_8_J_psi_1S   || kfc==kf_3P0_c_8_psi_2S    ||
          kfc==kf_3P0_b_8_Upsilon_1S || kfc==kf_3P0_b_8_Upsilon_2S|| kfc==kf_3P0_b_8_Upsilon_3S)  return new XS_gg_g3P0_oct(args);
    }
  }
  return NULL;
}

XS_gg_g3P0_oct::XS_gg_g3P0_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(1 * 3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  
  if (fl[m_S].Kfcode()==kf_3P0_c_8_J_psi_1S)   LDME = 0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_c_8_psi_2S)     LDME = 0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_1S) LDME = 0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_2S) LDME = 0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P0_b_8_Upsilon_3S) LDME = 0.001*pss["LDME_992200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_gg_g3P0_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2); double s2 = sqr(s), s4 = sqr(s2), s6 = pow(s,6), s8 = sqr(s4);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double z = sqrt(t*u); double z2 = sqr(z); double z4 = sqr(z2); double z6 = pow(z,6); double z8 = sqr(z4);
  double M = m_mass;
  double line1 = s2*z4*pow(s2-z2,4)+sqr(M)*s*z2*sqr(s2+z2)*(3.*s2-2.*z2)*(2.*s4-6.*s2*z2+3*z4);
  double line2 = pow(M,4)*(9.*pow(s,12)-84.*pow(s,10)*z2+265.*s8*z4-382.*s6*z6+276.*s4*z8-88.*s2*pow(z,10)+9.*pow(z,12));
  double line3 = -pow(M,6)*s*(54.*pow(s,10)-357.*s8*z2+844.*s6*z4-898.*s4*z6+439*s2*z8-81*pow(z,10));
  double line4 = pow(M,8)*(153.*pow(s,10)-798.*s8*z2+1415.*s6*z4-1041.*s4*z6+301*s2*z8-18*pow(z,10));
  double line5 = -pow(M,10)*s*(270.*s8-1089.*s6*z2+1365.*s4*z4-616.*s2*z6+87.*z8);
  double line6 = pow(M,12)*(324.*s8-951.*s6*z2+769.*s4*z4-189.*s2*z6+9.*z8);
  double line7 = -9.*pow(M,14)*s*(6.*s2-z2)*(5.*s4-9.*s2*z2+3.*z4);
  double line8 = 3.*pow(M,16)*s2*(51*s4-59*s2*z2+12*z4)-27.*pow(M,18)*pow(s,3)*(2.*s2-z2)+9.*pow(M,20)*s4;
  double dnom = s*z2*pow(s-sqr(M),4)*pow(s*sqr(M)+z2,4);
  double all = line1+line2+line3+line4+line5+line6+line7+line8;
  m_pref = 5./12.*pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0)/pow(m_mass,3);
  // std::cout<< "------gg--3P0-----"<< m_pref*all/dnom*LDME<< std::endl;
  return m_pref*all/dnom*LDME;
}

bool XS_gg_g3P0_oct::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Octet production of 3P1 states (3P1 = 9920443, 9920553)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P1_oct,"1XS_qg_q3P1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_qg_q3P1_oct(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_qg_q3P1_oct(args);
    }
  }
  return NULL;
}

XS_qg_q3P1_oct::XS_qg_q3P1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P1_c_8_J_psi_1S)   LDME = 3.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_c_8_psi_2S)     LDME = 3.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_1S) LDME = 3.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_2S) LDME = 3.*0.001*pss["LDME_992100553"].Get<double>();
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_3S) LDME = 3.*0.001*pss["LDME_992200553"].Get<double>();  
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qg_q3P1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(); double sM = s-m_mass2;
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double heq0 = (t*(sqr(s*(sM))+sqr(u*(s+m_mass2))));
  double heq1 = 4.*m_mass2*(s*u*(sqr(t)+t*u+sqr(u)));
  double dnom = (pow(t-m_mass2,4)*sqr(sM));
  m_pref   = -5./(27.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);

  return CouplingFactor(3,0)*m_pref*(heq0+heq1)/dnom*LDME;
}

bool XS_qg_q3P1_oct::SetColours(const Vec4D_Vector& mom) 
{
 if ((m_g+2) == m_S){
    size_t cross = ran->Get()<0.5 ? 0 : 1;
    if (cross==0){
      m_colours[m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_g][1-m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_S][1-m_a] = Flow::Counter();
  
    }
    if (cross==1){
      m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_g][m_a] = Flow::Counter();
    }
  }
  if ((m_g+2) != m_S){
    m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
    m_colours[m_g][m_a] = m_colours[5-m_S][m_a] = Flow::Counter();
    m_colours[1-m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
  }
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P1_oct,"1XS_qqbar_g3P1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsOctetMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_qqbar_g3P1_oct(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_qqbar_g3P1_oct(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P1_oct::XS_qqbar_g3P1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass);
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P1_c_8_J_psi_1S)   LDME = 3.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_c_8_psi_2S)     LDME = 3.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_1S) LDME = 3.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_2S) LDME = 3.*0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_3S) LDME = 3.*0.001*pss["LDME_992200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 40./(81.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qqbar_g3P1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double heq0 = s*(sqr(t)+sqr(u));
  double heq1 = 4.*m_mass2*t*u;
  double dnom = pow(s-m_mass2,4);
  return CouplingFactor(3,0)*m_pref*(heq0+heq1)/dnom*LDME;
}

bool XS_qqbar_g3P1_oct::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[3-m_a][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[2+m_a][1] = Flow::Counter();
  m_colours[2][m_a]   = m_colours[3][1-m_a] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P1_oct,"1XS_gg_g3P1_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P1_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_gg_g3P1_oct(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P1_c_8_J_psi_1S   || kfc==kf_3P1_c_8_psi_2S    ||
          kfc==kf_3P1_b_8_Upsilon_1S || kfc==kf_3P1_b_8_Upsilon_2S|| kfc==kf_3P1_b_8_Upsilon_3S)  return new XS_gg_g3P1_oct(args);
    }
  }
  return NULL;
}

XS_gg_g3P1_oct::XS_gg_g3P1_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P1_c_8_J_psi_1S)   LDME = 3.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_c_8_psi_2S)     LDME = 3.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_1S) LDME = 3.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_2S) LDME = 3.*0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P1_b_8_Upsilon_3S) LDME = 3.*0.001*pss["LDME_992200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_gg_g3P1_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2); double s2 = sqr(s), s4 = sqr(s2), s6 = pow(s,6), s8 = sqr(s4);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double M = m_mass, M2 = m_mass2;

  double z = sqrt(t*u), z2 = sqr(z), z4 = sqr(z2), z6 = pow(z,6), z8 = sqr(z4);
  double heq0 = s*z2*(sqr(s2-z2)-2.*M2*s*z2-sqr(M2)*(s2+2.*z2)+pow(M,8))*(sqr(s2-z2)-M2*s*(2.*s2-z2)+sqr(M2)*s2);
  double heq1l1 = M2*(2.*sqr(s2-z2)*(s6-4.*s4*z2+s2*z4-z6)-M2*s*(2.*s2-z2)*(5.*s6-17.*s4*z2+9.*s2*z4-z6)+sqr(M2)*(21.*s8-49.*s6*z2+21.*s4*z4-4.*s2*z6+z8));
  double heq1l2 = M2*(-pow(M,6)*s*(24.*s6-30.*s4*z2+6.*s2*z4-z6)+pow(M,8)*s2*(16.*s4-9.*s2*z2+2*z4)-pow(M,10)*pow(s,3)*(6.*s2-z2)+pow(M,12)*s4);
  double heq = heq0 + heq1l1 + heq1l2;
  double dnom = pow(s-M2,4)*pow(s*M2+z2,4);

  m_pref = 5./(6.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0);

  // std::cout<< "------gg-3P1----"<< m_pref*heq/dnom*LDME<<std::endl;
  return m_pref*heq/dnom*LDME;
}

bool XS_gg_g3P1_oct::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Octet production of 3P2 states (3P2 = 9900445, 9900555)
//
//////////////////////////////////////////////////////////////////////////

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qg_q3P2_oct,"1XS_qg_q3P2_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qg_q3P2_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[2]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[2]==fl[1])) &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_qg_q3P2_oct(args);
    }
    if ( ((fl[0].IsQuark() && fl[1].IsGluon() && fl[3]==fl[0]) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[3]==fl[1])) &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_qg_q3P2_oct(args);
    }
  }
  return NULL;
}

XS_qg_q3P2_oct::XS_qg_q3P2_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsGluon()) m_g = i;
    if (i<2 && fl[i].IsQuark()) m_q = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P2_c_8_J_psi_1S)   LDME = 5.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_c_8_psi_2S)     LDME = 5.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_1S) LDME = 5.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_2S) LDME = 5.*0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_3S) LDME = 5.*0.001*pss["LDME_992200553"].Get<double>(); 
  m_a      = fl[m_q].IsAnti() ? 1 : 0; 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = -1./(27.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }

}

double XS_qg_q3P2_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(); double s2 = sqr(s); double sM4 = pow(s-m_mass2,4);
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double m_mass4 = sqr(m_mass2);
  double heq0 = t*(s2+sqr(u)+12.*m_mass2*s*sqr(u)*(s2+m_mass2*s+m_mass4)/sM4);
  double heq1 = 12.*m_mass2*s*u*(sqr(s-m_mass2)*(s2+m_mass4)-sqr(s+m_mass2)*t*u)/sM4;
  double heq2 = (6.*m_mass4*(s2+sqr(u)+2.*s2*t*u*((s-m_mass2)*(2.*t+u)-sqr(u))/sM4))/t;
  double dnom = pow(t-m_mass2,4);

  return CouplingFactor(3,0)*m_pref*(heq0+heq1+heq2)/dnom*LDME;
}

bool XS_qg_q3P2_oct::SetColours(const Vec4D_Vector& mom) 
{
  if ((m_g+2) == m_S){
    size_t cross = ran->Get()<0.5 ? 0 : 1;
    if (cross==0){
      m_colours[m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_g][1-m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_S][1-m_a] = Flow::Counter();
  
    }
    if (cross==1){
      m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
      m_colours[1-m_g][m_a]   = m_colours[m_S][m_a] = Flow::Counter();
      m_colours[5-m_S][m_a]   = m_colours[m_g][m_a] = Flow::Counter();
    }
  }
  if ((m_g+2) != m_S){
    m_colours[m_g][1-m_a] = m_colours[m_S][1-m_a] = Flow::Counter();
    m_colours[m_g][m_a] = m_colours[5-m_S][m_a] = Flow::Counter();
    m_colours[1-m_g][m_a] = m_colours[m_S][m_a] = Flow::Counter();
  }
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqbar_g3P2_oct,"1XS_qqbar_g3P2_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_qqbar_g3P2_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[2].IsGluon() && fl[3].IsOctetMeson() ) {
	    kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_qqbar_g3P2_oct(args);
    }
    if (fl[0].IsQuark() && fl[1].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3].IsGluon() && fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_qqbar_g3P2_oct(args);
    }
  }
  return NULL;
}

XS_qqbar_g3P2_oct::XS_qqbar_g3P2_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i<2 && fl[i].IsAnti())  m_a = i;
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P2_c_8_J_psi_1S)   LDME = 5.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_c_8_psi_2S)     LDME = 5.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_1S) LDME = 5.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_2S) LDME = 5.*0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_3S) LDME = 5.*0.001*pss["LDME_992200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  m_pref   = 8./(81.*pow(m_mass,3))*pow(4.*M_PI*m_alphaS,3);
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_qqbar_g3P2_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2();
  double t  = (mom[0]-mom[2]).Abs2();
  double u  = (mom[0]-mom[3]).Abs2();
  double heq0 = s*(sqr(t)+sqr(u));
  double heq1 = 36.*m_mass2*t*u;
  double heq2 = 18.*sqr(m_mass2)*(sqr(t)+sqr(u))/s;
  double dnom = pow(s-m_mass2,4);

  return CouplingFactor(3,0)*m_pref*(heq0+heq1+heq2)/dnom*LDME;
}

bool XS_qqbar_g3P2_oct::SetColours(const Vec4D_Vector& mom) 
{
  m_colours[1-m_a][0] = m_colours[3-m_a][0] = Flow::Counter();
  m_colours[m_a][1]   = m_colours[2+m_a][1] = Flow::Counter();
  m_colours[2][m_a]   = m_colours[3][1-m_a] = Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_g3P2_oct,"1XS_gg_g3P2_oct")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::XS_gg_g3P2_oct>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()==4 && args.m_orders[0]==3 && args.m_orders[1]==0) {
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[2].IsGluon()  &&
	 fl[3].IsOctetMeson() ) {
      kf_code kfc = fl[3].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_gg_g3P2_oct(args);
    }
    if ( fl[0].IsGluon() && fl[1].IsGluon() && fl[3].IsGluon()  &&
	 fl[2].IsOctetMeson() ) {
      kf_code kfc = fl[2].Kfcode();
      if (kfc==kf_3P2_c_8_J_psi_1S   || kfc==kf_3P2_c_8_psi_2S    ||
          kfc==kf_3P2_b_8_Upsilon_1S || kfc==kf_3P2_b_8_Upsilon_2S|| kfc==kf_3P2_b_8_Upsilon_3S)  return new XS_gg_g3P2_oct(args);
    }
  }
  return NULL;
}

XS_gg_g3P2_oct::XS_gg_g3P2_oct(const External_ME_Args& args):
  ME2_Base(args)
{
  m_oqcd = 3; m_oew = 0;
  const ATOOLS::Flavour_Vector& fl = args.Flavours();
  for (short int i=0;i<4;i++) {
    if (i>1 && fl[i].IsOctetMeson()) m_S = i;
    m_colours[i][0] = m_colours[i][1] = 0;
  }
  m_mass   = fl[m_S].Mass(true);  m_mass2 = sqr(m_mass); LDME = 0.;
  double mc = ATOOLS::Flavour(kf_c).Mass();
  double mb = ATOOLS::Flavour(kf_b).Mass();
  if (mc == 0) mc = 1.4;
  if (mb == 0) mb = 4.88;

  auto pss = Settings::GetMainSettings()["EXTRA_XS"];
  pss["LDME_992000443"].SetDefault(3E-03 * mc * mc*1000);
  pss["LDME_992100443"].SetDefault(0.005 * mc * mc*1000);
  pss["LDME_992000553"].SetDefault(5 * 0.0121 * mb * mb*1000);
  pss["LDME_992100553"].SetDefault(5 * 0.0141 * mb * mb*1000);
  pss["LDME_992200553"].SetDefault(5 * 0.0164 * mb * mb*1000);
  if (fl[m_S].Kfcode()==kf_3P2_c_8_J_psi_1S)   LDME = 5.*0.001*pss["LDME_992000443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_c_8_psi_2S)     LDME = 5.*0.001*pss["LDME_992100443"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_1S) LDME = 5.*0.001*pss["LDME_992000553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_2S) LDME = 5.*0.001*pss["LDME_992100553"].Get<double>(); 
  if (fl[m_S].Kfcode()==kf_3P2_b_8_Upsilon_3S) LDME = 5.*0.001*pss["LDME_992200553"].Get<double>(); 
  m_alphaS = MODEL::s_model->ScalarConstant("alpha_S");
  if (LDME == 0.0){
    throw std::runtime_error("\033[31mCross section is zero. Cannot continue.\033[0m");
  }
}

double XS_gg_g3P2_oct::operator()(const Vec4D_Vector& mom) 
{
  double s  = (mom[0]+mom[1]).Abs2(), sM2 = sqr(s-m_mass2); double s2 = sqr(s), s4 = sqr(s2), s6 = pow(s,6), s8 = sqr(s4), s10 = pow(s,10), s12 = pow(s,12);
  double t  = (mom[0]-mom[2]).Abs2(), tM2 = sqr(t-m_mass2);
  double u  = (mom[0]-mom[3]).Abs2(), uM2 = sqr(u-m_mass2);

  double M = m_mass, M2 = m_mass2;

  double z = sqrt(t*u), z2 = sqr(z), z4 = sqr(z2), z6 = pow(z,6), z8 = sqr(z4), z10 = pow(z,10), z12 = pow(z,12);
  // double heq0 = s*z2*(sqr(s2-z2)-2.*M2*s*z2-sqr(M2)*(s2+2.*z2)+pow(M,8))*(sqr(s2-z2)-M2*s*(2.*s2-z2)+sqr(M2)*s2);
  // double heq1l1 = M2*(2.*sqr(s2-z2)*(s6-4.*s4*z2+s2*z4-z6)-M2*s*(2.*s2-z2)*(5.*s6-17.*s4*z2+9.*s2*z4-z6)+sqr(M2)*(21.*s8-49.*s6*z2+21.*s4*z4-4.*s2*z6+z8));
  // double heq1l2 = M2*(-pow(M,6)*s*(24.*s6-30.*s4*z2+6.*s2*z4-z6)+pow(M,8)*s2*(16.*s4-9.*s2*z2+2*z4)-pow(M,10)*pow(s,3)*(6.*s2-z2)+pow(M,12)*s4);
  // double heq = heq0 + heq1l1 + heq1l2;
  // double dnom = pow(s-M2,4)*pow(s*M2+z2,4);

  double heq0 = s*z2*(s2*pow(s2-z2,4)-M2*s*z2*sqr(s2-z2)*(11.*s2+2.*z2)
                +sqr(M2)*(s8-12.*s6*z2+41.*s4*z4-20.*s2*z6+z8)
                -pow(M,6)*(s8-12.*s6*z2+41.*s4*z4-20.*s2*z6+z8)
                +pow(M,8)*(29.*s6-114.*s4*z2+108.*s2*z4-10.*z6)
                -pow(M,10)*s*(65.*s4-104.*s2*z2-33.*z4)+ pow(M,12)*(54.*s4-20.*s2*z2+7.*z4)
                -pow(M,14)*s*(23.*s2+5.*z2)+7.*pow(M,16)*s2);
  double heq1 = M2*(2.*s2*sqr(s2-z2)*(s6-4.*s4*z2+s2*z4-z6)
                  - M2*s*(10.*s10-37.*s8*z2+19.*s6*z4+11.*s4*z6-s2*z8-4*z10)
                  + sqr(M2)*(25.*s10-61.*s8*z2+27.*s6*z4-34.*s4*z6+23.*s2*z8-2.*z10)
                  - pow(M,6)*s*(42.*s8-77.*s6*z2+41.*s4*z4-22.*s2*z6+17.*s8)
                  + pow(M,8)*(53.*s8-88*s6*z2+69.*s4*z4-68.*s2*z6+3*z8)
                  - pow(M,10)*s*(54.*s6 - 85.*s4*z2 + 60.*s2*z4 - 9.*z6) + pow(M,12)*s2*(43.*s4-47.*s2*z2+20*z4)
                  - pow(M,14)*pow(s,3)*(22.*s2-9.*z2)+5.*pow(M,16)*s4);
  double heq2 = sqr(M2)*(2.*s2*(s12-8.*s10*z2+22.*s8*z4-24.*s6*z6+10.*s4*z8-3.*s2*z10+z12)
                      - M2*s*(16.*s12-102.*s10*z2+210.*s8*z4-153.*s6*z6+36.*s4*z8-6.*s2*z10+4.*z12)
                      + sqr(M2)*(60.*s12-306.*s10*z2+482*s8*z4-271.*s6*z6+77.*s4*z8-18.*s2*z10+2.*z12)
                      - pow(M,6)*s*(140.*s10-573*s8*z2+710.*s6*z4-344.*s4*z6+91.*s2*z8-18.*z10)
                      + pow(M,8)*(226.*s10-741.*s8*z2+737.*s6*z4-310.*s4*z6+77.*s2*z8-4.*z10)
                      - pow(M,10)*s*(264.*s8-686.*s6*z2+541.*s4*z4-177.*s2*z6+25.*z8)
                      + pow(M,12)*(226.*s8-452*s6*z2+261*s4*z4-55.*s2*z6+2*z8)
                      - pow(M,14)*(140.*s6-201.*s4*z2+71.*s2*z4-6.*z6)
                      + pow(M,16)*s2*(60.*s4-53*s2*z2+8.*z4)-2.*pow(M,18)*pow(s,3)*(8.*s2-3.*z2)+2.*pow(M,20)*s4)/(s*z2);
  double dnom = pow(s-M2,6)*pow(s*M2+z2,4);
  m_pref = pow(4.*M_PI*m_alphaS,3)*CouplingFactor(3,0)/(6.*pow(m_mass,3));
  //  std::cout<<" ----- gg-3P2 -----" << LDME*m_pref*(heq0+3*heq1+3*heq2)/dnom << std::endl;
  return LDME*m_pref*(heq0+3*heq1+3*heq2)/dnom;
}

bool XS_gg_g3P2_oct::SetColours(const Vec4D_Vector& mom) 
{
  size_t bit = ran->Get()<0.5 ? 0 : 1;
  size_t cross = ran->Get()<0.5 ? 0 : 1;
  if (cross==0) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[1][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[m_S][1-bit] = m_colours[5-m_S][bit] = Flow::Counter();
  }
  if (cross==1) {
    m_colours[0][bit] = m_colours[m_S][bit] = Flow::Counter();
    m_colours[0][1-bit] = m_colours[5-m_S][1-bit] = Flow::Counter();
    m_colours[1][bit] = m_colours[5-m_S][bit] = Flow::Counter();
    m_colours[1][1-bit] = m_colours[m_S][1-bit] = Flow::Counter();
  }
  return true;
}
