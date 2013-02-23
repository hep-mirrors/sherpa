#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {
  class XS_f1f1_f1f1: public ME2_Base {  
  private:
    bool   m_Z_on, m_P_on;
    int    m_anti;
    double m_mz2, m_wz2;
    double m_sin2tw, m_cos2tw, m_eq, m_y3f, m_v, m_a;
    double m_aqed, m_pref_qed, m_pref_Z, m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f1_f1f1(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f1_f1f1,"0XS_f1f1_f1f1")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f1_f1f1>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  if (fl[0]!=fl[1] || fl[0]!=fl[2] || fl[0]!=fl[3]) return NULL;

  if (!(fl[0].IsFermion() && fl[1].IsFermion() && 
	fl[2].IsFermion() && fl[3].IsFermion())) return NULL;

  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!(ATOOLS::Flavour(kf_Z).IsOn() ||
	ATOOLS::Flavour(kf_photon).IsOn())) return NULL;
  if ((fl[0].Charge()!=0. && ATOOLS::Flavour(kf_photon).IsOn()) ||
      ATOOLS::Flavour(kf_Z).IsOn()) return new XS_f1f1_f1f1(pi,fl);
  return NULL;
}


XS_f1f1_f1f1::XS_f1f1_f1f1(const Process_Info& pi, const Flavour_Vector& fl):
  ME2_Base(pi, fl),
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_anti(int(fl[0].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq(fl[0].Charge()),
  m_y3f((2.*int(fl[0].IsUptype())-1)/2.),
  m_v(m_y3f-2.*m_eq*m_sin2tw), m_a(m_y3f),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
}

double XS_f1f1_f1f1::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  M_t = 0., M_u = 0., M_mix = 0.;
  if (m_P_on) {
    M_t   +=     sqr(m_pref_qed*m_eq*m_eq)    * (s*s+u*u)/(t*t);
    M_mix +=  2.*sqr(m_pref_qed*m_eq*m_eq)/3. * (s*s)/(t*u);
    M_u   +=     sqr(m_pref_qed*m_eq*m_eq)    * (s*s+t*t)/(u*u); 
  }
  if (m_Z_on) {
    M_t   +=     sqr(m_pref_Z)/((sqr(t-m_mz2)+m_mz2*m_wz2)) *
      (sqr(sqr(m_v)+sqr(m_a)) * (u*u+s*s) +
       4.*sqr(m_a*m_v) * (s*s-u*u));
    M_mix +=  2.*sqr(m_pref_Z)/3.  *
      ((t-m_mz2)*(u-m_mz2)-m_mz2*m_wz2)/(sqr((t-m_mz2)*(u-m_mz2)-m_mz2*m_wz2)+m_mz2*m_wz2*(t+s-2.*m_mz2)) *
      (sqr(sqr(m_v)+sqr(m_a))+ 4.*sqr(m_a*m_v)) * s*s;
    M_u   +=     sqr(m_pref_Z)/((sqr(u-m_mz2)+m_mz2*m_wz2)) *
      (sqr(sqr(m_v)+sqr(m_a)) * (s*s+t*t) +
       4.*sqr(m_a*m_v) * (s*s-t*t));
  }
  if (m_P_on && m_Z_on) {
    M_t   +=  m_pref_qed*m_pref_Z  *                   
      (t-m_mz2)/(t*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_v*m_eq) * (s*s+u*u) + sqr(m_a*m_eq) * (s*s-u*u));
    M_mix += 2.*m_pref_qed*m_pref_Z/3. * 
      ((u-m_mz2)/(t*(sqr(u-m_mz2)+m_mz2*m_wz2)) + 
       (t-m_mz2)/(u*(sqr(t-m_mz2)+m_mz2*m_wz2)) )  * 
      sqr(m_eq)*(sqr(m_v)+sqr(m_a)) * s*s;
    M_u   +=  m_pref_qed*m_pref_Z   * 
      (u-m_mz2)/(u*(sqr(u-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_v*m_eq) * (s*s+t*t) + sqr(m_a*m_eq) * (s*s-t*t));
  }
  return (M_t + M_u + M_mix)*CouplingFactor(0,2);
}

bool XS_f1f1_f1f1::SetColours(const Vec4D_Vector& mom) 
{
  (*this)(mom);
  if (SetColours()) {}
  //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
  //else
  //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
  return true;
}


bool XS_f1f1_f1f1::SetColours() 
{
  if (M_t > (M_t+M_u) * ran->Get()) {
    p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[3][m_anti] = p_colours[1][m_anti] = Flow::Counter();
    return true;
  }
  p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
  p_colours[2][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  return false;
}

double XS_f1f1_f1f1::KFactor(double scale) 
{ 
  return 1.; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f1b_f1f1b: public ME2_Base {  
  private:
    bool   m_Z_on, m_P_on;
    int    m_anti1, m_anti2;
    double m_mz2, m_wz2;
    double m_sin2tw, m_cos2tw, m_eq, m_y3f, m_v, m_a;
    double m_aqed, m_pref_qed, m_pref_Z, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    XS_f1f1b_f1f1b(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f1b_f1f1b,"0XS_f1f1b_f1f1b")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f1b_f1f1b>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  if (fl[0]!=fl[1].Bar() || 
      fl[0]!=fl[2] || fl[1]!=fl[3]) return NULL;
  
  if (!(fl[0].IsFermion() && fl[1].IsFermion() && 
	fl[2].IsFermion() && fl[3].IsFermion())) return NULL;
  
  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!(ATOOLS::Flavour(kf_Z).IsOn() ||
	ATOOLS::Flavour(kf_photon).IsOn())) return NULL;
  if ((fl[0].Charge()!=0. && ATOOLS::Flavour(kf_photon).IsOn()) ||
      ATOOLS::Flavour(kf_Z).IsOn()) return new XS_f1f1b_f1f1b(pi,fl);
  return NULL;
}

XS_f1f1b_f1f1b::XS_f1f1b_f1f1b(const Process_Info& pi, const Flavour_Vector& fl):
  ME2_Base(pi, fl), 
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_anti1(int(fl[0].IsAnti())),m_anti2(int(fl[2].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq(fl[0].Charge()),
  m_y3f((2.*int(fl[0].IsUptype())-1)/2.),
  m_v(m_y3f-2.*m_eq*m_sin2tw), m_a(m_y3f),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
}

double XS_f1f1b_f1f1b::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) { double help = u; u = t; t = help; }
  M_t = 0., M_s = 0., M_mix = 0.;
  if (m_P_on) {
    M_t   +=     sqr(m_pref_qed*m_eq*m_eq)    * (u*u+s*s)/(t*t);
    M_mix +=  2.*sqr(m_pref_qed*m_eq*m_eq)/3. * (u*u)/(t*s);
    M_s   +=     sqr(m_pref_qed*m_eq*m_eq)    * (u*u+t*t)/(s*s); 
  }
  if (m_Z_on) {
    M_s   +=     sqr(m_pref_Z)/((sqr(s-m_mz2)+m_mz2*m_wz2)) *
      ((sqr(m_v)+sqr(m_a)) * (sqr(m_v)+sqr(m_a)) * (u*u+t*t) +
       4.*(m_a*m_v*m_a*m_v) * (u*u-t*t));
    M_mix +=  2.*sqr(m_pref_Z)/3.  *
      ((t-m_mz2)*(s-m_mz2)-m_mz2*m_wz2)/(sqr((t-m_mz2)*(s-m_mz2)-m_mz2*m_wz2)+m_mz2*m_wz2*(t+s-2.*m_mz2)) *
      (sqr(m_v*m_v+m_a*m_a) + 4.*m_a*m_a*m_v*m_v) * u*u;
    M_t   +=     sqr(m_pref_Z)/((sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ((sqr(m_v)+sqr(m_a)) * (sqr(m_v)+sqr(m_a)) * (u*u+s*s) +
       4.*(m_a*m_v*m_a*m_v) * (u*u-s*s));
  }
  if (m_P_on && m_Z_on) {
    M_t   +=  m_pref_qed*m_pref_Z  *                   
      (t-m_mz2)/(t*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_v*m_eq) * (s*s+u*u) + sqr(m_a*m_eq) * (s*s-u*u));
    M_mix += 2.*m_pref_qed*m_pref_Z/3. * 
      ((s-m_mz2)/(t*(sqr(s-m_mz2)+m_mz2*m_wz2)) + 
       (t-m_mz2)/(s*(sqr(t-m_mz2)+m_mz2*m_wz2)) )  * 
      sqr(m_eq)*(sqr(m_v)+sqr(m_a)) * u*u;
    M_s   +=  m_pref_qed*m_pref_Z   * 
      (s-m_mz2)/(s*(sqr(s-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_v*m_eq) * (u*u+t*t) + sqr(m_a*m_eq) * (u*u-t*t));
  }
  return 2.*(M_t + M_s + M_mix)*CouplingFactor(0,2);
}

bool XS_f1f1b_f1f1b::SetColours(const Vec4D_Vector& mom) 
{
  (*this)(mom);
  if (SetColours()) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) {}
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
  }
  else {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = s;
  }
  return true;
}


bool XS_f1f1b_f1f1b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran->Get()) {
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) {
      p_colours[3][1-m_anti2] = p_colours[0][m_anti1] = Flow::Counter();
      p_colours[2][m_anti2]   = p_colours[1][1-m_anti1] = Flow::Counter();
    }
    else {
      p_colours[2][m_anti2]   = p_colours[0][m_anti1] = Flow::Counter();
      p_colours[3][1-m_anti2] = p_colours[1][1-m_anti1] = Flow::Counter();
    }
    return true;
  }
  p_colours[0][m_anti1] = p_colours[1][1-m_anti1] = Flow::Counter();
  p_colours[2][m_anti2] = p_colours[3][1-m_anti2] = Flow::Counter();
  return false;
}

double XS_f1f1b_f1f1b::KFactor(double scale) 
{ 
  return 1.; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f1b_f2f2b: public ME2_Base {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti1, m_anti2;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    // constructor
    XS_f1f1b_f2f2b(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f1b_f2f2b,"0XS_f1f1b_f2f2b")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f1b_f2f2b>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  kf_code  kfc1  = fl[0].Kfcode(), kfc2  = fl[2].Kfcode();
  if (fl[0]!=fl[1].Bar() || 
      fl[0]==fl[2] || fl[0]==fl[3] || 
      fl[2]!=fl[3].Bar()) return NULL;
  
  if (!(fl[0].IsFermion() && fl[1].IsFermion() && 
        fl[2].IsFermion() && fl[3].IsFermion())) return NULL;
  
  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if ((fl[0].Charge()!=0. && fl[2].Charge()!=0. &&
       ATOOLS::Flavour(kf_photon).IsOn()) ||
      (ATOOLS::Flavour(kf_Wplus).IsOn() && 
       fl[0].IsQuark() && fl[2].IsQuark() && 
       ((kfc1%2==0 && kfc2%2!=0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2))>0) ||
        (kfc1%2!=0 && kfc2%2==0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2))>0))) ||
      ATOOLS::Flavour(kf_Z).IsOn() )
    return new XS_f1f1b_f2f2b(pi,fl);
  return NULL;
}


XS_f1f1b_f2f2b::XS_f1f1b_f2f2b(const Process_Info& pi, const Flavour_Vector& fl) :
  ME2_Base(pi, fl), 
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_W_on(ATOOLS::Flavour(kf_Wplus).IsOn()),
  m_anti1(int(fl[0].IsAnti())),m_anti2(int(fl[2].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(fl[0].Charge()),
  m_eq2(fl[2].Charge()),
  m_y3f1((2.*int(fl[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(fl[2].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[2].Kfcode();
    if (fl[0].IsQuark() && fl[2].IsQuark()) { 
      if (kfc1%2==0 && kfc2%2!=0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2);
      if (kfc1%2!=0 && kfc2%2==0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
}

double XS_f1f1b_f2f2b::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) { double help = u; u = t; t = help; }
  M_s = M_t = M_mix = 0.;
  if (m_P_on) {
    M_s   +=     sqr(m_pref_qed*m_eq1*m_eq2)   *  (u*u+t*t)/(s*s); 
  }
  if (m_Z_on) {
    M_s   +=     sqr(m_pref_Z)/((sqr(s-m_mz2)+m_mz2*m_wz2)) *
      ((sqr(m_v1)+sqr(m_a1)) * (sqr(m_v2)+sqr(m_a2)) * (u*u+t*t) +
       4.*(m_a1*m_v1*m_a2*m_v2) * (u*u-t*t));
  }
  if (m_P_on && m_Z_on) {
    M_s   +=  m_pref_qed*m_pref_Z * 
      (s-m_mz2)/(s*(sqr(s-m_mz2)+m_mz2*m_wz2)) *
      (m_eq1*m_eq2*m_v1*m_v2  * (u*u+t*t) + m_eq1*m_eq2*m_a1*m_a2  * (u*u-t*t));
  }
  if (m_W_on) {
    M_t   +=     sqr(m_pref_W)/((sqr(t-m_mw2)+m_mw2*m_ww2)) * (2.*u*u);
  }
  if (m_W_on && m_P_on) {
    M_mix +=  2.*m_pref_qed*m_pref_W * (m_eq1*m_eq2)/3.   *
      (t-m_mw2)/(s*(sqr(t-m_mw2)+m_mw2*m_ww2)) * u*u;
  }
  if (m_W_on && m_Z_on) {
    double mixed = sqrt(m_mz2*m_wz2*m_mw2*m_ww2);
    M_mix +=  2.*m_pref_W*m_pref_Z /3.   *
      ((s-m_mz2)*(t-m_mw2)-mixed)/
      (sqr((s-m_mz2)*(t-m_mw2)-mixed)+sqr(sqrt(m_mz2*m_wz2)*(t-m_mw2)+sqrt(m_mw2*m_ww2)*(s-m_mz2))) *
      2.*(m_v1*m_v2+m_a1*m_a2) * u*u;
  }
  return 2.*(M_s+M_t+M_mix)*CouplingFactor(0,2);
}

bool XS_f1f1b_f2f2b::SetColours(const Vec4D_Vector& mom) 
{
  (*this)(mom);
  if (SetColours()) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) {}
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
  }
  else {} 
  //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = s;
  return true;
}


bool XS_f1f1b_f2f2b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran->Get()) {
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) {
      p_colours[3][1-m_anti2] = p_colours[0][m_anti1] = Flow::Counter();
      p_colours[2][m_anti2]   = p_colours[1][1-m_anti1] = Flow::Counter();
    }
    else {
      p_colours[2][m_anti2]   = p_colours[0][m_anti1] = Flow::Counter();
      p_colours[3][1-m_anti2] = p_colours[1][1-m_anti1] = Flow::Counter();
    }
    return true;
  }
  p_colours[0][m_anti1] = p_colours[1][1-m_anti1] = Flow::Counter();
  p_colours[2][m_anti2] = p_colours[3][1-m_anti2] = Flow::Counter();
  return false;
}

double XS_f1f1b_f2f2b::KFactor(double scale) 
{ 
  return 1.; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f2_f1f2: public ME2_Base {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti, m_rev;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2_f1f2(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f2_f1f2,"0XS_f1f2_f1f2")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f2_f1f2>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  if (!((fl[0]==fl[2] && fl[1]==fl[3]) || 
        (fl[0]==fl[3] && fl[1]==fl[2])) || 
      (fl[0].IsAnti() && !fl[1].IsAnti()) ||
      (!fl[0].IsAnti() && fl[1].IsAnti())) return NULL;

  if (!(fl[0].IsFermion() && fl[1].IsFermion() && 
        fl[2].IsFermion() && fl[3].IsFermion())) return NULL;

  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!(ATOOLS::Flavour(kf_Z).IsOn() ||
        ATOOLS::Flavour(kf_photon).IsOn() ||
        ATOOLS::Flavour(kf_Wplus).IsOn())) return NULL;
  if ((fl[0].Charge()!=0. && fl[1].Charge()!=0. &&
       ATOOLS::Flavour(kf_photon).IsOn()) ||
      (ATOOLS::Flavour(kf_Wplus).IsOn() && 
       fl[0].IsQuark() && fl[1].IsQuark() && 
       ((kfc1%2==0 && kfc2%2!=0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2))>0) ||
        (kfc1%2!=0 && kfc2%2==0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2))>0))) ||
      ATOOLS::Flavour(kf_Z).IsOn()) return new XS_f1f2_f1f2(pi,fl);
  return NULL;
}


XS_f1f2_f1f2::XS_f1f2_f1f2(const Process_Info& pi, const Flavour_Vector& fl) :
  ME2_Base(pi, fl), 
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_W_on(ATOOLS::Flavour(kf_Wplus).IsOn()),
  m_anti(int(fl[0].IsAnti())),m_rev(int(fl[0]!=fl[2])),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(std::pow(-1.,m_anti) * fl[0].Charge()),
  m_eq2(std::pow(-1.,m_anti) * fl[1].Charge()),
  m_y3f1((2.*int(fl[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(fl[1].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
    if (fl[0].IsQuark() && fl[1].IsQuark()) { 
      if (kfc1%2==0 && kfc2%2!=0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2);
      if (kfc1%2!=0 && kfc2%2==0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
}

double XS_f1f2_f1f2::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if (m_rev) { double help = u; u = t; t = help; }
  M_t = M_u = M_mix = 0.;
  if (m_P_on) {
    M_t   +=     sqr(m_pref_qed*m_eq1*m_eq2)   *  (u*u+s*s)/(t*t); 
  }
  if (m_Z_on) {
    M_t   +=     sqr(m_pref_Z)/((sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ((sqr(m_v1)+sqr(m_a1)) * (sqr(m_v2)+sqr(m_a2)) * (u*u+s*s) +
       4.*(m_a1*m_v1*m_a2*m_v2) * (s*s-u*u));
  }
  if (m_P_on && m_Z_on) {
    M_t   +=  m_pref_qed*m_pref_Z * 
      (t-m_mz2)/(t*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      (m_eq1*m_eq2*m_v1*m_v2  * (s*s+u*u) + m_eq1*m_eq2*m_a1*m_a2  * (s*s-u*u));
  }
  if (m_W_on) {
    M_u   +=     sqr(m_pref_W)/((sqr(u-m_mw2)+m_mw2*m_ww2)) * (2.*s*s);
  }
  if (m_W_on && m_P_on) {
    M_mix +=  2.*m_pref_qed*m_pref_W * (m_eq1*m_eq2)/3.   *
      (u-m_mw2)/(t*(sqr(u-m_mw2)+m_mw2*m_ww2)) * s*s;
  }
  if (m_W_on && m_Z_on) {
    double mixed = sqrt(m_mz2*m_wz2*m_mw2*m_ww2);
    M_mix +=  2.*m_pref_W*m_pref_Z /3.   *
      ((t-m_mz2)*(u-m_mw2)-mixed)/
      (sqr((t-m_mz2)*(u-m_mw2)-mixed)+sqr(sqrt(m_mz2*m_wz2)*(u-m_mw2)+sqrt(m_mw2*m_ww2)*(t-m_mz2))) *
      2.*(m_v1*m_v2+m_a1*m_a2) * s*s;
  }
  return 2.*(M_t+M_u+M_mix)*CouplingFactor(0,2);
}

bool XS_f1f2_f1f2::SetColours(const Vec4D_Vector& mom) 
{
  (*this)(mom);
  if (SetColours()) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    if (m_rev) {
      //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
    }
  }
  else {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
    if (m_rev) {
      //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    }
  }
  return true;
}


bool XS_f1f2_f1f2::SetColours() 
{
  if (M_t > (M_t+M_u) * ran->Get()) {
    if (m_rev) {
      p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
      p_colours[2][m_anti] = p_colours[1][m_anti] = Flow::Counter();
    }
    else {
      p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
      p_colours[3][m_anti] = p_colours[1][m_anti] = Flow::Counter();
    }
    return true;
  }
  if (!m_rev) {
    p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[2][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  else {
    p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[3][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  return false;
}

double XS_f1f2_f1f2::KFactor(double scale) 
{ 
  return 1.; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f2b_f1f2b: public ME2_Base {  
  private:
    bool   m_Z_on, m_P_on, m_W_on;
    int    m_anti, m_rev;
    double m_mz2, m_wz2, m_mw2, m_ww2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    Complex m_ckm;
    double m_aqed, m_pref_qed, m_pref_Z, m_pref_W, m_kappa, m_colfac;
    double M_t, M_s, M_mix;
  public:
    // constructor
    XS_f1f2b_f1f2b(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f2b_f1f2b,"0XS_f1f2b_f1f2b")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f2b_f1f2b>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  if (!((fl[0]==fl[2] && fl[1]==fl[3]) || 
        (fl[0]==fl[3] && fl[1]==fl[2])) || 
      (fl[0].IsAnti() && fl[1].IsAnti()) ||
      (!fl[0].IsAnti() && !fl[1].IsAnti())) return NULL;

  if (!(fl[0].IsFermion() && fl[1].IsFermion() && 
        fl[2].IsFermion() && fl[3].IsFermion())) return NULL;

  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!(ATOOLS::Flavour(kf_Z).IsOn() ||
        ATOOLS::Flavour(kf_photon).IsOn() ||
        ATOOLS::Flavour(kf_Wplus).IsOn())) return NULL;
  if ((fl[0].Charge()!=0. && fl[1].Charge()!=0. &&
       ATOOLS::Flavour(kf_photon).IsOn()) ||
      (ATOOLS::Flavour(kf_Wplus).IsOn() && 
       fl[0].IsQuark() && fl[1].IsQuark() && 
       ((kfc1%2==0 && kfc2%2!=0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2))>0) ||
        (kfc1%2!=0 && kfc2%2==0 && 
         abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2))>0))) ||
      ATOOLS::Flavour(kf_Z).IsOn()) return new XS_f1f2b_f1f2b(pi,fl);
  return NULL;
}

XS_f1f2b_f1f2b::XS_f1f2b_f1f2b(const Process_Info& pi, const Flavour_Vector& fl) :
  ME2_Base(pi, fl), 
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_W_on(ATOOLS::Flavour(kf_Wplus).IsOn()),
  m_anti(int(fl[0].IsAnti())),m_rev(int(fl[0]!=fl[2])),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(fl[0].Charge()),
  m_eq2(fl[1].Charge()),
  m_y3f1((2.*int(fl[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(fl[1].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
    if (fl[0].IsQuark() && fl[1].IsQuark()) { 
      if (kfc1%2==0 && kfc2%2!=0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2);
      if (kfc1%2!=0 && kfc2%2==0)
        m_ckm = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
}

double XS_f1f2b_f1f2b::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if (m_rev) { double help = u; u = t; t = help; }
  M_t = M_s = M_mix = 0.;
  if (m_P_on) {
    M_t   +=     sqr(m_pref_qed*m_eq1*m_eq2)   *  (u*u+s*s)/(t*t); 
  }
  if (m_Z_on) {
    M_t   +=     sqr(m_pref_Z)/((sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ((sqr(m_v1)+sqr(m_a1)) * (sqr(m_v2)+sqr(m_a2)) * (u*u+s*s) +
       4.*(m_a1*m_v1*m_a2*m_v2) * (s*s-u*u));
  }
  if (m_P_on && m_Z_on) {
    M_t   +=  m_pref_qed*m_pref_Z * 
      (t-m_mz2)/(t*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      (m_eq1*m_eq2*m_v1*m_v2  * (s*s+u*u) + m_eq1*m_eq2*m_a1*m_a2  * (s*s-u*u));
  }
  if (m_W_on) {
    M_s   +=     sqr(m_pref_W)/((sqr(s-m_mw2)+m_mw2*m_ww2)) * (2.*u*u);
  }
  if (m_W_on && m_P_on) {
    M_mix +=  2.*m_pref_qed*m_pref_W * (m_eq1*m_eq2)/3.   *
      (s-m_mw2)/(t*(sqr(s-m_mw2)+m_mw2*m_ww2)) * u*u;
  }
  if (m_W_on && m_Z_on) {
    double mixed = sqrt(m_mz2*m_wz2*m_mw2*m_ww2);
    M_mix +=  2.*m_pref_W*m_pref_Z /3.   *
      ((t-m_mz2)*(s-m_mw2)-mixed)/
      (sqr((t-m_mz2)*(s-m_mw2)-mixed)+sqr(sqrt(m_mz2*m_wz2)*(s-m_mw2)+sqrt(m_mw2*m_ww2)*(t-m_mz2))) *
      2.*(m_v1*m_v2+m_a1*m_a2) * u*u;
  }
  return 2.*(M_t+M_s+M_mix)*CouplingFactor(0,2);
}

bool XS_f1f2b_f1f2b::SetColours(const Vec4D_Vector& mom) 
{
  (*this)(mom);
  if (SetColours()) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    if (m_rev) {
      //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
    }
  }
  else {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = s;
  }
  return true;
}


bool XS_f1f2b_f1f2b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran->Get()) {
    if (m_rev) {
      p_colours[3][m_anti]   = p_colours[0][m_anti]   = Flow::Counter();
      p_colours[2][1-m_anti] = p_colours[1][1-m_anti] = Flow::Counter();
    }
    else {
      p_colours[2][m_anti]   = p_colours[0][m_anti]   = Flow::Counter();
      p_colours[3][1-m_anti] = p_colours[1][1-m_anti] = Flow::Counter();
    }
    return true;
  }
  p_colours[0][m_anti]       = p_colours[1][1-m_anti]       = Flow::Counter();
  p_colours[2+m_rev][m_anti] = p_colours[3-m_rev][1-m_anti] = Flow::Counter();
  return false;
}

double XS_f1f2b_f1f2b::KFactor(double scale) 
{ 
  return 1.; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f2_f3f4: public ME2_Base {  
  private:
    int    m_anti, m_rev;
    double m_mw2, m_ww2, m_pref_W;
    Complex m_ckm1, m_ckm2;
    double m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2_f3f4(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f2_f3f4,"0XS_f1f2_f3f4")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f2_f3f4>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  kf_code kfc3 = fl[2].Kfcode(), kfc4 = fl[3].Kfcode();
  if (!(fl[0].IsQuark() &&fl[1].IsQuark() &&
        fl[2].IsQuark() &&fl[3].IsQuark())) return NULL;
  if (!(((fl[0].IsUptype() && fl[1].IsDowntype()) ||
         (fl[0].IsDowntype() && fl[1].IsUptype())) &&
        ((fl[2].IsUptype() && fl[3].IsDowntype()) ||
         (fl[2].IsDowntype() && fl[3].IsUptype())))  ||
      (!(!fl[0].IsAnti() && !fl[1].IsAnti() && 
         !fl[2].IsAnti() && !fl[3].IsAnti()) &&
       !(fl[0].IsAnti() && fl[1].IsAnti() && 
         fl[2].IsAnti() && fl[3].IsAnti()))) return NULL;
  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!ATOOLS::Flavour(kf_Wplus).IsOn()) return NULL;
  if (fl[0].IsUptype()) {
    if (fl[2].IsDowntype() && 
        (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc3/2))>0) ||
         !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc2/2))>0)))
      return NULL;
    if (fl[2].IsUptype() &&
        (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc4/2))>0) ||
         !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc2/2))>0)))
      return NULL;
  }
  else {
    if (fl[2].IsDowntype() && 
        (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc1/2))>0) ||
         !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc3/2))>0)))
      return NULL;
    if (fl[2].IsUptype() &&
        (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc1/2))>0) ||
         !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc4/2))>0)))
      return NULL;
  }
  return new XS_f1f2_f3f4(pi,fl);
}

XS_f1f2_f3f4::XS_f1f2_f3f4(const Process_Info& pi, const Flavour_Vector& fl) :
  ME2_Base(pi,fl), 
  m_anti(int(fl[0].IsAnti())), m_rev(false),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Width())),
  m_pref_W((4.*M_PI*MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms())))/
           (2.*MODEL::s_model->ScalarConstant(std::string("sin2_thetaW")))),
  m_ckm1(Complex(0.,0.)), m_ckm2(Complex(0.,0.))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  kf_code kfc3 = fl[2].Kfcode(), kfc4 = fl[3].Kfcode();

  if (fl[0].IsUptype()) {
    if (fl[2].IsDowntype()) {
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc3/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc2/2);
    } 
    if (fl[2].IsUptype()) {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc4/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc2/2);
    }
  }
  else {
    if (fl[2].IsDowntype()) {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc3/2);
    }
    if (fl[2].IsUptype()) {
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc4/2);
    }
  }
}

double XS_f1f2_f3f4::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if (m_rev) { double help = u; u = t; t = help; }
  return sqr(m_pref_W)*CouplingFactor(0,2)/((sqr(t-m_mw2)+m_mw2*m_ww2)) * (s*s);
}

bool XS_f1f2_f3f4::SetColours(const Vec4D_Vector& mom) 
{
  if (m_rev) {
    p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[2][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  else {
    p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[3][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  if (m_rev) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
  }
  else {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
  }
  return true;
}


bool XS_f1f2_f3f4::SetColours() 
{
  return true;
}

double XS_f1f2_f3f4::KFactor(double scale) 
{ 
  return 1.; 
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace EXTRAXS {
  class XS_f1f2b_f3f4b: public ME2_Base {  
  private:
    int    m_anti, m_rev, m_schannel;
    double m_mw2, m_ww2, m_pref_W;
    Complex m_ckm1, m_ckm2;
    double m_kappa, m_colfac;
    double M_t, M_u, M_mix;
  public:
    // constructor
    XS_f1f2b_f3f4b(const Process_Info& pi, const Flavour_Vector& fl);
    
    // member functions
    double operator()(const Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   SetColours();

    double KFactor(const double scale);
  };
}

DECLARE_TREEME2_GETTER(XS_f1f2b_f3f4b,"0XS_f1f2b_f3f4b")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_f1f2b_f3f4b>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;

  if (MODEL::s_model->Name()=="QCD") return NULL;
  if (pi.m_oqcd!=0 || pi.m_oew!=2) return NULL;
  if (!ATOOLS::Flavour(kf_Wplus).IsOn()) return NULL;
  if (fl[0].IntCharge()+fl[1].IntCharge()!= fl[2].IntCharge()+fl[3].IntCharge())
    return NULL;
  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  kf_code kfc3 = fl[2].Kfcode(), kfc4 = fl[3].Kfcode();
  if (!(fl[0].IsQuark() &&fl[1].IsQuark() &&
        fl[2].IsQuark() &&fl[3].IsQuark())) return NULL;
  if (fl[0].IsUptype() && fl[1].IsDowntype()) {
    if ((fl[2].IsUptype() && 
         (!fl[3].IsDowntype() ||
          (fl[0].IsAnti()^fl[2].IsAnti()) || 
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc4/2))>0))) ||
        abs(fl[0].IntCharge()+fl[1].IntCharge())!=3) return NULL;
    if ((fl[3].IsUptype() && 
         (!fl[2].IsDowntype() ||
          (fl[0].IsAnti()^fl[3].IsAnti()) || 
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc3/2))>0))) ||
        abs(fl[0].IntCharge()+fl[1].IntCharge())!=3) return NULL;
  }
  if (fl[0].IsDowntype() && fl[1].IsUptype()) {
    if ((fl[2].IsDowntype() && 
         (!fl[3].IsUptype() ||
          (fl[0].IsAnti()^fl[2].IsAnti()) || 
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc3/2))>0))) ||
        abs(fl[0].IntCharge()+fl[1].IntCharge())!=3) return NULL;
    if ((fl[3].IsDowntype() && 
         (!fl[2].IsUptype() ||
          (fl[0].IsAnti()^fl[3].IsAnti()) || 
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc4/2))>0))) ||
        abs(fl[0].IntCharge()+fl[1].IntCharge())!=3) return NULL;
  }
  if (fl[0].IsUptype() && fl[1].IsUptype()) {
    if (!fl[2].IsDowntype() || !fl[3].IsDowntype() ||
        !(fl[0].IsAnti()^fl[1].IsAnti()) ||
        !(fl[2].IsAnti()^fl[3].IsAnti())) return NULL;
    if ((!(fl[0].IsAnti()^fl[2].IsAnti()) &&
         (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc3/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc4/2))>0))) ||
        abs(fl[0].IntCharge()-fl[2].IntCharge())!=3) return NULL;
    if ((!(fl[0].IsAnti()^fl[3].IsAnti()) &&
         (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc4/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc3/2))>0))) ||
        abs(fl[0].IntCharge()-fl[3].IntCharge())!=3) return NULL;
  }
  if (fl[0].IsDowntype() && fl[1].IsDowntype()) {
    if (!fl[2].IsUptype() || !fl[3].IsUptype() ||
        !(fl[0].IsAnti()^fl[1].IsAnti()) ||
        !(fl[2].IsAnti()^fl[3].IsAnti())) return NULL;
    if ((!(fl[0].IsAnti()^fl[2].IsAnti()) &&
         (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc1/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc2/2))>0))) ||
        abs(fl[0].IntCharge()-fl[2].IntCharge())!=3) return NULL;
    if ((!(fl[0].IsAnti()^fl[3].IsAnti()) &&
         (!(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc1/2))>0) ||
          !(abs(MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc2/2))>0))) ||
        abs(fl[0].IntCharge()-fl[3].IntCharge())!=3) return NULL;
  }
  return new XS_f1f2b_f3f4b(pi,fl); 
}

XS_f1f2b_f3f4b::XS_f1f2b_f3f4b(const Process_Info& pi, const Flavour_Vector& fl) :
  ME2_Base(pi,fl), 
  m_anti(int(fl[0].IsAnti())), m_rev(false), m_schannel(true),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(kf_Wplus).Width())),
  m_pref_W((4.*M_PI*MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms())))/
           (2.*MODEL::s_model->ScalarConstant(std::string("sin2_thetaW")))),
  m_ckm1(Complex(0.,0.)), m_ckm2(Complex(0.,0.))
{
  //for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  kf_code kfc1 = fl[0].Kfcode(), kfc2 = fl[1].Kfcode();
  kf_code kfc3 = fl[2].Kfcode(), kfc4 = fl[3].Kfcode();

  if (fl[0].IsUptype() && fl[1].IsDowntype()) {
    if (fl[2].IsUptype()) {
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc4/2);
    }
    else {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc2/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc3/2);
    }
    //p_colours[1][1-m_anti]       = p_colours[0][m_anti]       = 500;
    //p_colours[3-m_rev][1-m_anti] = p_colours[2+m_rev][m_anti] = 501;
  }
  else if (fl[0].IsDowntype() && fl[1].IsUptype()) {
    if (fl[2].IsDowntype()) {
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc3/2);
    }
    else {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc4/2);
    }
    //p_colours[1][1-m_anti]       = p_colours[0][m_anti]       = 500;
    //p_colours[3-m_rev][1-m_anti] = p_colours[2+m_rev][m_anti] = 501;
  }
  else if (fl[0].IsUptype() && fl[1].IsUptype()) {
    m_schannel = false;
    if (!(fl[0].IsAnti()^fl[2].IsAnti())) {  
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc3/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc4/2);
    }
    if (!(fl[0].IsAnti()^fl[3].IsAnti())) {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc1/2-1,kfc4/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc2/2-1,kfc3/2);
    }    
    //p_colours[2+m_rev][m_anti]   = p_colours[0][m_anti]   = 500;
    //p_colours[3-m_rev][1-m_anti] = p_colours[1][1-m_anti] = 501;
  }
  else if (fl[0].IsDowntype() && fl[1].IsDowntype()) {
    m_schannel = false;
    if (!(fl[0].IsAnti()^fl[2].IsAnti())) {  
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc2/2);
    }
    if (!(fl[0].IsAnti()^fl[3].IsAnti())) {
      m_rev  = true;
      m_ckm1 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc4/2-1,kfc1/2);
      m_ckm2 = MODEL::s_model->ComplexMatrixElement(string("CKM"),kfc3/2-1,kfc2/2);
    }    
    //p_colours[2+m_rev][m_anti]   = p_colours[0][m_anti]   = 500;
    //p_colours[3-m_rev][1-m_anti] = p_colours[1][1-m_anti] = 501;
  }
}

double XS_f1f2b_f3f4b::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if (m_rev) { double help = u; u = t; t = help; }
  if (m_schannel) return sqr(m_pref_W)/((sqr(s-m_mw2)+m_mw2*m_ww2)) * (u*u);
  return sqr(m_pref_W)*CouplingFactor(0,2)/((sqr(t-m_mw2)+m_mw2*m_ww2)) * (u*u);
}

bool XS_f1f2b_f3f4b::SetColours(const Vec4D_Vector& mom) 
{
  if (m_schannel) {
    //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(s);
  }
  else {
    if (m_rev) {
      //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(u);
    }
    else {
      //m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = dabs(t);
    }
  }
  return true;
}


bool XS_f1f2b_f3f4b::SetColours() 
{
  return true;
}

double XS_f1f2b_f3f4b::KFactor(double scale) 
{ 
  return 1.; 
}
