#include "XS_4F.H"
#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Flow.H"
#include "Random.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;



template <> 
Single_XS *Single_XS::GetProcess<XS_f1f1_f1f1>(const size_t nin,const size_t nout,
					       const ATOOLS::Flavour *flavours,
					       const size_t nqed, const size_t nqcd)
{
  if (flavours[0]!=flavours[1] || 
      flavours[0]!=flavours[2] || flavours[0]!=flavours[3])  return NULL;
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()))         return NULL;
  if ((flavours[0].Charge()!=0. && ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn())                 return new XS_f1f1_f1f1(nin,nout,flavours); 
  return NULL;
}

XS_f1f1_f1f1::XS_f1f1_f1f1(const size_t nin,const size_t nout,
			   const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours),
  m_Z_on(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()), m_P_on(ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()),
  m_anti(int(flavours[0].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width())),
  m_sin2tw(ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq(flavours[0].Charge()),
  m_y3f((2.*int(flavours[0].IsUptype())-1)/2.),
  m_v(m_y3f-2.*m_eq*m_sin2tw), m_a(m_y3f),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw))
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  //   std::cout<<"Init f1f1 -> f1f1 : anti = "<<m_anti<<" : Z_on = "<<m_Z_on<<", photon_on = "<<m_P_on<<std::endl
  // 	   <<"(pref_Z = "<<m_pref_Z<<", pref_QED = "<<m_pref_qed<<" -> "<<sqrt(m_pref_qed)
  // 	   <<", aqed("<<ATOOLS::rpa.gen.Ecms()<<" ) = "<<m_aqed<<", sin2tw, cos2tw = "
  // 	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl;
}

double XS_f1f1_f1f1::operator()(double s,double t,double u) 
{
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
  return M_t + M_u + M_mix;
}

bool XS_f1f1_f1f1::SetColours(double s, double t, double u) 
{
  bool swap  = m_swaped;
  RestoreInOrder();
  (*this)(s,t,u);
  if (SetColours()) m_scale[PHASIC::stp::as] = t;
               else m_scale[PHASIC::stp::as] = u;
  if (swap) SwapInOrder();
  return true;
}


bool XS_f1f1_f1f1::SetColours() 
{
  if (M_t > (M_t+M_u) * ran.Get()) {
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

template <> 
Single_XS *Single_XS::GetProcess<XS_f1f1b_f1f1b>(const size_t nin,const size_t nout,
						 const ATOOLS::Flavour *flavours,
						 const size_t nqed, const size_t nqcd)
{
  //std::cout<<"Test this : "<<flavours[0]<<" "<<flavours[1]<<" "<<flavours[2]<<" "<<flavours[3]<<std::endl;
  if (flavours[0]!=flavours[1].Bar() || 
      flavours[0]!=flavours[2] || flavours[1]!=flavours[3])  return NULL;
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()))         return NULL;
  if ((flavours[0].Charge()!=0. && ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn())                 return new XS_f1f1b_f1f1b(nin,nout,flavours); 
  return NULL;
}

XS_f1f1b_f1f1b::XS_f1f1b_f1f1b(const size_t nin,const size_t nout,
			   const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours), 
  m_Z_on(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()), m_P_on(ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()),
  m_anti1(int(flavours[0].IsAnti())),m_anti2(int(flavours[2].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width())),
  m_sin2tw(ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq(flavours[0].Charge()),
  m_y3f((2.*int(flavours[0].IsUptype())-1)/2.),
  m_v(m_y3f-2.*m_eq*m_sin2tw), m_a(m_y3f),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw))
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  // std::cout<<"Init f1f1b -> f1f1b : anti = "<<m_anti1<<m_anti2
  // 	   <<" : Z_on = "<<m_Z_on<<", photon_on = "<<m_P_on<<std::endl
  // 	   <<"(pref_Z = "<<m_pref_Z<<", pref_QED = "<<m_pref_qed<<" -> "<<sqrt(m_pref_qed)
  // 	   <<", aqed("<<ATOOLS::rpa.gen.Ecms()<<" ) = "<<m_aqed<<", sin2tw, cos2tw = "
  // 	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl;
}

double XS_f1f1b_f1f1b::operator()(double s,double t,double u) 
{
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
  return 2.*(M_t + M_s + M_mix);
}

bool XS_f1f1b_f1f1b::SetColours(double s, double t, double u) 
{
  bool swap  = m_swaped;
  RestoreInOrder();
  (*this)(s,t,u);
  if (SetColours()) {
    m_scale[PHASIC::stp::as] = t;
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) m_scale[PHASIC::stp::as] = u;
  }
  else m_scale[PHASIC::stp::as] = s;
  if (swap) SwapInOrder();
  return true;
}


bool XS_f1f1b_f1f1b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran.Get()) {
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

template <> 
Single_XS *Single_XS::GetProcess<XS_f1f1b_f2f2b>(const size_t nin,const size_t nout,
						 const ATOOLS::Flavour *flavours,
						 const size_t nqed, const size_t nqcd)
{
  //std::cout<<"Test this : "<<flavours[0]<<" "<<flavours[1]<<" "<<flavours[2]<<" "<<flavours[3]<<std::endl;
  int  kfc1  = abs(flavours[0].Kfcode()), kfc2  = abs(flavours[2].Kfcode());
  if (flavours[0]!=flavours[1].Bar() || 
      flavours[0]!=flavours[2] || flavours[0]==flavours[3] || 
      flavours[2]!=flavours[3].Bar())                        return NULL;
  if (!flavours[0].IsFermion() || !flavours[2].IsFermion())  return NULL;
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if ((flavours[0].Charge()!=0. && flavours[2].Charge()!=0. &&
       ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      (ATOOLS::Flavour(ATOOLS::kf::W).IsOn() && 
       flavours[0].IsQuark() && flavours[2].IsQuark() && 
       ((kfc1/2==kfc1/2. && kfc2/2!=kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1))>0) ||
	(kfc1/2!=kfc1/2. && kfc2/2==kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1))>0))) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() )                return new XS_f1f1b_f2f2b(nin,nout,flavours); 
  //   std::cout<<"Checks "<<ATOOLS::Flavour(ATOOLS::kf::W).IsOn()<<" "
  // 	   <<(flavours[0].IsQuark() && flavours[1].IsQuark())<<" "
  // 	   <<kfc1<<" "<<kfc2<<" "
  // 	   <<rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1)<<std::endl;
  return NULL;
}

XS_f1f1b_f2f2b::XS_f1f1b_f2f2b(const size_t nin,const size_t nout,
			   const ATOOLS::Flavour *flavours) :
  Single_XS(nin,nout,flavours), 
  m_Z_on(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()), m_P_on(ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()),
  m_W_on(ATOOLS::Flavour(ATOOLS::kf::W).IsOn()),
  m_anti1(int(flavours[0].IsAnti())),m_anti2(int(flavours[2].IsAnti())),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width())),
  m_sin2tw(ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(flavours[0].Charge()),
  m_eq2(flavours[2].Charge()),
  m_y3f1((2.*int(flavours[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(flavours[2].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    int kfc1 = abs(flavours[0].Kfcode()), kfc2 = abs(flavours[2].Kfcode());
    if (flavours[0].IsQuark() && flavours[2].IsQuark()) { 
      if (kfc1/2==kfc1/2. && kfc2/2!=kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1);
      if (kfc1/2!=kfc1/2. && kfc2/2==kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
  // std::cout<<"Init f1f1b -> f2f2b : anti = "<<m_anti1<<m_anti2
  // 	   <<" : Z_on = "<<m_Z_on<<", photon_on = "<<m_P_on<<std::endl
  // 	   <<"(pref_Z = "<<m_pref_Z<<", pref_QED = "<<m_pref_qed<<" -> "<<sqrt(m_pref_qed)
  // 	   <<", aqed("<<ATOOLS::rpa.gen.Ecms()<<" ) = "<<m_aqed<<", sin2tw, cos2tw = "
  // 	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl;
}

double XS_f1f1b_f2f2b::operator()(double s,double t,double u) 
{
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
  return 2.*(M_s+M_t+M_mix);
}

bool XS_f1f1b_f2f2b::SetColours(double s, double t, double u) 
{
  bool swap  = m_swaped;
  RestoreInOrder();
  (*this)(s,t,u);
  if (SetColours()) {
    m_scale[PHASIC::stp::as] = t;
    if ((m_anti1&&(!m_anti2)) || ((!m_anti1)&&m_anti2)) m_scale[PHASIC::stp::as] = u;
  }
  else m_scale[PHASIC::stp::as] = s;
  if (swap) SwapInOrder();
  return true;
}


bool XS_f1f1b_f2f2b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran.Get()) {
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

template <> 
Single_XS *Single_XS::GetProcess<XS_f1f2_f1f2>(const size_t nin,const size_t nout,
					       const ATOOLS::Flavour *flavours,
					       const size_t nqed, const size_t nqcd)
{
  //std::cout<<"Test this : "<<flavours[0]<<" "<<flavours[1]<<" "<<flavours[2]<<" "<<flavours[3]<<std::endl;
  int kfc1 = flavours[0].Kfcode(), kfc2 = flavours[1].Kfcode();
  if (!(flavours[0]==flavours[2] && flavours[1]==flavours[3]) || 
      !(flavours[0]==flavours[3] && flavours[1]==flavours[2]) || 
      (flavours[0].IsAnti() && !flavours[1].IsAnti()) ||
      (!flavours[0].IsAnti() && flavours[1].IsAnti()))       return NULL;
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::photon).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::W).IsOn()))              return NULL;
  if ((flavours[0].Charge()!=0. && flavours[1].Charge()!=0. &&
       ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      (ATOOLS::Flavour(ATOOLS::kf::W).IsOn() && 
       flavours[0].IsQuark() && flavours[1].IsQuark() && 
       ((kfc1/2==kfc1/2. && kfc2/2!=kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1))>0) ||
	(kfc1/2!=kfc1/2. && kfc2/2==kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1))>0))) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() )                return new XS_f1f2_f1f2(nin,nout,flavours); 
  //std::cout<<"Checks "<<ATOOLS::Flavour(ATOOLS::kf::W).IsOn()<<" "
  // 	   <<(flavours[0].IsQuark() && flavours[1].IsQuark())<<" "
  // 	   <<kfc1<<" "<<kfc2<<" "
  // 	   <<rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1)<<std::endl;
  return NULL;
}

XS_f1f2_f1f2::XS_f1f2_f1f2(const size_t nin,const size_t nout,
			   const ATOOLS::Flavour *flavours) :
  Single_XS(nin,nout,flavours), 
  m_Z_on(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()), m_P_on(ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()),
  m_W_on(ATOOLS::Flavour(ATOOLS::kf::W).IsOn()),
  m_anti(int(flavours[0].IsAnti())),m_rev(int(flavours[0]==flavours[2])),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width())),
  m_sin2tw(ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(m_anti * flavours[0].Charge()),
  m_eq2(m_anti * flavours[1].Charge()),
  m_y3f1((2.*int(flavours[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(flavours[1].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    int kfc1 = flavours[0].Kfcode(), kfc2 = flavours[1].Kfcode();
    if (flavours[0].IsQuark() && flavours[1].IsQuark()) { 
      if (kfc1/2==kfc1/2. && kfc2/2!=kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1);
      if (kfc1/2!=kfc1/2. && kfc2/2==kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
  // std::cout<<"Init f1f2 -> f1f2 : anti = "<<m_anti
  // 	   <<" : Z_on = "<<m_Z_on<<", photon_on = "<<m_P_on<<" : W_on = "<<m_W_on<<std::endl
  // 	   <<"(pref_Z = "<<m_pref_Z<<", pref_QED = "<<m_pref_qed<<" -> "<<sqrt(m_pref_qed)
  // 	   <<", aqed("<<ATOOLS::rpa.gen.Ecms()<<" ) = "<<m_aqed<<", sin2tw, cos2tw = "
  // 	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl
  // 	   <<" m_eq1,2 = "<<m_eq1<<", "<<m_eq2
  // 	   <<" ("<<flavours[0].Charge()<<", "<<flavours[1].Charge()<<")"
  // 	   <<", CKM = "<<m_ckm<<std::endl;
}

double XS_f1f2_f1f2::operator()(double s,double t,double u) 
{
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
  return 2.*(M_t+M_u+M_mix);
}

bool XS_f1f2_f1f2::SetColours(double s, double t, double u) 
{
  bool swap  = m_swaped;
  RestoreInOrder();
  (*this)(s,t,u);
  if (SetColours()) {
    m_scale[PHASIC::stp::as] = t;
    if (m_rev) m_scale[PHASIC::stp::as] = u;
  }
  else {
    m_scale[PHASIC::stp::as] = u;
    if (m_rev) m_scale[PHASIC::stp::as] = t;
  }
  if (swap) SwapInOrder();
  return true;
}


bool XS_f1f2_f1f2::SetColours() 
{
  if (M_t > (M_t+M_u) * ran.Get()) {
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

template <> 
Single_XS *Single_XS::GetProcess<XS_f1f2b_f1f2b>(const size_t nin,const size_t nout,
						 const ATOOLS::Flavour *flavours,
						 const size_t nqed, const size_t nqcd)
{
  int kfc1 = flavours[0].Kfcode(), kfc2 = flavours[1].Kfcode();
  //std::cout<<"Test this : "<<flavours[0]<<" "<<flavours[1]<<" "<<flavours[2]<<" "<<flavours[3]<<std::endl;
  if (!(flavours[0]==flavours[2] && flavours[1]==flavours[3]) || 
      !(flavours[0]==flavours[3] && flavours[1]==flavours[2]) || 
      (flavours[0].IsAnti() && flavours[1].IsAnti()) ||
      (!flavours[0].IsAnti() && !flavours[1].IsAnti()))      return NULL;
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::photon).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::W).IsOn()))              return NULL;
  if ((flavours[0].Charge()!=0. && flavours[1].Charge()!=0. &&
       ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      (ATOOLS::Flavour(ATOOLS::kf::W).IsOn() && 
       flavours[0].IsQuark() && flavours[1].IsQuark() && 
       ((kfc1/2==kfc1/2. && kfc2/2!=kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1))>0) ||
	(kfc1/2!=kfc1/2. && kfc2/2==kfc2/2. && 
	 abs(rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1))>0))) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() )                return new XS_f1f2b_f1f2b(nin,nout,flavours); 
  //   std::cout<<"Checks "<<ATOOLS::Flavour(ATOOLS::kf::W).IsOn()<<" "
  // 	   <<(flavours[0].IsQuark() && flavours[1].IsQuark())<<" "
  // 	   <<kfc1<<" "<<kfc2<<" "
  // 	   <<rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1)<<std::endl;
  return NULL;
}

XS_f1f2b_f1f2b::XS_f1f2b_f1f2b(const size_t nin,const size_t nout,
			       const ATOOLS::Flavour *flavours) :
  Single_XS(nin,nout,flavours), 
  m_Z_on(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()), m_P_on(ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()),
  m_W_on(ATOOLS::Flavour(ATOOLS::kf::W).IsOn()),
  m_anti(int(flavours[0].IsAnti())),m_rev(int(flavours[0]==flavours[2])),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width())),
  m_mw2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass())),
  m_ww2(ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width())),
  m_sin2tw(ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(flavours[0].Charge()),
  m_eq2(flavours[1].Charge()),
  m_y3f1((2.*int(flavours[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(flavours[1].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_ckm(Complex(0.,0.)),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_pref_W((4.*M_PI*m_aqed)/(4.*m_sin2tw))
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (m_W_on) {
    int kfc1 = flavours[0].Kfcode(), kfc2 = flavours[1].Kfcode();
    if (flavours[0].IsQuark() && flavours[1].IsQuark()) { 
      if (kfc1/2==kfc1/2. && kfc2/2!=kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc2/2,kfc1/2-1);
      if (kfc1/2!=kfc1/2. && kfc2/2==kfc2/2.)
	m_ckm = rpa.gen.ComplexMatrixElement(string("CKM"),kfc1/2,kfc2/2-1);
    }
    if (abs(m_ckm)==0.) m_W_on = false;
  }
  //   std::cout<<"Init f1f2 -> f1f2 : anti = "<<m_anti
  // 	   <<" : Z_on = "<<m_Z_on<<", photon_on = "<<m_P_on<<" : W_on = "<<m_W_on<<std::endl
  // 	   <<"(pref_Z = "<<m_pref_Z<<", pref_QED = "<<m_pref_qed<<" -> "<<sqrt(m_pref_qed)
  // 	   <<", aqed("<<ATOOLS::rpa.gen.Ecms()<<" ) = "<<m_aqed<<", sin2tw, cos2tw = "
  // 	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl
  // 	   <<" m_eq1,2 = "<<m_eq1<<", "<<m_eq2
  // 	   <<" ("<<flavours[0].Charge()<<", "<<flavours[1].Charge()<<")"
  // 	   <<", CKM = "<<m_ckm<<std::endl;
}

double XS_f1f2b_f1f2b::operator()(double s,double t,double u) 
{
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
  return 2.*(M_t+M_s+M_mix);
}

bool XS_f1f2b_f1f2b::SetColours(double s, double t, double u) 
{
  bool swap  = m_swaped;
  RestoreInOrder();
  (*this)(s,t,u);
  if (SetColours()) {
    m_scale[PHASIC::stp::as] = t;
    if (m_rev) m_scale[PHASIC::stp::as] = u;
  }
  else {
    m_scale[PHASIC::stp::as] = s;
  }
  if (swap) SwapInOrder();
  return true;
}


bool XS_f1f2b_f1f2b::SetColours() 
{
  if (M_t > (M_t+M_s) * ran.Get()) {
    if (m_rev) {
      p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
      p_colours[2][1-m_anti] = p_colours[1][1-m_anti] = Flow::Counter();
    }
    else {
      p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
      p_colours[3][1-m_anti] = p_colours[1][1-m_anti] = Flow::Counter();
    }
    return true;
  }
  p_colours[1][1-m_anti] = p_colours[0][m_anti] = Flow::Counter();
  if (!m_rev) p_colours[3][1-m_anti] = p_colours[2][m_anti] = Flow::Counter();
         else p_colours[2][1-m_anti] = p_colours[3][m_anti] = Flow::Counter();
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


template <> 
Single_XS *Single_XS::GetProcess<XS_f1f2b_f3f4b>(const size_t nin,const size_t nout,
						 const ATOOLS::Flavour *flavours,
						 const size_t nqed, const size_t nqcd)
{
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!flavours[0].IsQuark() || !flavours[1].IsQuark())      return NULL;
  if (!flavours[2].IsQuark() || !flavours[3].IsQuark())      return NULL;
  bool up[4], anti[4];
  for (short int i=0;i<4;++i) {
    if (flavours[i].IsUptype()) up[i] = true;
    else if (!flavours[i].IsDowntype())                      return NULL;
    else up[i] = false;
    anti[i]    = flavours[i].IsAnti();
  }
  if (anti[0]==anti[1] || anti[2]==anti[3])                  return NULL;
  if (!ATOOLS::Flavour(ATOOLS::kf::W).IsOn())                return NULL;
  if ((up[0] && !up[1] && up[2] && !up[3] && anti[0]==anti[2]) ||
      (!up[0] && up[1] && up[2] && !up[3] && anti[1]==anti[2]) ||
      (up[0] && !up[1] && !up[2] && up[3] && anti[1]==anti[2]) ||
      (!up[0] && up[1] && !up[2] && up[3] && anti[0]==anti[2])){ 
    // Check if CKM allowed.
    /*
    double m_ckm2_0 = 1., m_ckm2_1 = 1.;
    int ints[4];
    for (short int i=0;i<4;++i) ints[i] = ATOOLS::kf_table.ToInt(flavours[i].Kfcode());
    if (flavours[0].IsQuark())
      m_ckm2_0 = std::abs(sqr(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[0]/2-1,ints[1]/2)));
    if (flavours[2].IsQuark())
      m_ckm2_1 = std::abs(sqr(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[2]/2-1,ints[3]/2)));
    if (m_ckm_0<1.e-6 || m_ckm_1<1.e-6)                      return NULL;
    */
    return new XS_f1f2b_f3f4b(nin,nout,flavours); 
  }
  return NULL;
}

XS_f1f2b_f3f4b::XS_f1f2b_f3f4b(const size_t nin,const size_t nout,
			       const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  int ints[4];
  for (short int i=0;i<4;++i) ints[i] = ATOOLS::kf_table.ToInt(flavours[i].Kfcode());
  m_barred[0] = flavours[0].IsAnti();
  m_barred[1] = flavours[2].IsAnti();
  if (flavours[0].IsDowntype()) std::swap(ints[0],ints[1]);
  if (flavours[2].IsDowntype()) std::swap(ints[2],ints[3]);
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  p_colours[0][flavours[0].IsAnti()]=p_colours[1][1-flavours[0].IsAnti()]=ATOOLS::Flow::Counter();
  p_colours[2][flavours[2].IsAnti()]=p_colours[3][1-flavours[2].IsAnti()]=ATOOLS::Flow::Counter();

  /*
  m_ckm2[0]   = m_ckm2[1]=1.;
  m_W_on      = true;
  if (flavours[0].IsQuark())
    m_ckm2[0] = std::abs(sqr(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[0]/2-1,ints[1]/2)));
  if (flavours[2].IsQuark())
    m_ckm2[1] = std::abs(sqr(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[2]/2-1,ints[3]/2)));
  if (m_ckm[0]<1.e-6 || m_ckm[1]<1.e-6) m_W_on = false;
  */
  m_mw2       = ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass());
  m_ww2       = ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width());
  m_aqed      = MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  m_sin2tw    = ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"));
  m_resonances.push_back(ATOOLS::Flavour(ATOOLS::kf::W));
}

double XS_f1f2b_f3f4b::operator()(double s,double t,double u) 
{
  if (m_W_on) {
    double abst_2 = p_momenta[0]*p_momenta[2];
    if (m_swaped) abst_2 = p_momenta[1]*p_momenta[2];
    return ATOOLS::sqr(M_PI*m_aqed/m_sin2tw)*16/3.*m_ckm2[0]*m_ckm2[1]*
      ATOOLS::sqr(abst_2)/(ATOOLS::sqr(s-m_mw2)+m_mw2*m_ww2); 
  }
  return 0.;
}

bool XS_f1f2b_f3f4b::SetColours(double s,double t,double u) 
{ 
  m_scale[PHASIC::stp::fac]=s;
  p_colours[0][m_barred[0]]=ATOOLS::Flow::Counter();
  p_colours[1][1-m_barred[0]]=p_colours[0][m_barred[0]];
  p_colours[2][m_barred[1]]=ATOOLS::Flow::Counter();
  p_colours[3][1-m_barred[1]]=p_colours[2][m_barred[1]];
  return true; 
}

double XS_f1f2b_f3f4b::KFactor(const double scale) 
{
  return 1.;
}

