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
  if (ATOOLS::rpa.gen.Model()==ATOOLS::Model_Type::pure_QCD) return NULL;
  if (nqcd!=0 || nqed!=2)                                    return NULL;
  if (!(ATOOLS::Flavour(ATOOLS::kf::Z).IsOn() ||
	ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()))         return NULL;
  for (int i=1;i<4;i++) { if (flavours[0]!=flavours[i])      return NULL; }
  if ((flavours[0].Charge()!=0. && ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) ||
      ATOOLS::Flavour(ATOOLS::kf::Z).IsOn())                 return new XS_f1f1_f1f1(nin,nout,flavours); 
  return NULL;
}

XS_f1f1_f1f1::XS_f1f1_f1f1(const size_t nin,const size_t nout,
			   const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours), m_Z_on(true), m_P_on(true), 
  m_anti(-2*int(flavours[0].IsAnti())+1),
  m_y3f(2*int(flavours[0].IsUptype())-1)
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_aqed      = MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  m_eq        = m_anti * flavours[0].Charge();
  m_mz2       = ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass());
  m_wz2       = ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width());
  m_sin2tw    = ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"));
  m_cos2tw    = 1.-m_sin2tw;
  m_pref_qed  = (4.*M_PI*m_aqed*m_eq*m_eq);
  m_pref_Z    = (4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw);
  if (!ATOOLS::Flavour(ATOOLS::kf::Z).IsOn())      m_Z_on = false;
  if (!ATOOLS::Flavour(ATOOLS::kf::photon).IsOn()) m_P_on = false;
  std::cout<<"Init f1f1 -> f1f1 : "<<m_anti<<" : "<<m_Z_on<<", "<<m_P_on
	   <<"("<<m_pref_Z<<", "<<m_pref_qed<<" <- "<<sqrt(4.*M_PI*m_aqed)*m_eq<<","
	   <<m_sin2tw<<","<<m_cos2tw<<")"<<std::endl;
}

double XS_f1f1_f1f1::operator()(double s,double t,double u) 
{
  M_t = 0., M_u = 0., M_mix = 0.;
  if (m_P_on) {
    M_t   +=     sqr(m_pref_qed)    * (s*s+u*u)/(t*t);
    M_mix +=  2.*sqr(m_pref_qed)/3. * (s*s)/(t*u);
    M_u   +=     sqr(m_pref_qed)    * (s*s+t*t)/(u*u); 
  }
  if (m_Z_on) {
    M_t   +=     sqr(m_pref_Z) /(4.*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      (sqr(sqr(m_y3f-2.*m_eq*m_sin2tw)+sqr(2.*m_eq*m_sin2tw)) * (s*s+u*u) +
       sqr(sqr(m_y3f-2.*m_eq*m_sin2tw)-sqr(2.*m_eq*m_sin2tw)) * (s*s-u*u));
    M_mix +=  4.*sqr(m_pref_Z)/3. /4. *
      ((t-m_mz2)*(u-m_mz2)-m_mz2*m_wz2)/(sqr((t-m_mz2)*(u-m_mz2)-m_mz2*m_wz2)+m_mz2*m_wz2*(t+u-2.*m_mz2)) *
      (sqr(sqr(m_y3f-2.*m_eq*m_sin2tw))+sqr(sqr(2.*m_eq*m_sin2tw))) * s*s;
    M_u   +=     sqr(m_pref_Z) /(4.*(sqr(u-m_mz2)+m_mz2*m_wz2)) *
      (sqr(sqr(m_y3f-2.*m_eq*m_sin2tw)+sqr(2.*m_eq*m_sin2tw)) * (s*s+t*t) +
       sqr(sqr(m_y3f-2.*m_eq*m_sin2tw)-sqr(2.*m_eq*m_sin2tw)) * (s*s-t*t));
  }
  if (m_P_on && m_Z_on) {
    M_t   +=  m_pref_qed*m_pref_Z /2. *                   
      (t-m_mz2)/(t*(sqr(t-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_y3f-4.*m_eq*m_sin2tw) * (s*s+u*u) + 1. * (s*s-u*u));
    M_mix += 2.*m_pref_qed*m_pref_Z/3./2. * 
      ((u-m_mz2)/(t*(sqr(u-m_mz2)+m_mz2*m_wz2)) + 
       (t-m_mz2)/(u*(sqr(t-m_mz2)+m_mz2*m_wz2)) )  * 
      ( sqr(m_y3f-2.*m_eq*m_sin2tw) + sqr(2.*m_eq*m_sin2tw) ) * (s*s);
    M_u   +=  m_pref_qed*m_pref_Z /2.  * 
      (u-m_mz2)/(u*(sqr(u-m_mz2)+m_mz2*m_wz2)) *
      ( sqr(m_y3f-4.*m_eq*m_sin2tw) * (s*s+t*t) + 1. * (s*s-t*t));
  }
  return M_t + M_u + M_mix;
}

bool XS_f1f1_f1f1::SetColours(double s, double t, double u) 
{
  bool swap = m_swaped;
  RestoreInOrder();

  m_scale[PHASIC::stp::as] = (2.*s*t*u)/(s*s+t*t+u*u);
  
  M_t = 1. - 2.*(u*s) / (t*t);
  M_u = 1. - 2.*(s*t) / (u*u);
  
  bool result = SetColours();
  if (swap) SwapInOrder();
  return result;
}


bool XS_f1f1_f1f1::SetColours() 
{
  if (M_t > (M_t+M_u) * ran.Get()) {
    p_colours[2][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[3][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  else {
    p_colours[3][m_anti] = p_colours[0][m_anti] = Flow::Counter();
    p_colours[2][m_anti] = p_colours[1][m_anti] = Flow::Counter();
  }
  return 1;
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

double XS_f1f2b_f3f4b::KFactor(double scale) 
{
  return 1.;
}
