#include "Off_Shell_EW.H"

#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Running_AlphaS.H"
#include "Flow.H"

using namespace EXTRAXS;
using namespace MODEL;

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_qqb_llb>(const size_t nin,const size_t nout,
						    const ATOOLS::Flavour *flavours)
{
  if ((flavours[2].IsLepton() && flavours[3]==flavours[2].Bar() && 
       flavours[0].IsQuark() && flavours[1]==flavours[0].Bar()) ||
      (flavours[0].IsLepton() && flavours[1]==flavours[0].Bar() && 
       flavours[2].IsQuark() && flavours[3]==flavours[2].Bar())){ 
    return new Off_Shell_qqb_llb(nin,nout,flavours); 
  }
  return NULL;
}

Off_Shell_qqb_llb::Off_Shell_qqb_llb(const size_t nin,const size_t nout,
				     const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  MZ2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass());
  GZ2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width());
  alpha  = aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  sin2tw = ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"));
  if (ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()) kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  else kappa  = 0.;
  mass     = flavours[2].Mass();
  qe       = flavours[0].Charge();
  qf       = flavours[2].Charge();
  ae       = flavours[0].IsoWeak();      
  af       = flavours[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  colfac   = 1.;
  kswitch  = 0;  
  fac      = 2./(3.*M_PI);
  fin      = 2.*M_PI/9. - 7./(3.*M_PI) + 9./(3.*M_PI);
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (flavours[2].IsQuark()) {
    barred = flavours[2].IsAnti();
    p_colours[2][barred] = p_colours[3][1-barred] = 500;
    colfac = 3.;
  }
  if (flavours[0].IsQuark())  {
    barred = flavours[0].IsAnti();
    p_colours[0][barred] = p_colours[1][1-barred] = 500;
    colfac  = 1./3.;
    kswitch = 1;
  }
  m_resonances.push_back(ATOOLS::Flavour(ATOOLS::kf::Z));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
}

double Off_Shell_qqb_llb::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  chi1  = kappa * s * (s-MZ2)/(ATOOLS::sqr(s-MZ2) + GZ2*MZ2);
  chi2  = ATOOLS::sqr(kappa * s)/(ATOOLS::sqr(s-MZ2) + GZ2*MZ2);
  term1 = (1+ATOOLS::sqr(1.+2.*t/s)) * (ATOOLS::sqr(qf*qe) + 2.*(qf*qe*vf*ve) * chi1 +
				(ae*ae+ve*ve) * (af*af+vf*vf) * chi2);
  term2 = (1.+2.*t/s) * (4. * qe*qf*ae*af * chi1 + 8. * ae*ve*af*vf * chi2);
  return ATOOLS::sqr(4.*M_PI*alpha) * colfac * (term1+term2); 
}

bool Off_Shell_qqb_llb::SetColours(double s,double t,double u) 
{ 
  m_scale[PHASIC::stp::as]=s;
  return true; 
}

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_q1q2b_lnulb>(const size_t nin,const size_t nout,
							const ATOOLS::Flavour *flavours)
{
  if (!flavours[0].IsQuark() || !flavours[2].IsLepton()) return NULL;
  bool up[4], anti[4];
  for (short int i=0;i<4;++i) {
    if (flavours[i].IsUptype()) up[i]=true;
    else if (!flavours[i].IsDowntype()) return NULL;
    else up[i]=false;
    anti[i]=flavours[i].IsAnti();
  }
  if (anti[0]==anti[1] || anti[2]==anti[3]) return NULL;
  if ((up[0] && !up[1] && up[2] && !up[3] && anti[0]==anti[2]) ||
      (!up[0] && up[1] && up[2] && !up[3] && anti[1]==anti[2]) ||
      (up[0] && !up[1] && !up[2] && up[3] && anti[1]==anti[2]) ||
      (!up[0] && up[1] && !up[2] && up[3] && anti[0]==anti[2])){ 
    return new Off_Shell_q1q2b_q3q4b(nin,nout,flavours); 
  }
  return NULL;
}

Off_Shell_q1q2b_lnulb::Off_Shell_q1q2b_lnulb(const size_t nin,const size_t nout,
					     const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  int ints[4];
  for (short int i=0;i<4;++i) ints[i]=ATOOLS::kf_table.ToInt(flavours[i].Kfcode());
  if (flavours[0].IsDowntype()) std::swap(ints[0],ints[1]);
  if (flavours[2].IsDowntype()) std::swap(ints[2],ints[3]);
  m_ckm2[0]=std::abs(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[0]/2-1,ints[1]/2));
  m_ckm2[1]=1.;
  m_mw2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass());
  m_ww2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width());
  m_aqed=MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  m_sin2tw=ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  p_colours[0][flavours[0].IsAnti()]=p_colours[1][1-flavours[0].IsAnti()]=ATOOLS::Flow::Counter();
  m_resonances.push_back(ATOOLS::Flavour(ATOOLS::kf::W));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
}

double Off_Shell_q1q2b_lnulb::operator()(double s,double t,double u) 
{
  double sc=p_momenta[0]*p_momenta[2];
  if (m_swaped) sc=p_momenta[1]*p_momenta[2];
  return ATOOLS::sqr(M_PI*m_aqed/m_sin2tw)*16./3.*m_ckm2[0]*m_ckm2[1]*
    ATOOLS::sqr(sc)/(ATOOLS::sqr(s-m_mw2)+m_mw2*m_ww2); 
}

bool Off_Shell_q1q2b_lnulb::SetColours(double s,double t,double u) 
{ 
  m_scale[PHASIC::stp::as]=s;
  return true; 
}

template <> 
Single_XS *Single_XS::GetProcess<Off_Shell_q1q2b_q3q4b>(const size_t nin,const size_t nout,
							const ATOOLS::Flavour *flavours)
{
  if (!flavours[0].IsQuark() || !flavours[2].IsQuark()) return NULL;
  bool up[4], anti[4];
  for (short int i=0;i<4;++i) {
    if (flavours[i].IsUptype()) up[i]=true;
    else if (!flavours[i].IsDowntype()) return NULL;
    else up[i]=false;
    anti[i]=flavours[i].IsAnti();
  }
  if (anti[0]==anti[1] || anti[2]==anti[3]) return NULL;
  if ((up[0] && !up[1] && up[2] && !up[3] && anti[0]==anti[2]) ||
      (!up[0] && up[1] && up[2] && !up[3] && anti[1]==anti[2]) ||
      (up[0] && !up[1] && !up[2] && up[3] && anti[1]==anti[2]) ||
      (!up[0] && up[1] && !up[2] && up[3] && anti[0]==anti[2])){ 
    return new Off_Shell_q1q2b_q3q4b(nin,nout,flavours); 
  }
  return NULL;
}

Off_Shell_q1q2b_q3q4b::Off_Shell_q1q2b_q3q4b(const size_t nin,const size_t nout,
					     const ATOOLS::Flavour *flavours):
  Single_XS(nin,nout,flavours) 
{
  int ints[4];
  for (short int i=0;i<4;++i) ints[i]=ATOOLS::kf_table.ToInt(flavours[i].Kfcode());
  if (flavours[0].IsDowntype()) std::swap(ints[0],ints[1]);
  if (flavours[2].IsDowntype()) std::swap(ints[2],ints[3]);
  m_ckm2[0]=std::abs(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[0]/2-1,ints[1]/2));
  m_ckm2[1]=std::abs(ATOOLS::rpa.gen.ComplexMatrixElement("CKM",ints[2]/2-1,ints[3]/2));
  if (m_ckm2[1]==0.) m_ckm2[1]=1.;
  m_mw2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Mass());
  m_ww2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::W).Width());
  m_aqed=MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  m_sin2tw=ATOOLS::rpa.gen.ScalarConstant(std::string("sin2_thetaW"));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  p_colours[0][flavours[0].IsAnti()]=p_colours[1][1-flavours[0].IsAnti()]=ATOOLS::Flow::Counter();
  p_colours[2][flavours[2].IsAnti()]=p_colours[3][1-flavours[2].IsAnti()]=ATOOLS::Flow::Counter();
  m_resonances.push_back(ATOOLS::Flavour(ATOOLS::kf::W));
  m_nvector=m_nvector+2;
  CreateMomenta(m_nvector);
}

double Off_Shell_q1q2b_q3q4b::operator()(double s,double t,double u) 
{
  double sc=p_momenta[0]*p_momenta[2];
  if (m_swaped) sc=p_momenta[1]*p_momenta[2];
  return ATOOLS::sqr(M_PI*m_aqed/m_sin2tw)*16.*m_ckm2[0]*m_ckm2[1]*
    ATOOLS::sqr(sc)/(ATOOLS::sqr(s-m_mw2)+m_mw2*m_ww2); 
}

bool Off_Shell_q1q2b_q3q4b::SetColours(double s,double t,double u) 
{ 
  m_scale[PHASIC::stp::as]=s;
  return true; 
}

