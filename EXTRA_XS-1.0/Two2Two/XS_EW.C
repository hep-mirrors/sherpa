#include "XS_EW.H"
#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Flow.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

template <> 
Single_XS *Single_XS::GetProcess<XS_q1q2b_q3q4b>(const size_t nin,const size_t nout,
							const ATOOLS::Flavour *flavours)
{
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
    return new XS_q1q2b_q3q4b(nin,nout,flavours); 
  }
  return NULL;
}

XS_q1q2b_q3q4b::XS_q1q2b_q3q4b(const size_t nin,const size_t nout,
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
  p_colours[0][flavours[2].IsAnti()]=p_colours[1][1-flavours[2].IsAnti()]=ATOOLS::Flow::Counter();
  m_resonances.push_back(ATOOLS::Flavour(ATOOLS::kf::W));
}

double XS_q1q2b_q3q4b::operator()(double s,double t,double u) 
{
  double sc=p_momenta[0]*p_momenta[2];
  if (m_swaped) sc=p_momenta[1]*p_momenta[2];
  return ATOOLS::sqr(M_PI*m_aqed/m_sin2tw)*16*m_ckm2[0]*m_ckm2[1]*
    ATOOLS::sqr(sc)/(ATOOLS::sqr(s-m_mw2)+m_mw2*m_ww2); 
}

bool XS_q1q2b_q3q4b::SetColours(double s,double t,double u) 
{ 
  m_scale=s;
  return true; 
}

double XS_q1q2b_q3q4b::KFactor(double scale) 
{
  return 1.;
}
