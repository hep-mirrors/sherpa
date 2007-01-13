#include "BFKL_Sudakov.H"

#include "LO_LEV.H"
#include "Random.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace ATOOLS;

BFKL_Sudakov::BFKL_Sudakov():
  p_levs(new LEV_Group()),
  m_weight(0.0) 
{
  p_levs->Add(new G_GG());
}

BFKL_Sudakov::~BFKL_Sudakov()
{
  delete p_levs;
}

bool BFKL_Sudakov::Initialize()
{
  m_asmax=(*MODEL::as)(m_kt2min);
  p_levs->SetKT2Min(m_kt2min);
  m_y=m_ya;
  return true;
}

bool BFKL_Sudakov::Dice()
{
  if (m_kt2max<=m_kt2min) return false;
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // kt limit
  p_levs->SetKT2Max(m_kt2max);
  // get last integral
  m_integral=p_levs->MajorIntegral();
  // set branching probability
  double qt2(m_kt2min);
  switch (m_ktscheme) {
  case 0:
    break;
  case 1: {
    qt2=m_q.PPerp2();
    break;
  }
  default:
    THROW(fatal_error,"Invalid scale scheme.");
  }
  m_gamma=(*MODEL::as)(qt2)/(2.0*M_PI)*
    log(m_q.PPerp2()/m_kt2min)*m_integral;
  // dice new y
  if (m_ya>m_yb) m_y+=log(rn[0])/m_gamma;
  else m_y-=log(rn[0])/m_gamma;
  // dice new kt2
  m_kt2=m_kt2min*pow(m_kt2max/m_kt2min,rn[1]);
  // dice new phi
  m_phi=2.0*M_PI*rn[2];
  // select splitting
  p_levs->SelectLEV(rn[3]);
  return m_ya>m_yb?m_y>m_yb:m_y<m_yb;
}

bool BFKL_Sudakov::Approve(const ATOOLS::Vec4D &k1,
			   const ATOOLS::Vec4D &q1,
			   const ATOOLS::Vec4D &k2,
			   const ATOOLS::Vec4D &q2)
{
  // lower qt and kt limits
  if (q2.PPerp2()<m_kt2min || 
      k2.PPerp2()<m_kt2min) return false;
  // rejection weight
  double weight(1.0);
  // splitting weight
  weight*=p_levs->Selected()->Value(k1,q1,k2,q2)/
    p_levs->Selected()->MajorValue(k1,q1,k2,q2);
  // veto emission
  if (weight<ran.Get()) return false;
  // reweighting weight
  m_weight=1.0;
  // divide by \frac{d}{dy}\ln\Delta(y,\bar y)
  m_weight/=m_gamma;
  // set scale
  double kt2(m_kt2min);
  switch (m_ktscheme) {
  case 0:
    break;
  case 1: {
    kt2=k2.PPerp2();
    break;
  }
  default:
    THROW(fatal_error,"Invalid scale scheme.");
  }
  // coupling & splitting weight
  m_weight*=(*MODEL::as)(kt2)/(2.0*M_PI)*m_integral;
  // kt integration domain weight
  m_weight*=log(m_kt2max/m_kt2min);
  return true;
}

