#include "BFKL_Sudakov.H"

#include "LO_LEV.H"
#include "Random.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace ATOOLS;

BFKL_Sudakov::BFKL_Sudakov():
  p_levs(new LEV_Group()) 
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
  m_oldy=m_y=m_ya;
  return true;
}

bool BFKL_Sudakov::Dice()
{
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // kt limit
  p_levs->SetKT2Max(m_kt2max);
  // get last integral
  m_integral=p_levs->MajorIntegral();
  // set branching probability
  m_gamma=m_asmax/(2.0*M_PI)*log(m_kt2max/m_kt2min)*m_integral;
  // dice new y
  m_y+=log(rn[0])/m_gamma;
  // dice new kt2
  m_kt2=m_kt2min*pow(m_kt2max/m_kt2min,rn[1]);
  // dice new phi
  m_phi=2.0*M_PI*rn[2];
  // select splitting
  p_levs->SelectLEV(rn[3]);
  return m_y>m_yb;
}

bool BFKL_Sudakov::Approve(const Vector_Vector &k,
			   const ATOOLS::Vec4D &q1,
			   const ATOOLS::Vec4D &q2)
{
  double qt22(q2.PPerp2()), kt2(k.back().PPerp2());
  // lower qt and kt limits
  if (qt22<m_kt2min || kt2<m_kt2min) return false;
  double weight(1.0);
  // coupling weight
  weight*=(*MODEL::as)(kt2)/m_asmax;
  // splitting weight
  weight*=p_levs->Selected()->Value(m_y,m_kt2)/
    p_levs->Selected()->MajorValue(m_y,m_kt2);
  // propagator weight
  weight*=kt2/qt22;
  // veto emission
  if (weight<ran.Get()) return false;
  // store old y
  m_oldy=m_y;
  return true;
}

