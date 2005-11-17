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
  // set naive kt limit
  p_levs->SetKT2Max(m_kt2max/4.0);
  // get last integral
  m_integral=p_levs->MajorIntegral();
  // dice new y
  double kt2max(Min(m_kt2min*exp(dabs(m_oldy-m_yb)),m_kt2max));
  m_y+=log(rn[0])/(m_asmax/(2.0*M_PI)*log(m_kt2max/m_kt2min)*m_integral);
  // dice new kt2
  m_kt2=m_kt2min*pow(m_kt2max/m_kt2min,rn[1]);
  // dice new phi
  m_phi=2.0*M_PI*rn[2];
  // select splitting
  p_levs->SelectLEV(rn[3]);
  return m_y>m_yb;
}

bool BFKL_Sudakov::CalculateWeight(const Vector_Vector &k,
				   const ATOOLS::Vec4D &q)
{
  m_weight=1.0;
  double qt2(q.PPerp2());
  // qt limits
  if (qt2<m_kt2min) return false;
  // coupling weight
  if ((*MODEL::as)(m_kt2)<m_asmax*ran.Get()) return false;
  // kt limits
  double kt2max(Min(m_kt2min*exp(dabs(m_oldy-m_y)),m_kt2max));
  if (m_kt2>kt2max) return false;
//   // qt weight
//   if (log(qt2/m_kt2min)<log(m_kt2max/m_kt2min)*ran.Get()) return false;
  // splitting weight
  if (p_levs->Selected()->Value(m_y,m_kt2)<
      p_levs->Selected()->Value(m_y,m_kt2)*ran.Get()) return false;
  m_oldy=m_y;
  return true;
}
