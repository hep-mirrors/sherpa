#include "BFKL_Sudakov.H"

#include "LO_LEV.H"
#include "Random.H"
#include "Running_AlphaS.H"

#define N_C 3.0

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

double BFKL_Sudakov::AlphaSBar(const double &kt2) const
{
  return (*MODEL::as)(kt2)*N_C/M_PI;
}

bool BFKL_Sudakov::Initialize()
{
  m_asbarmax=AlphaSBar(m_kt2min);
  m_y=m_ya;
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
  double qt2max(m_qt2+m_kt2max+2.0*sqrt(m_qt2*m_kt2max));
  m_y+=log(rn[0])/(m_asbarmax*log(qt2max/m_kt2min)*m_integral);
  // dice new qt2
  m_qt2=m_kt2min*pow(qt2max/m_kt2min,rn[1]);
  // dice new phi
  m_phi=2.0*M_PI*rn[2];
  // select splitting
  p_levs->SelectLEV(rn[3]);
  return m_y>m_yb;
}

bool BFKL_Sudakov::CalculateWeight(const Vector_Vector &moms)
{
  m_weight=1.0;
  double kt2(moms.back().PPerp2());
  // kt limits
  if (kt2<m_kt2min) return false;
  // coupling weight
  if (AlphaSBar(m_qt2)<m_asbarmax*ran.Get()) return false;
  // angular ordering weight
  if (sqrt(m_kt2min*m_kt2max)*exp(-dabs(m_y))<kt2) return false;
  // splitting weight
  if (p_levs->Selected()->Value(m_y,kt2)<
      p_levs->Selected()->Value(m_y,kt2)*ran.Get()) return false;
  return true;
}
