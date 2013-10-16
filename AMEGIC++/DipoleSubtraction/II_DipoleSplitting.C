#include "AMEGIC++/DipoleSubtraction/II_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

void II_DipoleSplitting::SetMomenta(const Vec4D *mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = (m_pk*m_pi-m_pi*m_pj-m_pj*m_pk)/(m_pk*m_pi);

  m_ptk  = m_pk;
  m_ptij = m_xijk*m_pi;

  Vec4D K  = m_pi-m_pj+m_pk;
  Vec4D Kt = m_ptij+m_pk;
  Vec4D KKt = K+Kt;
  for(int i=2;i<=m_m;i++) m_mom[i]-=2.*(m_mom[i]*KKt/KKt.Abs2()*KKt-m_mom[i]*K/K.Abs2()*Kt);

  m_vi   = (m_pi*m_pj)/(m_pi*m_pk);
  m_a = m_vi;

  m_Q2 = (-m_pi+m_pj-m_pk).Abs2();
  if (m_es==0) {
    m_kt2 = m_Q2*(1.-m_xijk)/m_xijk*m_vi;
  }
  else {
  m_kt2 = m_Q2/m_xijk*m_vi;
  switch (m_ft) {
  case 1:
    m_kt2*=(1.-m_xijk);
    break;
  case 4:
    m_kt2*=(1.-m_xijk);
    break;
  }
  }

//   m_pt1  =    m_pj;
//   m_pt2  =-1.*m_vi*m_pk;
  m_pt1  =    m_pj-m_vi*m_pk;
  m_pt2  =    m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_xijk)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 1.-2.*m_xijk*(1.-m_xijk);
    m_av  = m_sff;
    break;
  case 3:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk;
    break;
  case 4:
    m_sff = m_xijk/(1.-m_xijk)+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk;
  }
}

double II_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;  
  switch (m_ft) {
  case 1:
    h*=m_sff;
    return h;
  case 2:
    h*=m_sff;   
    return h;
  case 3:
    return h*m_sff*CSC.TR/CSC.CA;
  case 4:
    h*=2.*m_sff;
    return h;
  }
  return 0.;
}

void II_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xijk/(1.-m_xijk)/4.);
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xijk/(1.-m_xijk)/2.);
    break;
  }
}
