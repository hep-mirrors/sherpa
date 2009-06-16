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

//   m_pt1  =    m_pj;
//   m_pt2  =-1.*m_vi*m_pk;
  m_pt1  =    m_pj-m_vi*m_pk;
  m_pt2  =    m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_xijk)-(1.+m_xijk);
    break;
  case 2:
    m_sff = 1.-2.*m_xijk*(1.-m_xijk);
    break;
  case 3:
    m_sff = m_xijk;
    break;
  case 4:
    m_sff = m_xijk/(1.-m_xijk)+m_xijk*(1.-m_xijk);
  }
  if (1) {
    double pf = (m_vi)/m_alpha;
    switch (m_ft) {
    case 1:
      m_sff -= pf*Vcijf(1);
      break;
    case 2:
      break;
    case 3:
      break;
    case 4:
      m_sff -= pf*(0.5*Vcijf(4)+m_nf*Vcijf(3));
      break;
    }
  }
}

double II_DipoleSplitting::GetF()
{
   if (m_vi>m_alpha) return 0.;
   if (m_vi<=m_amin) {
     return nan;
   }
  double h=m_spfac/(2.*m_pi*m_pj)/m_xijk;  
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
