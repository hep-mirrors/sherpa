#include "AMEGIC++/DipoleSubtraction/FI_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

void FI_DipoleSplitting::SetMomenta(const Vec4D* mom )
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = 1.-m_pi*m_pj/(m_pj*m_pk+m_pk*m_pi);

  m_ptk  = m_xijk*m_pk;
  m_ptij = m_pi+m_pj-(1.-m_xijk)*m_pk;

  m_zi   = (m_pi*m_ptk)/(m_ptij*m_ptk);
  m_zj   = 1.-m_zi;

  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_zi+(1.-m_xijk))-(1.+m_zi);
    break;
  case 2:
    m_sff = 2./(1.-m_zj+(1.-m_xijk))-(1.+m_zj);
    break;
  case 3:
    m_sff = 1.;
    break;
  case 4:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
  }
}

double FI_DipoleSplitting::GetF()
{
   if (1.-m_xijk>m_alpha) return 0.;
   if ((1.-m_xijk)<=m_amin) {
     return nan;
   }
  double h=SPFac()/(2.*m_pi*m_pj)/m_xijk;  
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

void FI_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,m_sff/(4.*m_zi*m_zj));
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,-m_sff/(2.*m_zi*m_zj));
    break;
  }
}


void FI_MassiveDipoleSplitting::SetMomenta(const Vec4D* mom )
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = 1.-(m_pi*m_pj-0.5*(m_mij-m_mi-m_mj))/(m_pj*m_pk+m_pk*m_pi);
  
  m_ptk  = m_xijk*m_pk;
  m_ptij = m_pi+m_pj-(1.-m_xijk)*m_pk;

  m_zi   = (m_pi*m_pk)/(m_pj*m_pk+m_pk*m_pi);
  m_zj   = 1.-m_zi;

//   m_pt1   =     m_zi*m_pi;
//   m_pt2   = -1.*m_zj*m_pj;
  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(2.-m_zi-m_xijk)-(1.+m_zi)-m_mij/(m_pi*m_pj);
    break;
  case 2:
    m_sff = 2./(2.-m_zj-m_xijk)-(1.+m_zj)-m_mij/(m_pi*m_pj);
    break;
  case 3:
    m_sff = 1.;
    break;
  case 4:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
  }
}

double FI_MassiveDipoleSplitting::GetF()
{
   if (1.-m_xijk>m_alpha) return 0.;
   if ((1.-m_xijk)<=m_amin) {
     return nan;
   }
  double h=SPFac()/((m_pi+m_pj).Abs2()-m_mij)/m_xijk;
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

void FI_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,-m_sff*(m_pi+m_pj).Abs2()/(4.*m_pt1.Abs2()));
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,-m_sff/(2.*m_zi*m_zj));
    break;
  }
}
