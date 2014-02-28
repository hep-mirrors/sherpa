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
  m_a = 1.0-m_xijk;

  m_ptk  = m_xijk*m_pk;
  m_ptij = m_pi+m_pj-(1.-m_xijk)*m_pk;

  m_zi   = (m_pi*m_ptk)/(m_ptij*m_ptk);
  m_zj   = 1.-m_zi;

  m_Q2 = (m_pi+m_pj-m_pk).Abs2();
  if (m_es==0) {
    m_kt2 = -m_Q2*(1.-m_xijk)/m_xijk*m_zi*m_zj;
  }
  else {
  m_kt2 = -m_Q2*(1.-m_xijk)/m_xijk;
  switch (m_ft) {
  case 1:
    m_kt2*=m_zj;
    break;
  case 2:
    m_kt2*=m_zi;
    break;
  case 4:
    m_kt2*=m_zi*m_zj;
    break;
  }
  }

  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_zi+(1.-m_xijk))-(1.+m_zi);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 2./(1.-m_zj+(1.-m_xijk))-(1.+m_zj);
    m_av  = m_sff;
    break;
  case 3:
    m_sff = 1.;
    m_av  = m_sff - 2.0*m_zi*m_zj;
    break;
  case 4:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
    m_av  = m_sff + m_zi*m_zj;
  }
  if (m_kt2<m_k0sqf) m_av=1.0;
}

double FI_DipoleSplitting::GetValue()
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
  m_a = 1.0-m_xijk;
  
  m_ptk  = m_xijk*m_pk;
  m_ptij = m_pi+m_pj-(1.-m_xijk)*m_pk;

  m_zi   = (m_pi*m_pk)/(m_pj*m_pk+m_pk*m_pi);
  m_zj   = 1.-m_zi;

  if (m_es==0) {
    m_kt2 = 2.0*m_pi*m_pj*m_zi*m_zj-sqr(m_zi)*m_mj-sqr(m_zj)*m_mi;
  }
  else {
  m_kt2 = 2.0*m_pi*m_pj;
  switch (m_ft) {
  case 1:
    m_kt2*=m_zj;
    break;
  case 2:
    m_kt2*=m_zi;
    break;
  case 4:
    m_kt2*=m_zi*m_zj;
    break;
  }
  }

//   m_pt1   =     m_zi*m_pi;
//   m_pt2   = -1.*m_zj*m_pj;
  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(2.-m_zi-m_xijk)-(1.+m_zi)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 2./(2.-m_zj-m_xijk)-(1.+m_zj)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case 3: {
    m_sff = 1.;
    double Q2=2.0*m_ptij*m_pk, mui2=m_mi/Q2;
    double eps=sqrt(sqr(1.0-m_xijk-2.0*mui2)-4.0*mui2*mui2)/(1.0-m_xijk);
    m_av  = m_sff - 2.0*(0.5*(1.0+eps)-m_zi)*(m_zi-0.5*(1.0-eps));
    break;
  }
  case 4:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
    m_av  = m_sff + m_zi*m_zj;
    break;
  case 5:
    m_sff = 2./(2.-m_zi-m_xijk)-(1.+m_zi)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case 6:
    m_sff = 2./(2.-m_zj-m_xijk)-(1.+m_zj)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case 7:
    m_sff = 2./(2.-m_zi-m_xijk)-2.-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case 8:
    m_sff = 2./(2.-m_zj-m_xijk)-2.-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  }
  if (m_kt2<m_k0sqf) m_av=1.0;
}

double FI_MassiveDipoleSplitting::GetValue()
{
  double h=1.0/((m_pi+m_pj).Abs2()-m_mij)/m_xijk;
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
  case 5:
    h*=m_sff;
    return h;
  case 6:
    h*=m_sff;   
    return h;
  case 7:
    h*=m_sff;
    return h;
  case 8:
    h*=m_sff;   
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
