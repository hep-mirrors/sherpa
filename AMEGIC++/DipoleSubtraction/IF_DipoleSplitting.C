#include "AMEGIC++/DipoleSubtraction/IF_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

void IF_DipoleSplitting::SetMomenta(const Vec4D *mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = 1.-m_pj*m_pk/(m_pi*m_pj+m_pk*m_pi);
  m_ptk  = m_pk+m_pj-(1.-m_xijk)*m_pi;
  m_ptij = m_xijk*m_pi;

  m_uj   = (m_pi*m_pj)/(m_pi*m_pj+m_pk*m_pi);
  m_uk   = 1.-m_uj;
  m_a = m_uj;

  m_Q2 = (-m_pi+m_pj+m_pk).Abs2();
  m_kt2 = -m_Q2*(1.-m_xijk)/m_xijk*m_uj;

//   m_pt1  =    m_pj/m_uj;
//   m_pt2  =-1.*m_pk/m_uk;
  m_pt1  =    m_pj/m_uj-m_pk/m_uk;
  m_pt2  =    m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_xijk+m_uj)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = (1.-2.*m_xijk*(1.-m_xijk));
    m_av  = m_sff;
    break;
  case 3:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk;
    break;
  case 4:
    m_sff = 1./(1.-m_xijk+m_uj)-1.+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk;
  }
}

double IF_DipoleSplitting::GetF()
{
  if (Reject(m_uj)) return 0.;
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

void IF_DipoleSplitting::CalcDiPolarizations()
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

void IF_MassiveDipoleSplitting::SetMomenta(const Vec4D *mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = 1.-m_pj*m_pk/(m_pi*m_pj+m_pk*m_pi);
  m_ptk  = m_pk+m_pj-(1.-m_xijk)*m_pi;
  m_ptij = m_xijk*m_pi;

  m_uj   = (m_pi*m_pj)/(m_pi*m_pj+m_pk*m_pi);
  m_uk   = 1.-m_uj;
  m_a = m_uj;

  m_kt2 = 2.0*m_pj*m_pk*m_uj;

//   m_pt1  =    m_pj/m_uj;
//   m_pt2  =-1.*m_pk/m_uk;
  m_pt1  =    m_pj/m_uj-m_pk/m_uk;
  m_pt2  =    m_ptij;
  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_xijk+m_uj)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = (1.-2.*m_xijk*(1.-m_xijk));
    m_av  = m_sff;
    break;
  case 3:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk - m_pk.Abs2()/(m_ptk*m_ptij)*m_uj/m_uk;
    break;
  case 4:
    m_sff = 1./(1.-m_xijk+m_uj)-1.+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk - m_pk.Abs2()/(2.0*m_ptk*m_ptij)*m_uj/m_uk;
  }
}

double IF_MassiveDipoleSplitting::GetF()
{
  if (Reject(m_uj)) return 0.;
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

void IF_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,m_sff*m_xijk/(1.-m_xijk)/2.*(m_pk*m_pj)/(m_uj*m_uk)/m_pt1.Abs2());
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,m_sff*m_xijk/(1.-m_xijk)*(m_pk*m_pj)/(m_uj*m_uk)/m_pt1.Abs2());
    break;
  }
}
