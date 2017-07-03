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
  if (m_es==0) {
    m_kt2 = -m_Q2*(1.-m_xijk)/m_xijk*m_uj;
  }
  else {
    m_kt2 = -m_Q2/m_xijk*m_uj;
    switch (m_ftype) {
    case spt::q2qg:
      m_kt2*=(1.-m_xijk);
      break;
    case spt::g2gg:
      m_kt2*=(1.-m_xijk);
      break;
    }
  }

//   m_pt1  =    m_pj/m_uj;
//   m_pt2  =-1.*m_pk/m_uk;
  m_pt1  =    m_pj/m_uj-m_pk/m_uk;
  m_pt2  =    m_ptij;

  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(1.-m_xijk+m_uj)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = (1.-2.*m_xijk*(1.-m_xijk));
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk;
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_xijk+m_uj)-1.+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk;
  }
  if (m_kt2<m_k0sqi) m_av=1.0;
}

double IF_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;
  return h*m_fac*m_sff;
}

void IF_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_ftype) {
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xijk/(1.-m_xijk)/4.);
    break;
  case spt::g2gg:
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

  if (m_es==0) {
    m_kt2 = 2.0*m_pj*m_pk*m_uj;
  }
  else {
    m_kt2 = 2.0*m_pj*m_pk*m_uj/(1.-m_xijk);
    switch (m_ftype) {
    case spt::q2qg:
      m_kt2*=(1.-m_xijk);
      break;
    case spt::g2gg:
      m_kt2*=(1.-m_xijk);
      break;
    }
  }

  m_pt1  =    m_pj/m_uj-m_pk/m_uk;
  m_pt2  =    m_ptij;
  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(1.-m_xijk+m_uj)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = (1.-2.*m_xijk*(1.-m_xijk));
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk
                  - m_pk.Abs2()/(m_ptk*m_ptij)*m_uj/m_uk;
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_xijk+m_uj)-1.+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk
                  - m_pk.Abs2()/(2.0*m_ptk*m_ptij)*m_uj/m_uk;
  }
  if (m_kt2<m_k0sqi) m_av=1.0;
}

double IF_MassiveDipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;
  return h*m_fac*m_sff;
}

void IF_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ftype) {
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,m_sff*m_xijk/(1.-m_xijk)/2.*(m_pk*m_pj)/(m_uj*m_uk)/m_pt1.Abs2());
    break;
  case spt::g2gg:
    CalcVectors(m_pt1,m_pt2,m_sff*m_xijk/(1.-m_xijk)*(m_pk*m_pj)/(m_uj*m_uk)/m_pt1.Abs2());
    break;
  }
}
