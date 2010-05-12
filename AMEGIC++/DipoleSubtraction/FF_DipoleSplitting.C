#include "AMEGIC++/DipoleSubtraction/FF_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

double Lambda(double x,double y, double z)
{
  return sqr(x)+sqr(y)+sqr(z)-2.*(x*y+x*z+y*z);
}

double Vrel(Vec4D p1,Vec4D& p2)
{
  double m21=(p1+p2).Abs2();
  double m1=p1.Abs2();
  double m2=p2.Abs2();
  return sqrt(Lambda(m21,m1,m2))/(m21-m1-m2);
}

void FF_DipoleSplitting::SetMomenta(const Vec4D* mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_yijk = m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);

  m_ptk  = 1./(1.-m_yijk)*m_pk;
  m_ptij = m_pi+m_pj-m_yijk/(1.-m_yijk)*m_pk;

  m_zi   = (m_pi*m_ptk)/(m_ptij*m_ptk);
  m_zj   = 1.-m_zi;

//   m_pt1   =     m_zi*m_pi;
//   m_pt2   = -1.*m_zj*m_pj;
  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-(1.+m_zi);
    break;
  case 2:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-(1.+m_zj);
    break;
  case 3:
    m_sff = 1.;
    break;
  case 4:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-2.;
  }
  if (1) {
    double pf = m_yijk/(m_alpha*(1.-0.5*m_alpha));
    switch (m_ft) {
    case 1:
    case 2:
      m_sff -= pf*Vcijf(m_ft);
      break;
    case 3:
      m_sff -= pf*Vcijf(m_ft);
    break;
    case 4:
      m_sff -= pf*Vcijf(m_ft);
    }
  }
}

double FF_DipoleSplitting::GetF()
{
   if (m_yijk>m_alpha) return 0.;

   if (m_yijk<=m_amin) {
      return nan;
   }

  double h=SPFac()/(2.*m_pi*m_pj);  
  switch (m_ft) {
  case 1:
    h*= m_sff;
    return h;
  case 2:
    h*= m_sff;   
    return h;
  case 3:
    return h*m_sff*CSC.TR/CSC.CA;
  case 4:
    h*=2.*m_sff;
    return h;
  }
  return 0.;
}

void FF_DipoleSplitting::CalcDiPolarizations()
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


void FF_MassiveDipoleSplitting::SetMomenta(const Vec4D* mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];
  
  m_Q = m_pi+m_pj+m_pk;
  m_Q2 = m_Q.Abs2();

  m_ptk  = sqrt(Lambda(m_Q2,m_mij,m_mk)/Lambda(m_Q2,(m_pi+m_pj).Abs2(),m_mk))*(m_pk-(m_Q*m_pk)/m_Q2*m_Q)+
    (m_Q2+m_mk-m_mij)/(2.*m_Q2)*m_Q;
  m_ptij = m_Q-m_ptk;

  m_yijk = m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);

  m_zi   = (m_pi*m_pk)/(m_pi*m_pk+m_pj*m_pk);
  m_zj   = 1.-m_zi;

  m_vijk = Vrel(m_pi+m_pj,m_pk);
  
  m_zim  = m_zi-0.5*(1.-m_vijk);
  m_zjm  = m_zj-0.5*(1.-m_vijk);
  m_zpm  = 0.;
  if (m_ft>2) {
    m_zpm = (2.*m_mi+(m_Q2-m_mi-m_mj-m_mk)*m_yijk)/(2.*(m_mi+m_mj+(m_Q2-m_mi-m_mj-m_mk)*m_yijk))
      *(1.-sqr(m_vijk*Vrel(m_pi+m_pj,m_pi)));
  }

  m_pt1   =     m_zim*m_pi-m_zjm*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zi+m_mij/(m_pi*m_pj));
    break;
  case 2:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zj+m_mij/(m_pi*m_pj));
    break;
  case 3:
    m_sff = 1.-2.*m_kappa*(m_zpm-m_mi/(m_pi+m_pj).Abs2());
    break;
  case 4:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-(2.-m_kappa*m_zpm)/m_vijk;
  }
}

double FF_MassiveDipoleSplitting::GetF()
{
   if (m_yijk>m_alpha) return 0.;
   if (m_yijk<=m_amin) {
      return nan;
   }

  double h=SPFac()/((m_pi+m_pj).Abs2()-m_mij);  
  switch (m_ft) {
  case 1:
    h*= m_sff;
    return h;
  case 2:
    h*= m_sff;   
    return h;
  case 3:
    return h*m_sff*CSC.TR/CSC.CA;
  case 4:
    h*=2.*m_sff;
    return h;
  }
  return 0.;
}

void FF_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,-m_sff*(m_pi+m_pj).Abs2()/(4.*m_pt1.Abs2()));
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_vijk/(2.*m_zim*m_zjm));
    break;
  }
}
