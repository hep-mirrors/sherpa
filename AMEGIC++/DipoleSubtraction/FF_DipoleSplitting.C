#include "AMEGIC++/DipoleSubtraction/FF_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

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
  m_a = m_yijk;

  m_ptk  = 1./(1.-m_yijk)*m_pk;
  m_ptij = m_pi+m_pj-m_yijk/(1.-m_yijk)*m_pk;

  m_zi   = (m_pi*m_ptk)/(m_ptij*m_ptk);
  m_zj   = 1.-m_zi;

  m_Q2 = (m_pi+m_pj+m_pk).Abs2();
  if (m_es==0) {
    m_kt2 = m_Q2*m_yijk*m_zi*m_zj;
  }
  else {
  m_kt2 = m_Q2*m_yijk;
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
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-(1.+m_zi);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-(1.+m_zj);
    m_av  = m_sff;
    break;
  case 3:
    m_sff = 1.;
    m_av  = m_sff - 2.0*m_zi*m_zj;
    break;
  case 4:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-2.;
    m_av  = m_sff + m_zi*m_zj;
  }
}

double FF_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj);  
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
  m_yp = 1. - 2.*(sqrt(m_mk*m_Q2) - m_mk)/(m_Q2 - m_mi - m_mj - m_mk);
  m_a = m_yijk/m_yp;
  if ((m_mij && !IsEqual(m_ptij.Abs2(),m_mij,1.0e-3)) ||
      (m_mk && !IsEqual(m_ptk.Abs2(),m_mk,1.0e-3))) {
    msg_Tracking()<<METHOD<<"(): Kinematics unstable in {\n"
		  <<"  p_i = "<<m_pi<<" "<<sqrt(dabs(m_pi.Abs2()))<<"\n"
		  <<"  p_j = "<<m_pj<<" "<<sqrt(dabs(m_pj.Abs2()))<<"\n"
		  <<"  p_k = "<<m_pk<<" "<<sqrt(dabs(m_pk.Abs2()))<<"\n}"<<std::endl;
    m_a=0.0;
  }

  m_zi   = (m_pi*m_pk)/(m_pi*m_pk+m_pj*m_pk);
  m_zj   = 1.-m_zi;

  if (m_es==0) {
    m_kt2  = 2.0*m_pi*m_pj*m_zi*m_zj-sqr(m_zi)*m_mj-sqr(m_zj)*m_mi;
  }
  else {
  m_kt2  = 2.0*m_pi*m_pj;
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

  m_vijk = Vrel(m_pi+m_pj,m_pk);
  
  m_zim  = m_zi-0.5*(1.-m_vijk);
  m_zjm  = m_zj-0.5*(1.-m_vijk);
  m_zpm  = 0.;
  if (m_ft==3 || m_ft==4) {
    m_zpm = sqr((2.*m_mi+(m_Q2-m_mi-m_mj-m_mk)*m_yijk)/(2.*(m_mi+m_mj+(m_Q2-m_mi-m_mj-m_mk)*m_yijk)))
      *(1.-sqr(m_vijk*Vrel(m_pi+m_pj,m_pi)));
  }

  m_pt1   =     m_zim*m_pi-m_zjm*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zi+m_mij/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zj+m_mij/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case 3:
    m_sff = 1.-2.*m_kappa*(m_zpm-m_mi/(m_pi+m_pj).Abs2());
    m_av  = m_sff - 2.0 * ( m_zim*m_zjm - m_zpm );
    break;
  case 4:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-(2.-m_kappa*m_zpm)/m_vijk;
    m_av  = m_sff + ( m_zim*m_zjm - m_zpm )/m_vijk;
    break;
  case 5: //gluino
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zi+m_mij/(m_pi*m_pj)); 
    m_av  = m_sff;
    break;
  case 6: //gluino
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zj+m_mij/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case 7: //squark
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(2.+m_mij/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case 8: //squark
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(2.+m_mij/(m_pi*m_pj));
    m_av  = m_sff;
  }
}

double FF_MassiveDipoleSplitting::GetValue()
{
  double h=1.0/((m_pi+m_pj).Abs2()-m_mij);  
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
  case 5:
    h*= m_sff;
    return h;
  case 6:
    h*= m_sff;   
    return h;
  case 7:
    h*= m_sff;
    return h;
  case 8:
    h*= m_sff;   
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
  default:
    return;
  }
}
