#include "AMEGIC++/DipoleSubtraction/FI_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

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
    switch (m_ftype) {
    case spt::q2qg:
      m_kt2*=m_zj;
      break;
    case spt::q2gq:
      m_kt2*=m_zi;
      break;
    case spt::g2qq:
      break;
    case spt::g2gg:
      m_kt2*=m_zi*m_zj;
      break;
    case spt::none:
      THROW(fatal_error, "Splitting type not set.");
    case spt::s2sg:
    case spt::s2gs:
    case spt::G2Gg:
    case spt::G2gG:
    case spt::V2Vg:
    case spt::V2gV:
      THROW(fatal_error, "DipoleSplitting can not handle splitting type "
          + ToString(m_ftype) + ".");
    }
  }

  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(1.-m_zi+(1.-m_xijk))-(1.+m_zi);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = 2./(1.-m_zj+(1.-m_xijk))-(1.+m_zj);
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = 1.;
    m_av  = m_sff - 2.0*m_zi*m_zj;
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
    m_av  = m_sff + m_zi*m_zj;
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    THROW(fatal_error, "DipoleSplitting can not handle splitting type "
        + ToString(m_ftype) + ".");
  }
  if (m_kt2<m_k0sqf) m_av=1.0;
}

double FI_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;
  return h*m_fac*m_sff;
}

void FI_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_ftype) {
  case spt::q2qg:
  case spt::q2gq:
    return;
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,m_sff/(4.*m_zi*m_zj));
    break;
  case spt::g2gg:
    CalcVectors(m_pt1,m_pt2,-m_sff/(2.*m_zi*m_zj));
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    THROW(fatal_error, "DipoleSplitting can not handle splitting type "
        + ToString(m_ftype) + ".");
  }
}


void FI_MassiveDipoleSplitting::SetMomenta(const Vec4D* mom )
{
  DEBUG_FUNC("m_ij^2="<<m_mij<<", m_i^2="<<m_mi<<", m_j^2="<<m_mj);
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
    switch (m_ftype) {
    case spt::q2qg:
      m_kt2*=m_zj;
      break;
    case spt::q2gq:
      m_kt2*=m_zi;
      break;
    case spt::g2qq:
      break;
    case spt::g2gg:
      m_kt2*=m_zi*m_zj;
      break;
    case spt::none:
      THROW(fatal_error, "Splitting type not set.");
    case spt::s2sg:
    case spt::s2gs:
    case spt::G2Gg:
    case spt::G2gG:
    case spt::V2Vg:
    case spt::V2gV:
      THROW(fatal_error, "DipoleSplitting can not handle splitting type "
          + ToString(m_ftype) + ".");
    }
  }

  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(2.-m_zi-m_xijk)-(1.+m_zi)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = 2./(2.-m_zj-m_xijk)-(1.+m_zj)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::g2qq: {
    m_sff = 1.;
    double Q2=2.0*m_ptij*m_pk, mui2=m_mi/Q2;
    double eps=sqrt(sqr(1.0-m_xijk-2.0*mui2)-4.0*mui2*mui2)/(1.0-m_xijk);
    m_av  = m_sff - 2.0*(0.5*(1.0+eps)-m_zi)*(m_zi-0.5*(1.0-eps));
    break;
  }
  case spt::g2gg:
    m_sff = 1./(1.-m_zi+(1.-m_xijk))+1./(1.-m_zj+(1.-m_xijk))-2.;
    m_av  = m_sff + m_zi*m_zj;
    break;
  case spt::s2sg:
    m_sff = 2./(2.-m_zi-m_xijk)-2.-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::s2gs:
    m_sff = 2./(2.-m_zj-m_xijk)-2.-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::G2Gg:
    m_sff = 2./(2.-m_zi-m_xijk)-(1.+m_zi)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::G2gG:
    m_sff = 2./(2.-m_zj-m_xijk)-(1.+m_zj)-m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::V2Vg:
    msg_Debugging()<<"Vsubmode="<<m_Vsubmode
                   <<", zi="<<m_zi<<", 1-xijk="<<1.-m_xijk
                   <<", a="<<2./(2.-m_zi-m_xijk)
                   <<", a'="<<2.*(m_pi*m_pk+m_pj*m_pk)/(m_pi*m_pj+m_pj*m_pk)
                   <<", b="<<-(1.+m_zi)<<", c="<<-m_mij/(m_pi*m_pj)<<std::endl;
    if      (m_Vsubmode==0) m_sff = 2./(2.-m_zi-m_xijk)
                                    -2.-m_mij/(m_pi*m_pj);
    else if (m_Vsubmode==1) m_sff = 2./(2.-m_zi-m_xijk)
                                    -(1.+m_zi)-m_mij/(m_pi*m_pj);
    else if (m_Vsubmode==2) m_sff = m_zi/(1.-m_zi)
                                    -m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::V2gV:
    if      (m_Vsubmode==0) m_sff = 2./(2.-m_zj-m_xijk)
                                    -2.-m_mij/(m_pi*m_pj);
    else if (m_Vsubmode==1) m_sff = 2./(2.-m_zj-m_xijk)
                                    -(1.+m_zj)-m_mij/(m_pi*m_pj);
    else if (m_Vsubmode==2) m_sff = m_zj/(1.-m_zj)
                                    -m_mij/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  }
  if (m_kt2<m_k0sqf) m_av=1.0;
}

double FI_MassiveDipoleSplitting::GetValue()
{
  double h=1.0/((m_pi+m_pj).Abs2()-m_mij)/m_xijk;
  return h*m_fac*m_sff;
}

void FI_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ftype) {
  case spt::q2qg:
  case spt::q2gq:
    return;
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,-m_sff*(m_pi+m_pj).Abs2()/(4.*m_pt1.Abs2()));
    return;
  case spt::g2gg:
    CalcVectors(m_pt1,m_pt2,-m_sff/(2.*m_zi*m_zj));
    return;
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    return;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  }
}
